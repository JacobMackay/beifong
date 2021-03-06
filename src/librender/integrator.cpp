#include <thread>
#include <mutex>

#include <enoki/morton.h>
#include <mitsuba/core/profiler.h>
#include <mitsuba/core/progress.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/film.h>
#include <mitsuba/render/adc.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/receiver.h>
#include <mitsuba/render/spiral.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <mutex>

NAMESPACE_BEGIN(mitsuba)

// -----------------------------------------------------------------------------

MTS_VARIANT SamplingIntegrator<Float, Spectrum>::
  SamplingIntegrator(const Properties &props): Base(props) {
      m_block_size = (uint32_t) props.size_("block_size", 0);
      uint32_t block_size = math::round_to_power_of_two(m_block_size);
      if (m_block_size > 0 && block_size != m_block_size) {
          Log(Warn,
              "Setting block size from %i to next higher power of two: %i",
              m_block_size, block_size);
          m_block_size = block_size;
      }

      m_samples_per_pass = (uint32_t) props.size_("samples_per_pass",
                                                    (size_t) -1);
      m_timeout = props.float_("timeout", -1.f);

      /// Disable direct visibility of emitters if needed
      m_hide_emitters = props.bool_("hide_emitters", false);
}

MTS_VARIANT SamplingIntegrator<Float, Spectrum>::
  ~SamplingIntegrator() { }

MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  cancel() {
      m_stop = true;
}

MTS_VARIANT std::vector<std::string> SamplingIntegrator<Float, Spectrum>::
  aov_names() const {
      return { };
}

MTS_VARIANT bool SamplingIntegrator<Float, Spectrum>::
  render(Scene *scene, Sensor *sensor) {
      ScopedPhase sp(ProfilerPhase::Render);
      m_stop = false;

      ref<Film> film = sensor->film();
      ScalarVector2i film_size = film->crop_size();

      size_t total_spp = sensor->sampler()->sample_count();
      size_t samples_per_pass = (m_samples_per_pass == (size_t) -1)
                               ? total_spp
                               : std::min(
                                   (size_t) m_samples_per_pass, total_spp);
      if ((total_spp % samples_per_pass) != 0) {
          Throw("sample_count (%d) must be a multiple of samples_per_pass (%d)."
          , total_spp, samples_per_pass);
      }

      size_t n_passes = (total_spp + samples_per_pass - 1) / samples_per_pass;

      std::vector<std::string> channels = aov_names();
      bool has_aovs = !channels.empty();

      // Insert default channels and set up the film
      for (size_t i = 0; i < 5; ++i) {
          channels.insert(channels.begin() + i, std::string(1, "XYZAW"[i]));
          film->prepare(channels);
      }

      m_render_timer.reset();
      if constexpr (!is_cuda_array_v<Float>) {
          /// Render on the CPU using a spiral pattern
          size_t n_threads = __global_thread_count;
          Log(Info, "Starting render job (%ix%i, %i sample%s,%s %i thread%s)",
            film_size.x(), film_size.y(),
            total_spp, total_spp == 1 ? "" : "s",
            n_passes > 1 ? tfm::format(" %d passes,", n_passes) : "",
            n_threads, n_threads == 1 ? "" : "s");

          if (m_timeout > 0.f) {
              Log(Info, "Timeout specified: %.2f seconds.", m_timeout);
          }

          // Find a good block size to use for splitting up the total workload.
          if (m_block_size == 0) {
              uint32_t block_size = MTS_BLOCK_SIZE;
              while (true) {
                  if (block_size == 1 ||
                      hprod((film_size + block_size - 1) / block_size)
                      >= n_threads) {
                          break;
                  }
                  block_size /= 2;
              }
              m_block_size = block_size;
          }

          Spiral spiral(film, m_block_size, n_passes);

          ThreadEnvironment env;
          ref<ProgressReporter> progress = new ProgressReporter("Rendering");
          std::mutex mutex;

          // Total number of blocks to be handled, including multiple passes.
          size_t total_blocks = spiral.block_count() * n_passes,
            blocks_done = 0;

          tbb::parallel_for(
            tbb::blocked_range<size_t>(0, total_blocks, 1),
            [&](const tbb::blocked_range<size_t> &range) {
                ScopedSetThreadEnvironment set_env(env);
                ref<Sampler> sampler = sensor->sampler()->clone();
                // ref<ImageBlock> block =
                //     new ImageBlock(m_block_size, channels.size(),
                //                     film->reconstruction_filter(), !has_aovs);
                ref<ImageBlock> block =
                    new ImageBlock(m_block_size, channels.size(),
                                    film->reconstruction_filter(), false);
                scoped_flush_denormals flush_denormals(true);
                std::unique_ptr<Float[]> aovs(new Float[channels.size()]);

                // For each block
                for (auto i = range.begin(); i != range.end() && !should_stop();
                    ++i) {
                    auto [offset, size, block_id] = spiral.next_block();
                    Assert(hprod(size) != 0);
                    block->set_size(size);
                    block->set_offset(offset);

                    render_block(scene, sensor, sampler, block,
                             aovs.get(), samples_per_pass, block_id);

                    film->put(block);

                    /* Critical section: update progress bar */ {
                        std::lock_guard<std::mutex> lock(mutex);
                        blocks_done++;
                        progress->update(blocks_done /
                            (ScalarFloat) total_blocks);
                    }
                }
          });
      } else {
          Log(Info, "Start rendering...");

          ref<Sampler> sampler = sensor->sampler();
          sampler->set_samples_per_wavefront((uint32_t) samples_per_pass);

          ScalarFloat diff_scale_factor =
            rsqrt((ScalarFloat) sampler->sample_count());
          ScalarUInt32 wavefront_size =
            hprod(film_size) * (uint32_t) samples_per_pass;
          if (sampler->wavefront_size() != wavefront_size) {
            sampler->seed(0, wavefront_size);
          }

          UInt32 idx = arange<UInt32>(wavefront_size);
          if (samples_per_pass != 1) {
              idx /= (uint32_t) samples_per_pass;
          }

          // ref<ImageBlock> block = new ImageBlock(film_size, channels.size(),
          //                                      film->reconstruction_filter(),
          //                                      !has_aovs);
          ref<ImageBlock> block = new ImageBlock(film_size, channels.size(),
                                               film->reconstruction_filter(),
                                               false);
          block->clear();
          Vector2f pos = Vector2f(Float(idx % uint32_t(film_size[0])),
                                Float(idx / uint32_t(film_size[0])));
          std::vector<Float> aovs(channels.size());

          for (size_t i = 0; i < n_passes; i++) {
              render_sample(scene, sensor, sampler, block, aovs.data(),
                          pos, diff_scale_factor);
          }

          film->put(block);
      }

      if (!m_stop) {
        Log(Info, "Rendering finished. (took %s)",
            util::time_string(m_render_timer.value(), true));
      }

      return !m_stop;
}

MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  render_block(const Scene *scene, const Sensor *sensor, Sampler *sampler,
      ImageBlock *block, Float *aovs, size_t sample_count_, size_t block_id)
      const {
    block->clear();
    uint32_t pixel_count  = (uint32_t)(m_block_size * m_block_size),
             sample_count = (uint32_t)(sample_count_ == (size_t) -1
                                           ? sampler->sample_count()
                                           : sample_count_);

    ScalarFloat diff_scale_factor =
        rsqrt((ScalarFloat) sampler->sample_count());

    if constexpr (!is_array_v<Float>) {
        for (uint32_t i = 0; i < pixel_count && !should_stop(); ++i) {
            sampler->seed(block_id * pixel_count + i);

            ScalarPoint2u pos = enoki::morton_decode<ScalarPoint2u>(i);
            if (any(pos >= block->size()))
                continue;

            pos += block->offset();
            for (uint32_t j = 0; j < sample_count && !should_stop(); ++j) {
                render_sample(scene, sensor, sampler, block, aovs,
                              pos, diff_scale_factor);
            }
        }
    } else if constexpr (is_array_v<Float> && !is_cuda_array_v<Float>) {
        // Ensure that the sample generation is fully deterministic
        sampler->seed(block_id);

        for (auto [index, active] : range<UInt32>(pixel_count * sample_count)) {
            if (should_stop()) {
                break;
            }
            Point2u pos = enoki::morton_decode<Point2u>
                (index / UInt32(sample_count));
            active &= !any(pos >= block->size());
            pos += block->offset();
            render_sample(scene, sensor, sampler, block, aovs,
                pos, diff_scale_factor, active);
        }
    } else {
        ENOKI_MARK_USED(scene);
        ENOKI_MARK_USED(sensor);
        ENOKI_MARK_USED(aovs);
        ENOKI_MARK_USED(diff_scale_factor);
        ENOKI_MARK_USED(pixel_count);
        ENOKI_MARK_USED(sample_count);
        Throw("Not implemented for CUDA arrays.");
    }
}

MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  render_sample(const Scene *scene, const Sensor *sensor, Sampler *sampler,
                    ImageBlock *block, Float *aovs, const Vector2f &pos,
                    ScalarFloat diff_scale_factor, Mask active) const {
        Vector2f position_sample = pos + sampler->next_2d(active);

        Point2f aperture_sample(.5f);
        if (sensor->needs_aperture_sample())
            aperture_sample = sampler->next_2d(active);

        Float time = sensor->shutter_open();
        if (sensor->shutter_open_time() > 0.f)
            time += sampler->next_1d(active) * sensor->shutter_open_time();

        Float wavelength_sample = sampler->next_1d(active);

        Vector2f adjusted_position =
            (position_sample - sensor->film()->crop_offset()) /
            sensor->film()->crop_size();

        auto [ray, ray_weight] = sensor->sample_ray_differential(
            time, wavelength_sample, adjusted_position, aperture_sample);

        ray.scale_differential(diff_scale_factor);

        const Medium *medium = sensor->medium();
        std::pair<Spectrum, Mask> result = sample(scene, sampler, ray, medium, aovs + 5, active);
        // std::pair<Spectrum, Mask> result = sample(scene, sampler, &ray, medium, aovs + 5, active);
        result.first = ray_weight * result.first;

        UnpolarizedSpectrum spec_u = depolarize(result.first);

        Color3f xyz;
        if constexpr (is_monochromatic_v<Spectrum>) {
            xyz = spec_u.x();
        } else if constexpr (is_rgb_v<Spectrum>) {
            xyz = srgb_to_xyz(spec_u, active);
        } else {
            static_assert(is_spectral_v<Spectrum>);
            xyz = spectrum_to_xyz(spec_u, ray.wavelengths, active);
        }

        aovs[0] = xyz.x();
        aovs[1] = xyz.y();
        aovs[2] = xyz.z();
        aovs[3] = select(result.second, Float(1.f), Float(0.f));
        aovs[4] = 1.f;

        block->put(position_sample, aovs, active);

        sampler->advance();
}

// ============================================================================
// Radar Receive Section
// ============================================================================
MTS_VARIANT bool SamplingIntegrator<Float, Spectrum>::
  receive(Scene *scene, Receiver *receiver) {
      ScopedPhase sp(ProfilerPhase::Receive);
      m_stop = false;

      // Changing the film and receive characteristics is a deep fundamental
      // change and will require a lot.

      // How can I do it differently? Fuck ton of aovs.

      // What is a 'film'? A film is an object which is a carrier for data
      // storage, typically in 2d corresponding to planar coordinates.

      // I will have a matrix...but this is a conflicting name...how about tessy
      // as short for tesseract? Or just a tess.

      // Maybe a 6 bucket is better? Start with 4

      // What is a sensor? A sensor is an object which occupies some space, and
      // through interactions collects signals.

      ref<ADC> adc = receiver->adc();
      ScalarVector2i adc_rd_size = adc->window_size();
      ScalarVector2i adc_rd_offset = adc->window_offset();
      ScalarVector2i adc_px_size = (1, 1);
      ScalarVector2i adc_px_offset = (0, 0);

      // I'll either need to make a different sensor or modify it. I guess This
      // is a deep change.
      // ref<Tess> tess = gamjigi->tess();
      // ScalarVector4i tess_size = tess->crop_size();

      size_t total_spp = receiver->sampler()->sample_count();
      size_t samples_per_pass = (m_samples_per_pass == (size_t) -1)
                               ? total_spp
                               : std::min(
                                   (size_t) m_samples_per_pass, total_spp);
      if ((total_spp % samples_per_pass) != 0) {
          Throw("sample_count (%d) must be a multiple of samples_per_pass (%d)."
          , total_spp, samples_per_pass);
      }

      size_t n_passes = (total_spp + samples_per_pass - 1) / samples_per_pass;

      std::vector<std::string> channels = aov_names();
      bool has_aovs = !channels.empty();

      // Insert default channels and set up the film
      // for (size_t i = 0; i < 5; ++i) {
      //     channels.insert(channels.begin() + i, std::string(1, "XYZAW"[i]));
      //     film->prepare(channels);
      // }
      for (size_t i = 0; i < 3; ++i) {
          channels.insert(channels.begin() + i, std::string(1, "YAW"[i]));
          adc->prepare(channels);
      }

      // Actually...it makes more sense to be a time-frequency method...maybe..
      // Maybe that only makes sense at baseband.
      // For now, let's not do the spiral method for range-doppler
      // We would have to do something like:

      // Can consider blocks now.

      m_render_timer.reset();
      // if constexpr (!is_cuda_array_v<Float>) {
      // This doesn't make sense for radar atm...
      if constexpr (false) {
          /// Render on the CPU using a spiral pattern

          // std::cout << "Hello" << std::endl;

          size_t n_threads = __global_thread_count;
          Log(Info, "Starting signal render job (%ix%i rd, %ix%i xy,"
            " %i sample%s,%s %i thread%s)",
            adc_rd_size.x(), adc_rd_size.y(),
            adc_px_size.x(), adc_px_size.y(),
            total_spp, total_spp == 1 ? "" : "s",
            n_passes > 1 ? tfm::format(" %d passes,", n_passes) : "",
            n_threads, n_threads == 1 ? "" : "s");

          if (m_timeout > 0.f) {
              Log(Info, "Timeout specified: %.2f seconds.", m_timeout);
          }

          // Blocking and spiralling doesn't make sense for traditional
          // renderingxtemporals. Unless the render pipeline explicitly only
          // renders range/doppler for a given interval.
          // With that in mind, I need to do something different...for now.
          // This is partly bc the pixels are independent in xy
          // This could be really good for ensuring incoherence between distant
          // blocks.

          // Find a good block size to use for splitting up the total workload.
          if (m_block_size == 0) {
              uint32_t block_size = MTS_BLOCK_SIZE;
              while (true) {
                  // Our film size pixel-wise should now be 1
                  // Old:
                  // if (block_size == 1 ||
                  //     hprod((film_size + block_size - 1) / block_size)
                  //     >= n_threads) {
                  //         break;
                  // }
                  // New:
                  if (block_size == 1 ||
                      hprod((adc_px_size + block_size - 1) / block_size)
                      >= n_threads) {
                          break;
                  }
                  block_size /= 2;
              }
              m_block_size = block_size;
          }

          // Here's another problem, our block size is now tied to pixels
          std::cout << m_block_size << std::endl;

          // I'm worried about this in that it may use film size to apply
          // offsets
          // This spiral would be spiralling over the range and doppler
          // Old:
          // Spiral spiral(film, m_block_size, n_passes);
          // New:
          Spiral spiral(adc_px_size, adc_px_offset, m_block_size, n_passes);

          ThreadEnvironment env;
          ref<ProgressReporter> progress =
            new ProgressReporter("Rendering Signal");
          std::mutex mutex;

          // Total number of blocks to be handled, including multiple passes.
          size_t total_blocks = spiral.block_count() * n_passes,
            blocks_done = 0;

          tbb::parallel_for(
            tbb::blocked_range<size_t>(0, total_blocks, 1),
            [&](const tbb::blocked_range<size_t> &range) {
                ScopedSetThreadEnvironment set_env(env);
                ref<Sampler> sampler = receiver->sampler()->clone();
                ref<SignalBlock> block =
                    // new SignalBlock(m_block_size, channels.size(),
                    //                 adc->reconstruction_filter(), !has_aovs);
                    new SignalBlock(m_block_size, channels.size(),
                                    adc->reconstruction_filter(), false);
                scoped_flush_denormals flush_denormals(true);
                std::unique_ptr<Float[]> aovs(new Float[channels.size()]);

                // For each block
                for (auto i = range.begin(); i != range.end() && !should_stop();
                    ++i) {
                    auto [offset, size, block_id] = spiral.next_block();
                    Assert(hprod(size) != 0);
                    block->set_size(size);
                    block->set_offset(offset);

                    receive_block(scene, receiver, sampler, block,
                             aovs.get(), samples_per_pass, block_id);
                    // std::cout << "Hello2" << std::endl;
                    adc->put(block);

                    /* Critical section: update progress bar */ {
                        std::lock_guard<std::mutex> lock(mutex);
                        blocks_done++;
                        progress->update(blocks_done /
                            (ScalarFloat) total_blocks);
                    }
                }
          });
      } else if constexpr (!is_cuda_array_v<Float>) {
          /// Render on the CPU using no pattern

          // size_t n_threads = __global_thread_count;
          // Log(Info, "Starting render job (%ix%i rd, %ix%i xy,"
          //   " %i sample%s,%s %i thread%s)",
          //   film_rd_size.x(), film_rd_size.y(),
          //   film_px_size.x(), film_px_size.y(),
          //   total_spp, total_spp == 1 ? "" : "s",
          //   n_passes > 1 ? tfm::format(" %d passes,", n_passes) : "",
          //   n_threads, n_threads == 1 ? "" : "s");
          //
          // if (m_timeout > 0.f) {
          //     Log(Info, "Timeout specified: %.2f seconds.", m_timeout);
          // }
          //
          // // Blocking and spiralling doesn't make sense for traditional
          // // renderingxtemporals. Unless the render pipeline explicitly only
          // // renders range/doppler for a given interval.
          // // With that in mind, I need to do something different...for now.
          // // This is partly bc the pixels are independent in xy
          // // This could be really good for ensuring incoherence between distant
          // // blocks.
          //
          // // Find a good block size to use for splitting up the total workload.
          // if (m_block_size == 0) {
          //     uint32_t block_size = MTS_BLOCK_SIZE;
          //     while (true) {
          //         // Our film size pixel-wise should now be 1
          //         // Old:
          //         // if (block_size == 1 ||
          //         //     hprod((film_size + block_size - 1) / block_size)
          //         //     >= n_threads) {
          //         //         break;
          //         // }
          //         // New:
          //         if (block_size == 1 ||
          //             hprod((film_px_size + block_size - 1) / block_size)
          //             >= n_threads) {
          //                 break;
          //         }
          //         block_size /= 2;
          //     }
          //     m_block_size = block_size;
          // }
          //
          //
          // // I'm worried about this in that it may use film size to apply
          // // offsets
          // // This spiral would be spiralling over the range and doppler
          // // Old:
          // // Spiral spiral(film, m_block_size, n_passes);
          // // New:
          // Spiral spiral(film_px_size, film_px_offset, m_block_size, n_passes);
          //
          // ThreadEnvironment env;
          // ref<ProgressReporter> progress =
          //   new ProgressReporter("Rendering Signal");
          // std::mutex mutex;
          //
          // // Total number of blocks to be handled, including multiple passes.
          // size_t total_blocks = spiral.block_count() * n_passes,
          //   blocks_done = 0;
          //
          // tbb::parallel_for(
          //   tbb::blocked_range<size_t>(0, total_blocks, 1),
          //   [&](const tbb::blocked_range<size_t> &range) {
          //       ScopedSetThreadEnvironment set_env(env);
          //       ref<Sampler> sampler = sensor->sampler()->clone();
          //       ref<ImageBlock> block =
          //           new ImageBlock(m_block_size, channels.size(),
          //                           film->reconstruction_filter(), !has_aovs);
          //       scoped_flush_denormals flush_denormals(true);
          //       std::unique_ptr<Float[]> aovs(new Float[channels.size()]);
          //
          //       // For each block
          //       for (auto i = range.begin(); i != range.end() && !should_stop();
          //           ++i) {
          //           auto [offset, size, block_id] = spiral.next_block();
          //           Assert(hprod(size) != 0);
          //           block->set_size(size);
          //           block->set_offset(offset);
          //
          //           receive_block(scene, sensor, sampler, block,
          //                    aovs.get(), samples_per_pass, block_id);
          //           std::cout << "Hello2" << std::endl;
          //           film->put(block);
          //
          //           /* Critical section: update progress bar */ {
          //               std::lock_guard<std::mutex> lock(mutex);
          //               blocks_done++;
          //               progress->update(blocks_done /
          //                   (ScalarFloat) total_blocks);
          //           }
          //       }
          // });

          // std::cout << "Hello_rd" << std::endl;

          // // Original:
          // Log(Info, "Start signal rendering...");
          // Supress Logs:
          Log(Debug, "Start signal rendering...");

          ref<Sampler> sampler = receiver->sampler();
          // sampler->set_samples_per_wavefront((uint32_t) samples_per_pass);

          ScalarFloat diff_scale_factor =
            rsqrt((ScalarFloat) sampler->sample_count());

          // This needs to change. Film_size now means range/doppler bins.
          // Old:
          // ScalarUInt32 wavefront_size =
          //   hprod(film_size) * (uint32_t) samples_per_pass;
          // if (sampler->wavefront_size() != wavefront_size) {
          //   sampler->seed(0, wavefront_size);
          // }

          // New:
          // ScalarUInt32 wavefront_size = (uint32_t) samples_per_pass;
          // if (sampler->wavefront_size() != wavefront_size) {
          //   sampler->seed(0, wavefront_size);
          // }

          // idx from 0:nsamples
          // UInt32 idx = arange<UInt32>(wavefront_size);
          // if (samples_per_pass != 1) {
          //     // idx from 0:1, how does this work as uint32?
          //     idx /= (uint32_t) samples_per_pass;
          // }
          // std::cout<<idx<<std::endl;

          // Old:
          // ref<ImageBlock> block = new ImageBlock(film_size, channels.size(),
          //                                      film->reconstruction_filter(),
          //                                      !has_aovs);
            // New:
          // ref<SignalBlock> block = new SignalBlock(adc_rd_size, channels.size(),
          //                                      adc->reconstruction_filter(),
          //                                      !has_aovs);
          ref<SignalBlock> block = new SignalBlock(adc_rd_size, channels.size(),
                                               adc->reconstruction_filter(),
                                               false);
           block->set_offset(receiver->adc()->window_offset());
           block->set_size(receiver->adc()->window_size());
          block->clear();

          // Maybe add noise

          // Modify something depending on film type.
          // Vector2f pos = Vector2f(Float(idx % uint32_t(film_size[0])),
          //                       Float(idx / uint32_t(film_size[0])));

          // When considering range doppler, there should be only 1
          // 'positionpixel'
          // We will now move across the 'pixel'
          // Pos should be of size (2,numsampes*numpixels)

          // Throws a malloc error if film size isn't 1x1
          // Old:
          // Vector2f pos = Vector2f(Float(idx % uint32_t(1)),
          //                       Float(idx / uint32_t(1)));
          // New:
          // Vector2f pos = Vector2f(1.f, 1.f);
          // Vector2f pos = Vector2f(0.5f, 0.5f);
          Vector2f pos = Vector2f(0.f, 0.f);
          std::vector<Float> aovs(channels.size());

          // Could make this parallel?
          // Old:
          // for (size_t i = 0; i < n_passes; i++) {
          //     receive_sample(scene, sensor, sampler, block, aovs.data(),
          //                 pos, diff_scale_factor);
          // }
          // New:
          for (size_t i = 0; i < total_spp; i++) {
              // std::cout << i << std::endl;
              receive_sample(scene, receiver, sampler, block, aovs.data(),
                          pos, diff_scale_factor);
          }


          adc->put(block);













      } else {
          std::cout << "Hello_c" << std::endl;

          // // Original:
          // Log(Info, "Start signal rendering...");
          // Supress Logs:
          Log(Debug, "Start signal rendering...");


          ref<Sampler> sampler = receiver->sampler();
          sampler->set_samples_per_wavefront((uint32_t) samples_per_pass);

          ScalarFloat diff_scale_factor =
            rsqrt((ScalarFloat) sampler->sample_count());

          // This needs to change. Film_size now means range/doppler bins.
          // Old:
          // ScalarUInt32 wavefront_size =
          //   hprod(film_size) * (uint32_t) samples_per_pass;
          // if (sampler->wavefront_size() != wavefront_size) {
          //   sampler->seed(0, wavefront_size);
          // }

          // New:
          ScalarUInt32 wavefront_size = (uint32_t) samples_per_pass;
          if (sampler->wavefront_size() != wavefront_size) {
            sampler->seed(0, wavefront_size);
          }

          // idx from 0:nsamples
          UInt32 idx = arange<UInt32>(wavefront_size);
          if (samples_per_pass != 1) {
              // idx from 0:1, how does this work as uint32?
              idx /= (uint32_t) samples_per_pass;
          }
          // std::cout<<idx<<std::endl;

          // Old:
          // ref<ImageBlock> block = new ImageBlock(film_size, channels.size(),
          //                                      film->reconstruction_filter(),
          //                                      !has_aovs);
            // New:
          // ref<SignalBlock> block = new SignalBlock(adc_rd_size, channels.size(),
          //                                      adc->reconstruction_filter(),
          //                                      !has_aovs);
          ref<SignalBlock> block = new SignalBlock(adc_rd_size, channels.size(),
                                               adc->reconstruction_filter(),
                                               false);
          block->clear();

          // Maybe put noise here?

          // Modify something depending on film type.
          // Vector2f pos = Vector2f(Float(idx % uint32_t(film_size[0])),
          //                       Float(idx / uint32_t(film_size[0])));

          // When considering range doppler, there should be only 1
          // 'positionpixel'
          // We will now move across the 'pixel'
          // Pos should be of size (2,numsampes*numpixels)

          // Throws a malloc error if film size isn't 1x1
          Vector2f pos = Vector2f(Float(idx % uint32_t(1)),
                                Float(idx / uint32_t(1)));
          std::vector<Float> aovs(channels.size());

          // Perhaps after each slow time we update the scene?
          for (size_t i = 0; i < n_passes; i++) {
              receive_sample(scene, receiver, sampler, block, aovs.data(),
                          pos, diff_scale_factor);
          }

          adc->put(block);
      }

      // // Original:
      // if (!m_stop) {
      //   Log(Info, "Signal Rendering finished. (took %s)",
      //       util::time_string(m_render_timer.value(), true));
      // }
      // Supress Logs:
      if (!m_stop) {
        Log(Debug, "Signal Rendering finished. (took %s)",
            util::time_string(m_render_timer.value(), true));
      }


      return !m_stop;
}

MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  receive_block(const Scene *scene, const Receiver *receiver, Sampler *sampler,
      SignalBlock *block, Float *aovs, size_t sample_count_, size_t block_id)
      const {
    block->clear();
    uint32_t pixel_count  = (uint32_t)(m_block_size * m_block_size),
             sample_count = (uint32_t)(sample_count_ == (size_t) -1
                                           ? sampler->sample_count()
                                           : sample_count_);

    ScalarFloat diff_scale_factor =
        rsqrt((ScalarFloat) sampler->sample_count());

    if constexpr (!is_array_v<Float>) {
        for (uint32_t i = 0; i < pixel_count && !should_stop(); ++i) {
            sampler->seed(block_id * pixel_count + i);

            ScalarPoint2u pos = enoki::morton_decode<ScalarPoint2u>(i);
            if (any(pos >= block->size()))
                continue;

            pos += block->offset();
            for (uint32_t j = 0; j < sample_count && !should_stop(); ++j) {
                receive_sample(scene, receiver, sampler, block, aovs,
                              pos, diff_scale_factor);
            }
        }
    } else if constexpr (is_array_v<Float> && !is_cuda_array_v<Float>) {
        // Ensure that the sample generation is fully deterministic
        sampler->seed(block_id);

        for (auto [index, active] : range<UInt32>(pixel_count * sample_count)) {
            if (should_stop()) {
                break;
            }
            Point2u pos = enoki::morton_decode<Point2u>
                (index / UInt32(sample_count));
            active &= !any(pos >= block->size());
            pos += block->offset();
            // render_sample(scene, sensor, sampler, block, aovs,
            //     pos, diff_scale_factor, active);
            receive_sample(scene, receiver, sampler, block, aovs,
                          pos, diff_scale_factor);
        }
    } else {
        ENOKI_MARK_USED(scene);
        ENOKI_MARK_USED(receiver);
        ENOKI_MARK_USED(aovs);
        ENOKI_MARK_USED(diff_scale_factor);
        ENOKI_MARK_USED(pixel_count);
        ENOKI_MARK_USED(sample_count);
        Throw("Not implemented for CUDA arrays.");
    }
}

// The idea here is that you guarantee that you are getting samples where you
// want so that they aren't wasted. We say we are starting rays at range r and
// doppler d. Cameras say they they start at x and y. If they happen to not get
// any hits, so be it, but they tried.
// For radar that means that we start with a ray with range and doppler value,
// then must find a path that achieves this.


// // ============================================================================
// // This is the code for fmcw. I will leave commented and safe for now.
// // Doing this so that I can work on pulsed.
// // ============================================================================
// MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
//   receive_sample(const Scene *scene, const Receiver *receiver, Sampler *sampler,
//                     SignalBlock *block, Float *aovs, const Vector2f &pos,
//                     ScalarFloat diff_scale_factor, Mask active) const {
//       // We've been given a location on the pixel grid, now we jitter it
//       // In terms of signal rendering, the position should be a time pos
//       // I would have a spiralling rangedoppler block, with a sampler random
//       // within that block; then a simple random direction and position.
//       Vector2f position_sample = pos + sampler->next_2d(active);
//
//       // std::cout << position_sample << std::endl;
//
//       // If we are forcing rays to have a certain range-doppler, we don't need
//       // a tuple because the result naturally falls in the wanted bin.
//
//       // This number goes from 0-1. Keep in mind for wigner and transforms
//       Point2f aperture_sample(.5f);
//       if (receiver->needs_aperture_sample()) {
//           aperture_sample = sampler->next_2d(active);
//       }
//
//       // Random test
//
//       // Currently the renderer chooses a random time based on shutter open.
//
//       // Two approaches: If we choose a random receive time, then we need to
//       // find a path that obeys this. kinda ok for cw radar. choose a time and
//       // whatever path you find is ok. But not for pulsed radar.
//
//       // Other approach: random time based on transmit signal, time simply
//       // increases with propagation.
//
//       // block->put will have to have phase aware routines unless the receiver
//       // encodes it into the power first. I think this is better.
//
//       // As a first foray, I'd choose a random time
//
//       // Radar perception in phase-space
//       // A flex on/nod toward graham: Sensors and Signals in Phase Space
//       // Note: We always have an envelope detection with wigner..
//       // Can I pose fmcw in the phase space domain? Why not? It should work.
//       // Make a note of this in the paper~~ ^^
//
//       // What do I want my output/data struct to be?
//       // u,v,t,w...l
//       // t,w...l
//       // u,v,t...l
//       // t...l
//
//       // I want a new signal block. has position and two more. (this could be
//       // useful later for capturing 4d snapshots of wigner/light fields)
//
//       // Instead of breaking up pixels and sending blocks to the gpu, can we
//       // have multiple scene descriptions and send those blockwise to gpu?
//
//       // Is the problem many inactive lanes?
//
//       // std::cout << receiver->shutter_open() << " " << receiver->shutter_open_time() << std::endl;
//
//       // Choose a random time for the sample to be received at between receive
//       // start and end
//       Float time = receiver->shutter_open();
//       if (receiver->shutter_open_time() > 0.f) {
//           time += sampler->next_1d(active) * receiver->shutter_open_time();
//       }
//       else {
//           time = 0.f;
//       }
//
//       // time = receiver-start + block-start...lets not forget, but put on pause
//
//       // Float time = receiver->shutter_open() + sampler->next_1d(active) * receiver->shutter_open_time()
//
//       // std::cout << receiver->shutter_open_time()<< receiver->shutter_open()<< time << std::endl;
//
//       // Receive time should be limited by the adc.
//       // No...receive time is our up or down chirp time.
//       // Then bins is our fft bins.
//       // Following davids picture, time can also be fft time
//       // I think I def need to do tx and rx time.
//       // It really should be a light tracer. leave time is from tx.
//
//       // chirp t = 250 e-6 s * 16k fft
//       // Maybe I do crop window?
//       // adc_samp_rate (eg 250 MSPS (2**18)), n_bins(eg 16k fft (2**14)/16*2**10)
//       // how long is out time period of collection?
//       // Output a time frequency mix
//
//       // The adc is before the fft...obvs.
//       // But that's not what I'm doing. I need power. I do wigner.
//
//       // digital in
//
//       // n samples vs fft bins
//       // 250 MSPS * 250 us per chirp = 62.5 k Samples...much higher than 16k fft
//       // 10k samps/chirp <--- this is what goes into fpga, 2**14 = 16k
//       // Dunno how they get this number, it says complex samples, maybe?
//
//       // High res version has 0.03m res, low res (what I must have been using)
//       // should be 0.1?
//
//       // chirp time (receive time), sample rate
//       // Might need to consider 'complex' sample rate?? I guess we can take
//       // sample rate for things like pulse, and fft bins for mixed. output bins
//       // for 16k, we have sample rate 80MSPS with 200us time. For the fft do
//       // we simply decimate? Or do we use true sample rate?
//       // I think we use true samples, then it's further software for fft.
//
//       // We can get truth range doppler easily, but should we simulate bins?
//       // Perhaps calculate range/doppler bin widths from true characteristics,
//       // then mess it up based on signals? If I want a comb, I can use the bin
//       // centre.
//
//       // Alternatively, is it so difficult to get the true time?
//
//       // From the integrator, rays come back with a time/wavelength/pathlength
//
//       // Float time = receiver->receive_start();
//       // if (receiver->receive_time() > 0.f) {
//       //     time += sampler->next_1d(active) * receiver->receive_time();
//       // } else {
//       //     time = 0.f;
//       // }
//
//       // Float time = transmitter->transmit_start();
//       // if (transmitter->transmit_time() > 0.f) {
//       //     time += sampler->next_1d(active) * transmitter->transmit_time();
//       // } else {
//       //     time = 0.f;
//       // }
//
//       // Leaves with a specified time and frequency from tx clock.
//       // Goes through scene, time evolving backward. When it hits the tx, mix
//       // the signals and get result. Not exactly correct wavelength wise, but
//       // an ok guess.
//
//       // I will have to stick to fmcw. If it's a pulse I'll need a light tracer.
//
//       // The sensor and emitter are distinct from transmitter and receiver, But
//       // closely related.
//       // Does an emitter have a transmitter? Does a sensor have a receiver?
//
//       // A receiver should have an input signal/member signal.
//       // A transmitter should have an output/member signal.
//
//       // First, lets try and marry these to the current framework.
//
//       // For pulsed signals this is clear. For CW signals, we take a window of
//       // the cw signal.
//       // What do these lines mean? We need a class transmitter that has a sub
//       // class, signal.
//       // We need a class receiver.
//       // Float time = transmitter->sig_start();
//       // if (transmitter->window_time() > 0.f) {
//       //     // time += sampler->next_1d(active) * receiver->window_time();
//       //     time += sampler->next_1d(active) * transmitter->window_time();
//       // } else {
//       //     time = 0.f;
//       // }
//
//       // two different ranges
//       // 1) the overall range of the sim, ie wideband/94GHZ +- BW
//       // 2) the ranges of frequency that each element has their radiance/
//       // reflectivity/absorption defined at.
//
//       // Float frequency_sample = sampler->next_1d(active);
//       // Wavelength/time is not random. Do I do a mapping, or a sampling.
//       // I guess I need to do a ray-weight as well.
//       // This will become deterministic plus an offset. It will become a beat
//       // frequency, and then a wavelength.
//       Float wavelength_sample = sampler->next_1d(active);
//
//       // The wavelength should be the tx wavelength at this time
//       // Float f(0.f);
//       // Float f0 = 94e9 - 6e9/2;
//       // Float f1 = 94e9 + 6e9/2;
//       //
//       // Float t1 = 240e-6;
//       // Float t2 = t1 + 10e-6;
//       // Float t3 = t2 + 240e-6;
//       // Float t4 = t3 + 10e-6;
//       //
//       // // Float tn = math::modulo(time, t4);
//       // Float tn = time;
//       //
//       // // These are definitely the culprit
//       // if (all(tn < t1)) {
//       //     // f = 2*((f1 - f0)/(t2 - t1))*tn + f0;
//       //     // f = 2*((6e9)/(240e-6))*tn + f0;
//       //     f = ((6e9)/(240e-6))*tn + f0;
//       // } else if (all(tn < t2)) {
//       //     f = f1;
//       // } else if (all(tn < t3)){
//       //     // f = 2*((f0 - f1)/(t3 - t2))*tn + f1;
//       //     // f = 2*((-6e9)/(240e-6))*tn + f1;
//       //     f = ((-6e9)/(240e-6))*tn + f1;
//       // } else {
//       //     f = f0;
//       // }
//       //
//       //
//       // // Float wavelength_sample = math::CVac<float>/f;
//       // // Float wavelength_sample = (f-f0)/(f1-f0);
//       // // Wavelength sample is reversed from frequency.
//       // Float wavelength_sample = 1.f - (f-f0)/(f1-f0);
//
//
//       // std::cout << tn *1e6 << " " << wavelength_sample <<  " " << f << std::endl;
//
//       // This is just a real (local) position.
//       // Get the sample from 0-1, take away the
//       // Old:
//       // Vector2f adjusted_position =
//       //   (position_sample - sensor->film()->crop_offset()) /
//       //   sensor->film()->crop_size();
//         // Perhaps sensor can have film & ADC?
//       // New:
//       Vector2f adjusted_position = position_sample;
//
//       // ray_weight will be like the bsdf of the sensor. I can also use it to
//       // apply my phase shift. Maybe not, this more relates to the wavelength.
//       // Could make things more complicated, pass a frequency sample, then
//       // convert to wavelength inside using medium and bandwidth conversion.
//       // std::cout << time << wavelength_sample << adjusted_position << aperture_sample << std::cout;
//
//       // std::cout << aperture_sample << std::endl;
//
//       // In the context of time, this differential should be which t,r/f,d bin
//       // This should probably take a frequency sample, then change it based
//       // on receiver characteristics, eg matched filter or receiver bandwidth
//       auto [ray, ray_weight] = receiver->sample_ray_differential(
//           time, wavelength_sample, adjusted_position, aperture_sample);
//       // this says that the incoming ray has a wavelength and time already.
//       // These are the rays landing on the sensor, we propagate back into the
//       // scene, lets interpret t and ?? as tx values.
//
//       // I can't know the wavelength/freq until I know the travel time.
//
//       // This version guarentees a pixel hit, but time/freq is unconstrained
//
//       // Here we make a copy of ray.
//
//       ray.scale_differential(diff_scale_factor);
//
//       // std::cout << ray.wavelengths << std::endl;
//
//       const Medium *medium = receiver->medium();
//       // std::tuple<Spectrum, Mask, Float> result =
//       //   sample(scene, sampler, ray, medium, aovs+1, active);
//       // When doing this, we only render the 'yellow'
//       // std::tuple<Spectrum, Mask, Float> result =
//       // sample(scene, sampler, ray, medium, aovs+5, active);
//       // std::tuple<Spectrum, Mask, Float> result =
//       // sample(scene, sampler, ray, medium, aovs+3, active);
//
//       // std::cout<<ray.time<<std::endl;
//
//       // std::pair<Spectrum, Mask> result = sample(scene, sampler, ray, medium,
//       //     aovs + 3, active);
//
//       // ray* rp = &ray;
//
//       // std::pair<Spectrum, Mask> result = sample(scene, sampler, &ray, medium,
//       //     aovs + 3, active);
//
//       Wavelength f_rx(0.f);
//       Wavelength f0 = 94e9 - 6e9/2;
//       Wavelength f1 = 94e9 + 6e9/2;
//
//       Float t1 = 240e-6;
//       Float t2 = t1 + 10e-6;
//       Float t3 = t2 + 240e-6;
//       Float t4 = t3 + 10e-6;
//
//       // std::cout << "int in" << std::endl;
//       Float tn = math::fmodulo(time, t4);
//       // Float tn = math::modulo(time, t4);
//       // std::cout << "int out" << std::endl;
//       // Float tn = time;
//
//       // These are definitely the culprit
//       if (all(tn < t1)) {
//           // f = 2*((f1 - f0)/(t2 - t1))*tn + f0;
//           // f = 2*((6e9)/(240e-6))*tn + f0;
//           f_rx = ((6e9)/(240e-6))*tn + f0;
//           // std::cout << "t1: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else if (all(tn < t2)) {
//           f_rx = f1;
//           // std::cout << "t2: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else if (all(tn < t3)){
//           // f = 2*((f0 - f1)/(t3 - t2))*tn + f1;
//           // f = 2*((-6e9)/(240e-6))*tn + f1;
//           f_rx = ((-6e9)/(240e-6))*(tn - t2) + f1;
//           // std::cout << "t3: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else {
//           f_rx = f0;
//           // std::cout << "t4: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       }
//
//
//       // std::cout << time << " " << tn << std::endl;
//
//       // Float wavelength_sample = math::CVac<float>/f;
//       // Float wavelength_sample = (f-f0)/(f1-f0);
//       // Wavelength sample is reversed from frequency.
//       //
//       Wavelength lambda_rx = math::CVac<double>/f_rx * 1e9;
//       ray.wavelengths = lambda_rx;
//
//       // // Now evaluate doppler, this sometimes crashes
//       // if(any_or<true>(neq(receiver, nullptr))) {
//       //     SurfaceInteraction3f si_rx;
//       //     si_rx.p = ray.o;
//       //     si_rx.wavelengths = ray.wavelengths;
//       //     si_rx.wi = ray.d;
//       //     si_rx.shape = receiver->shape();
//       //     ray.wavelengths += select(neq(receiver, nullptr), receiver->doppler(si_rx, active), 0.f);
//       // }
//
//       // ray.wavelengths = math::CVac<double>/f_rx * 1e9;
//
//       // std::cout<<"??: " << ray.wavelengths << "f: " << f<<std::endl;
//
//       // For the paper I've proposed, I need to address interference, and
//       // appropriate sampling. Not range. That is later.
//
//       // Example for pulse. Because of the way we sample the block, we already
//       // know the output frequency and time.
//       // From the resulting frequency we calculate the wavelength, and this is
//       // what we pass to sample. Time is similarly updated due to line delays.
//       // The ray bounces around and eventually hits the tx. The wavelength is
//       // converted to a frequency, and the arrival time is also updated from
//       // delay line. The tx power is then queried at time t, frequency f and the
//       // ray gain updated. If we're at a time when no pulse, sucked in.
//
//
//       // For an fmcw radar. We sample a beat frequency and time. We get the
//       // frequency from the transmitter via time and sampling. Receive freq =
//       // tx_freq - beat. We can now find the wavelength and toa at the receiver.
//       // Ray leaves rx at time and with ??, and propagates, eventually landing
//       // on the transmitter at t, ??. Tx power is queried at t, ??. If no value,
//       // sucked in.
//
//       // As part of this problem, we enter the sampling routine with a ray
//       // wavelength and time. For now assume no doppler shifts. We can find the
//       // round trip time t by querying the tx signal with wavelength of rx.
//       // This will give us a time in the signal domain. Naively lets start with
//       // a non-repeating signal. The delay is now the difference between signal
//       // time and ray initial time (plus line delays). This gives us a path
//       // length. We must now find a path that satisfies these constraints.
//       // Sounds like an optimisation problem, find a path/set of nodes such
//       // that return is maximised. Beware: there are many zeros.
//
//       // In both cases, our challenge is to find a set of paths which maximise
//       // the return and exploration of search space. Our challenge is that it is
//       // kinda constrained, kinda not.
//
//       // We currently have an unconstrained version. In this version (for fmcw)
//       // the sampling routine gets the transmitter wavelength at t. We can give
//       // it any wavelength, and still grab from tx, only t/lambda matters.
//
//       // ATM we randomly sample t, and directly sample ??. The next level of
//       // smart sampling will be to have a minimum path length.
//
//       // Either way, the receiver time t, and wavelength ?? should be true; or at
//       // least one of them.
//
//       // Next challenge: Create a signal plugin for transmitter. After that, use
//       // it in integrator. Update receiver sample ray differential. Let receiver
//       // get signal from tx as needed. Then, spin my radar around and make a
//       // pretty picture.
//
//       // How should I represent this data? There is a convenient 2d projection,
//       // slow-time/range.
//       // One simulation is a quasistatic collection. In the simplest form, it is
//       // a series of fft blocks (these are coherently summed for a chirp anyway)
//       // Lets say we have the high res, high noise block, and our quasistatic
//       // simulation includes n such blocks...why bother?
//
//       // Aovs have multiple layers. 1st: the multiple bins of fft from a chirp.
//       // 2nd: Up and down chirps. A cycle.
//       // 3rd: Multi-chirps.
//
//       // Finally we have the data structure:
//       // x->slow time
//       // y->beat (fast time->range)(range: average of all aovs)
//       // u->slow doppler
//       // v->fast doppler (difference/correlation of up downs)
//
//
//
//       // 5d function representations: x, y position. r, g direction. b intensity
//       // Naturally these have to be normalised.
//       // This way a flat wdf can be easily represented. Actually no: remember
//       // the big challenge: A single position can have many directions.
//
//       // Like our 4d wdf (4 axes + 1 value = 5) and projections:
//
//       // topleft: slow time vs total power
//       // middleleft: slow time vs range
//       // bottomleft: slow time vs macro doppler
//       // topmiddle: micro doppler vs total power
//       // middlemiddle: micro doppler vs range
//       // bottommiddle: micro doppler vs macro doppler
//       // rightmiddle: total power vs range
//       // rightbottom: total power vs macro doppler
//
//       // std::cout << ray.wavelengths << std::endl;
//       // It knows where to start propagation via ray location.
//
//       // Mixed:
//       // Output decides beat freq and time.
//       // Call receiver->sample_ray_differential. Inside this we will have some
//       // cool tx and mixing shit. It will return a ray with time t, and rx
//       // wavelength. At the moment the rx ?? means shit all-> because the
//       // integrator will overwrite the wavelength.
//
//       // Pulsed:
//       // Call receiver->sample_ray_differential. Doesn't know about tx signal.
//       // Happily gives a ray wavelength and time in band bounds.
//
//       std::pair<Spectrum, Mask> result = sample(scene, sampler, ray, medium,
//           aovs + 3, active);
//
//           // this is returning an intensity, but also a modified ray. This ray
//           // will have an updated time and wavelength.
//
//           // Then do receiver-> process.
//           // receiver MAY have a signal, gets it from transmitter.
//           // Now, what do we do if we don't have a signal?
//
//           // Each 'block' should be its own render
//           // receiver bandwidth = sps
//           // receiver df = bandwidth / bins
//           // receiver dt = 1 / bandwidth
//           // receiver t = bins / bandwidth
//
//           // t0 should come from the external render block.
//
//           // (ray.time - t0)
//           //
//           // tf_pos /= bandwidth * bins
//
//       // Make a new sample routine: one which allows the ray to be changed.
//       // In it, the ray's wavelength can change by hitting moving targets. It's
//       // time changes with propagation time and length (ray.t) changes with
//       // propagation distance. ray.t can be also updated to include phase
//       // changes at boundaries, or simply be phase.
//       // Maybe I don't even need...
//
//       // ScalarVector2i rd = (std::get<2>(result)/
//       //   (sensor->far_clip()-sensor->near_clip()), 0);
//
//       // We'll have to get this from our signal/film prop.
//       // Alternatively from our integrator.
//       // ScalarVector2i rd = (std::get<2>(result)/
//       //   (10-0.1), 0);
//       // Vector2f rd = {ray.time, ray.wavelengths[0]};
//
//       Wavelength f_tx(0.f);
//       Float tn2 = math::fmodulo(ray.time, t4);
//       // Float tn = math::modulo(time, t4);
//       // std::cout << "int out" << std::endl;
//       // Float tn = time;
//
//       // These are definitely the culprit
//       if (all(tn2 < t1)) {
//           // f = 2*((f1 - f0)/(t2 - t1))*tn + f0;
//           // f = 2*((6e9)/(240e-6))*tn + f0;
//           f_tx = ((6e9)/(240e-6))*tn2 + f0;
//           // std::cout << "t1: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else if (all(tn2 < t2)) {
//           f_tx = f1;
//           // std::cout << "t2: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else if (all(tn2 < t3)){
//           // f = 2*((f0 - f1)/(t3 - t2))*tn + f1;
//           // f = 2*((-6e9)/(240e-6))*tn + f1;
//           f_tx = ((-6e9)/(240e-6))*(tn2 - t2) + f1;
//           // std::cout << "t3: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       } else {
//           f_tx = f0;
//           // std::cout << "t4: tn: " << tn << " t: " << time << " f: " << f_rx[0] << std::endl;
//       }
//
//       // f_tx += (ray.wavelengths - lambda_rx);
//       // std::cout << ray.wavelengths/lambda_rx << std::endl;
//       f_tx *= (1/(ray.wavelengths/lambda_rx));
//
//       Vector2f rd = {time - receiver->shutter_open(), abs(f_tx[0]-f_rx[0])};
//
//
//       // Vector2f rd = {time - receiver->shutter_open(), math::CVac<float>/ray.wavelengths[0]*1e-9};
//       // Vector2f rd = {time, math::CVac<float>/ray.wavelengths[0] *1e-9};
//       // Vector2f rd = {time, ray.wavelengths[0]};
//
//       // rd = (rd - (receiver->adc()->centres() - receiver->adc()->bandwidth()/2))/receiver->adc()->bandwidth() * receiver->adc()->size();
//
//       // We get real value, eg time, range, doppler, freq
//       // rd = (rd - (receiver->adc()->centres() - receiver->adc()->bandwidth()/2))
//       //   /receiver->adc()->bandwidth() * receiver->adc()->size();
//
//       // std::cout << rd << " " << rd[0]* receiver->adc()->size().x() / receiver->adc()->bandwidth().x() << std::endl;
//       // std::cout << receiver->adc()->size().x() << receiver->adc()->bandwidth().x() << std::endl;
//
//       // std::cout << rd << std::endl;
//
//       rd *= receiver->adc()->size() / receiver->adc()->bandwidth();
//       // rd[0] *= receiver->adc()->size().x() / receiver->adc()->bandwidth().x();
//       // rd[1] *= receiver->adc()->size().y() / receiver->adc()->bandwidth().y();
//
//       // std::cout << rd << std::endl;
//
//         // val/range * bins
//
//         // Get time/freq from sim
//         // call receiver->process
//
//         // Basically I need a method which optionally takes 2 rays and converts
//         // their time/wavelength into bin indices.
//
//         // receiver->process()
//
//
//       // rd = (rd - (receiver->adc()->centres() - receiver->adc()->bandwidth()/2))/receiver->adc()->bandwidth() * receiver->adc()->window_size();
//       // rd = (rd - (receiver->adc()->centres() - receiver->adc()->bandwidth()/2))/receiver->adc()->bandwidth();
//
//       // 0-1 location
//       // rd[0] = (ray.time - (receiver->adc()->centres().x() - receiver->adc()->bandwidth().x()/2))/receiver->adc()->bandwidth().x();
//       // rd[1] = (ray.wavelengths[0] - (receiver->adc()->centres().y() - receiver->adc()->bandwidth().y()/2))/receiver->adc()->bandwidth().y();
//
//       // // 0-full location
//       // rd[0] = (ray.time - (receiver->adc()->centres().x() - receiver->adc()->bandwidth().x()/2)) * receiver->adc()->size().x()/receiver->adc()->bandwidth().x();
//       // rd[1] = (ray.wavelengths[0] - (receiver->adc()->centres().y() - receiver->adc()->bandwidth().y()/2)) * receiver->adc()->size().y()/receiver->adc()->bandwidth().y();
//
//       // rd*=receiver->adc()->size();
//       // rd*=receiver->adc()->window_size();
//       // rd = (rd - receiver->adc()->window_offset()) / receiver->adc()->window_size();
//       // rd = (rd + receiver->adc()->window_offset());
//       // rd[0] = (rd[0] - receiver->adc()->window_offset().x());
//       // rd[1] = (rd[1] - receiver->adc()->window_offset().y());
//       // rd[1] = (rd[1] + receiver->adc()->window_offset().y());
//       // rd[1] = rd[1] +  receiver->adc()->window_size().y() - receiver->adc()->window_offset().y();
//       // rd[1] = rd[1] - receiver->adc()->window_offset().y();
//       // rd[1] = rd[1];
//
//       // rd -= receiver->adc()->window_offset();
//       // rd[0] -= receiver->adc()->window_offset().y();
//       // rd[1] -= receiver->adc()->window_offset().y();
//
//       // This work?
//       // block->set_offset(receiver->adc()->window_offset());
//       // block->set_size(receiver->adc()->window_size());
//
//
//
//
//
//       // std::cout << rd << ray.wavelengths[0] << receiver->adc()->window_size() << receiver->adc()->size() <<std::endl;
//
//       // Vector2f adjusted_position =
//       //     (position_sample - sensor->film()->crop_offset()) /
//       //     sensor->film()->crop_size();
//
//       // std::cout << rd << std::endl;
//
//       // The point of this function is to get a real value and convert it to a
//       // bin idx. This should be similar to the hardware and include things
//       // like mixing, filtering etc. Note this is a single return value, the
//       // effects of multiple returns comes later.
//       //
//       // Fuck this, I think it needs to be t/f
//       // The internal representation is tf. My only concern is high res.
//       // Or is it just always magically delay time?
//       // Receiver n_bins
//       // Can I skip ahead or do i lose important info? I need to put things
//       // into the correct time bins because I need to account for interference.
//
//       // if we have a pulse, it should come in with a spread frequency.
//       // receiver->process_signal()
//
//       // Change this back, but allow rays to be modified.
//       // std::pair<Spectrum, Mask> result =
//       // sample(scene, sampler, ray, medium, aovs+5, active);
//       // Gives power, wavelength.
//
//         // std::cout << std::get<0>(result) << std::endl;
//
//         // We currently have a path tracer. ie ray time is at rx.
//
//         // Unconstrained radar path tracer.
//         // t_rx is receive time
//         // w_rx is transmit/clock frequency now.
//         // We don't know what the wavelength will be, because the path time is
//         // unknown. Lets send the tx signal NOW, and do everything backward.
//
//         // propagate through scene.
//         // when it hits the transmitter, get the time, and tx freq.
//         // apply 'mixnfilter'. just do incoherent atm.  to do coherent, we'd
//         // need to append each ray landing in a block then do combination sum.
//
//         // Each ray gets multiplied by a wigner layer wrt the path length..or t
//         // Do we include bounces? Or does that already happen when we hit the
//         // ground? It is simply a -1 after all. If we do path length, that
//         // assumes that all rays left with 0 phase. Phase will need to be
//         // calculated based on time, we can say 0 phase at time 0, it's all
//         // relative anyway....maybe this will even give interference between
//         // targets!!
//
//         // I could put phase calc here, but that would force it to be always
//         // be wigner. Maybe that's ok considering we're saying that signals are
//         // tf anyway. How does this limit my ability to do diff rendering?
//         // Let's build the code properly. Keep classic mods so that maybe we
//         // can go back.
//       result.first = ray_weight * result.first;
//       UnpolarizedSpectrum spec_u = depolarize(result.first);
//
//       // std::cout<<spec_u<<std::endl;
//       // spec_u should be like radiances at wavelengths
//
//       // This should return a simple power at 1 distance. The accumulated power
//       // of the entire sim
//
//       Color3f xyz;
//       if constexpr (is_monochromatic_v<Spectrum>) {
//           xyz = spec_u.x();
//       } else if constexpr (is_rgb_v<Spectrum>) {
//           xyz = srgb_to_xyz(spec_u, active);
//       } else {
//           static_assert(is_spectral_v<Spectrum>);
//           // This is in receive sample, so should be only used with adc film
//           // xyz = spectrum_to_xyz(spec_u, ray.wavelengths, active);
//           xyz = spec_u.x();
//           // std::cout<<xyz<<std::endl;
//       }
//       // if constexpr (is_wigner_v<Spectrum>){}
//
//       // What should I output?
//       // I have wavelength
//
//       // float freq = 94E9;
//       // float lambda = math::CVac/freq;
//
//       // aovs[0] = xyz.x();
//
//       // aovs[1] = lambda;
//       // aovs[1] = select(std::get<1>(result), Float(1.f), Float(0.f));
//       // aovs[1] = xyz.y();
//       // aovs[2] = xyz.z();
//       // aovs[3] = select(std::get<1>(result), Float(1.f), Float(0.f));
//       // aovs[4] = 1.f;
//
//       // aovs[0] = xyz.x();
//       // aovs[1] = xyz.y();
//       // aovs[2] = xyz.z();
//       // aovs[3] = select(std::get<1>(result), Float(1.f), Float(0.f));
//       // aovs[4] = 1.f;
//
//       // aovs[0] = xyz.x();
//       // aovs[1] = select(std::get<1>(result), Float(1.f), Float(0.f));
//       // aovs[2] = 1.f;
//
//       aovs[0] = xyz.x();
//       aovs[1] = select(result.second, Float(1.f), Float(0.f));
//       aovs[2] = 1.f;
//       // aovs[3] = cos(phase)
//
//       // For interference, what about a weighted sum of phase? Each ray
//       // contributes its phase * power. Or If all rays were equal, phase
//       // What about tracking phase variance?
//       // Or is it another aov?
//       // Normalised power*phase.
//
//       // Modify sample so that ray is not const.The sample routine can now
//       // change the ray.time and ray.wavelength.
//       // aovs become result(basically power), time, wavelength....or can I
//       // encode these factors into what they already are?
//
//       // tf_sample = ray.time, ray.wavelengths/frequency
//
//       // block->put(tf_sample, aovs, active);
//
//       // Old:
//       // block->put(position_sample, aovs, active);
//       // New:
//       block->put(rd, aovs, active);
//
//       // This should be for coherent sum
//       // block->append(rd, aovs, active);
//
//       // block should have tx channel and rx channel. These are aovs...maybe
//       // and nah, cause too fast.
//
//       // block->put(position_sample, time_bin (as a fraction of the span), frequency_bin, aovs, active);
//
//       sampler->advance();
// }

// ============================================================================
// Code for pulsed radar
// ============================================================================
// remember in proper framework, pos encodes t, f.
MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  receive_sample(const Scene *scene, const Receiver *receiver, Sampler *sampler,
                    SignalBlock *block, Float *aovs, const Vector2f &pos,
                    ScalarFloat diff_scale_factor, Mask active) const {

      // Spatial sampling -------------------------
      Vector2f position_sample = pos + sampler->next_2d(active);
      Vector2f adjusted_position = position_sample;
      // ==========================================

      // Direction sampling -----------------------
      Point2f aperture_sample(.5f);
      if (receiver->needs_aperture_sample()) {
          aperture_sample = sampler->next_2d(active);
      }
      // ==========================================

      // Time sampling ----------------------------
      Float time = receiver->adc_sampling_start();
      if (receiver->adc_sampling_time() > 0.f) {
          time += sampler->next_1d(active) * receiver->adc_sampling_time();
      } else {
          time = 0.f;
      }
      // ==========================================

      // Wavelength sampling ----------------------
      Float wavelength_sample = sampler->next_1d(active);
      // ==========================================

      // Generate ray -----------------------------
      auto [ray, ray_weight] = receiver->sample_ray_differential(
          time, wavelength_sample, adjusted_position, aperture_sample);
      ray.scale_differential(diff_scale_factor);
      // Ray should have phase = 0, unless it comes from mixer,
      // then it is whatever the tx is
      // ==========================================

      // Assign medium ----------------------------
      const Medium *medium = receiver->medium();
      // ==========================================

      // Prepare result and location --------------
      std::pair<Spectrum, Mask> result;
      Vector2f tf;
      // ==========================================

     // ray_weight = 1.f;
     // ray.maxt = ray_weight[0];

      if (receiver->receive_type() == "mix_resample") {
          // Save the receive frequency -----------
          Wavelength f_rx = MTS_C*rcp(ray.wavelengths*1e-9);
          // ======================================
          // Propagate through the scene ----------
          result = sample(scene, sampler, ray, medium,
              aovs + 3, active);
          // if (ray.phase != 0){
          //     std::cout << ray.phase << std::endl;
          // }
          // ======================================
          // Find the tf/beat freq ----------------
          tf[0] = time - receiver->adc_sampling_start();
          tf[1] = abs(MTS_C*rcp(ray.wavelengths[0]*1e-9)-f_rx[0]);
          // std::cout << tf << std::endl;
          // ======================================

      } else if (receiver->receive_type() == "raw") {
          // Propagate through the scene ----------
          result = sample(scene, sampler, ray, medium,
              aovs + 3, active);
          // ray goes through scene, phase decrementing
          // ray eventually hits transmitter, transmitter has some phase at time t.
          // phase -= phasetx
          // ======================================
          // Extract the doppler shifted freq -----
          tf[0] = time - receiver->adc_sampling_start();
          tf[1] = MTS_C*rcp(ray.wavelengths[0]*1e-9);
          // ======================================
      } else if (receiver->receive_type() == "raw_resample") {
          // Propagate through the scene ----------
          result = sample(scene, sampler, ray, medium,
              aovs + 3, active);
          // ======================================
          // Extract the doppler shifted freq -----
          tf[0] = time - receiver->adc_sampling_start();
          tf[1] = MTS_C*rcp(ray.wavelengths[0]*1e-9);
          // ======================================
      } else if (receiver->receive_type() == "mixer") {
          // if receiver->type == mixer
          // unconstrained:
          // choose random time
          // choose random beat frequency
          // find the corresponding tx frequency at time trx (can select from delta)
          // get the rx frequency = abs(txf - beat)
          // propagate randomly, arrives at tx at t_tx, f_tx~f_rx
          // evaluate tx gain at time t_tx, f_tx
          // will likely be 0 because we did not control pathlength
      }

      // Adjust t/f to pixels ---------------------
      tf *= receiver->adc()->size() / receiver->adc()->bandwidth();
      // ==========================================

      // Ignore receiver weight for now
      result.first = abs(ray_weight) * result.first;
      // result.first = ray_weight * result.first;

      // result.first = abs(ray_weight * result.first);
      // result.first = 1.f * result.first;

      // UnpolarizedSpectrum spec_u = depolarize(result.first);
      // Color3f xyz;
      // if constexpr (is_monochromatic_v<Spectrum>) {
      //     xyz = spec_u.x();
      // } else if constexpr (is_rgb_v<Spectrum>) {
      //     xyz = srgb_to_xyz(spec_u, active);
      // } else {
      //     static_assert(is_spectral_v<Spectrum>);
      //     xyz = spec_u.x();
      // }
      // aovs[0] = xyz.x();

      aovs[0] = select(result.second, hsum(result.first), Float(0.f));
      aovs[1] = select(result.second, Float(1.f), Float(0.f));
      aovs[2] = 1.f;

      block->put(tf, aovs, active);
      sampler->advance();
}

MTS_VARIANT std::pair<Spectrum, typename SamplingIntegrator<Float, Spectrum>::Mask>
SamplingIntegrator<Float, Spectrum>::sample(const Scene * /* scene */,
                                            Sampler * /* sampler */,
                                            const RayDifferential3f & /* ray */,
                                            const Medium * /* medium */,
                                            Float * /* aovs */,
                                            Mask /* active */) const {
    NotImplementedError("sample");
}

// MTS_VARIANT std::pair<Spectrum, typename SamplingIntegrator<Float, Spectrum>::Mask>
// SamplingIntegrator<Float, Spectrum>::sample(const Scene * /* scene */,
//                                             Sampler * /* sampler */,
//                                             RayDifferential3f & /* ray */,
//                                             const Medium * /* medium */,
//                                             Float * /* aovs */,
//                                             Mask /* active */) const {
//     NotImplementedError("sample");
// }

// MTS_VARIANT std::pair<Spectrum, typename SamplingIntegrator<Float, Spectrum>::Mask>
// SamplingIntegrator<Float, Spectrum>::sample(const Scene * /* scene */,
//                                             Sampler * /* sampler */,
//                                             RayDifferential3f * /* ray */,
//                                             const Medium * /* medium */,
//                                             Float * /* aovs */,
//                                             Mask /* active */) const {
//     NotImplementedError("sample");
// }

// MTS_VARIANT std::tuple<Spectrum,
//                         typename SamplingIntegrator<Float, Spectrum>::Mask,
//                         Float>SamplingIntegrator<Float, Spectrum>::
//   sample(const Scene * /* scene */,
//             Sampler * /* sampler */,
//             const RayDifferential3f & /* ray */,
//             const Medium * /* medium */,
//             Float * /* aovs */,
//             Mask /* active */) const {
//       NotImplementedError("sample");
// }

// -----------------------------------------------------------------------------

MTS_VARIANT MonteCarloIntegrator<Float, Spectrum>::
  MonteCarloIntegrator(const Properties &props) : Base(props) {
      /// Depth to begin using russian roulette
      m_rr_depth = props.int_("rr_depth", 5);
      if (m_rr_depth <= 0) {
        Throw("\"rr_depth\" must be set to a value greater than zero!");
      }

      /*  Longest visualized path depth (``-1 = infinite``). A value of \c 1 will
          visualize only directly visible light sources. \c 2 will lead to
          single-bounce (direct-only) illumination, and so on. */
      m_max_depth = props.int_("max_depth", -1);
      if (m_max_depth < 0 && m_max_depth != -1) {
          Throw("\"max_depth\" must be set to -1 (infinite) or a value >= 0");
      }
}

MTS_VARIANT MonteCarloIntegrator<Float, Spectrum>::
  ~MonteCarloIntegrator() {}

MTS_IMPLEMENT_CLASS_VARIANT(Integrator, Object, "integrator")
MTS_IMPLEMENT_CLASS_VARIANT(SamplingIntegrator, Integrator)
MTS_IMPLEMENT_CLASS_VARIANT(MonteCarloIntegrator, SamplingIntegrator)

MTS_INSTANTIATE_CLASS(Integrator)
MTS_INSTANTIATE_CLASS(SamplingIntegrator)
MTS_INSTANTIATE_CLASS(MonteCarloIntegrator)
NAMESPACE_END(mitsuba)
