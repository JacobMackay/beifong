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
          block->clear();
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

MTS_VARIANT void SamplingIntegrator<Float, Spectrum>::
  receive_sample(const Scene *scene, const Receiver *receiver, Sampler *sampler,
                    SignalBlock *block, Float *aovs, const Vector2f &pos,
                    ScalarFloat diff_scale_factor, Mask active) const {
      // We've been given a location on the pixel grid, now we jitter it
      // In terms of signal rendering, the position should be a time pos
      // I would have a spiralling rangedoppler block, with a sampler random
      // within that block; then a simple random direction and position.
      Vector2f position_sample = pos + sampler->next_2d(active);

      // std::cout << position_sample << std::endl;

      // If we are forcing rays to have a certain range-doppler, we don't need
      // a tuple because the result naturally falls in the wanted bin.

      // This number goes from 0-1. Keep in mind for wigner and transforms
      Point2f aperture_sample(.5f);
      if (receiver->needs_aperture_sample()) {
          aperture_sample = sampler->next_2d(active);
      }

      // Random test

      // Currently the renderer chooses a random time based on shutter open.

      // Two approaches: If we choose a random receive time, then we need to
      // find a path that obeys this. kinda ok for cw radar. choose a time and
      // whatever path you find is ok. But not for pulsed radar.

      // Other approach: random time based on transmit signal, time simply
      // increases with propagation.

      // block->put will have to have phase aware routines unless the receiver
      // encodes it into the power first. I think this is better.

      // As a first foray, I'd choose a random time

      // Radar perception in phase-space
      // A flex on/nod toward graham: Sensors and Signals in Phase Space
      // Note: We always have an envelope detection with wigner..
      // Can I pose fmcw in the phase space domain? Why not? It should work.
      // Make a note of this in the paper~~ ^^

      // What do I want my output/data struct to be?
      // u,v,t,w...l
      // t,w...l
      // u,v,t...l
      // t...l

      // I want a new signal block. has position and two more. (this could be
      // useful later for capturing 4d snapshots of wigner/light fields)

      // Instead of breaking up pixels and sending blocks to the gpu, can we
      // have multiple scene descriptions and send those blockwise to gpu?

      Float time = receiver->shutter_open();
      if (receiver->shutter_open_time() > 0.f) {
          time += sampler->next_1d(active) * receiver->shutter_open_time();
      } else {
          time = 0.f;
      }

      // Receive time should be limited by the adc.
      // No...receive time is our up or down chirp time.
      // Then bins is our fft bins.
      // Following davids picture, time can also be fft time
      // I think I def need to do tx and rx time.
      // It really should be a light tracer. leave time is from tx.

      // chirp t = 250 e-6 s * 16k fft
      // Maybe I do crop window?
      // adc_samp_rate (eg 250 MSPS (2**18)), n_bins(eg 16k fft (2**14)/16*2**10)
      // how long is out time period of collection?
      // Output a time frequency mix

      // The adc is before the fft...obvs.
      // But that's not what I'm doing. I need power. I do wigner.

      // digital in

      // n samples vs fft bins
      // 250 MSPS * 250 us per chirp = 62.5 k Samples...much higher than 16k fft
      // 10k samps/chirp <--- this is what goes into fpga, 2**14 = 16k
      // Dunno how they get this number, it says complex samples, maybe?

      // High res version has 0.03m res, low res (what I must have been using)
      // should be 0.1?

      // chirp time (receive time), sample rate
      // Might need to consider 'complex' sample rate?? I guess we can take
      // sample rate for things like pulse, and fft bins for mixed. output bins
      // for 16k, we have sample rate 80MSPS with 200us time. For the fft do
      // we simply decimate? Or do we use true sample rate?
      // I think we use true samples, then it's further software for fft.

      // We can get truth range doppler easily, but should we simulate bins?
      // Perhaps calculate range/doppler bin widths from true characteristics,
      // then mess it up based on signals? If I want a comb, I can use the bin
      // centre.

      // Alternatively, is it so difficult to get the true time?

      // From the integrator, rays come back with a time/wavelength/pathlength

      // Float time = receiver->receive_start();
      // if (receiver->receive_time() > 0.f) {
      //     time += sampler->next_1d(active) * receiver->receive_time();
      // } else {
      //     time = 0.f;
      // }

      // Float time = transmitter->transmit_start();
      // if (transmitter->transmit_time() > 0.f) {
      //     time += sampler->next_1d(active) * transmitter->transmit_time();
      // } else {
      //     time = 0.f;
      // }

      // Leaves with a specified time and frequency from tx clock.
      // Goes through scene, time evolving backward. When it hits the tx, mix
      // the signals and get result. Not exactly correct wavelength wise, but
      // an ok guess.

      // I will have to stick to fmcw. If it's a pulse I'll need a light tracer.

      // The sensor and emitter are distinct from transmitter and receiver, But
      // closely related.
      // Does an emitter have a transmitter? Does a sensor have a receiver?

      // A receiver should have an input signal/member signal.
      // A transmitter should have an output/member signal.

      // First, lets try and marry these to the current framework.

      // For pulsed signals this is clear. For CW signals, we take a window of
      // the cw signal.
      // What do these lines mean? We need a class transmitter that has a sub
      // class, signal.
      // We need a class receiver.
      // Float time = transmitter->sig_start();
      // if (transmitter->window_time() > 0.f) {
      //     // time += sampler->next_1d(active) * receiver->window_time();
      //     time += sampler->next_1d(active) * transmitter->window_time();
      // } else {
      //     time = 0.f;
      // }

      // two different ranges
      // 1) the overall range of the sim, ie wideband/94GHZ +- BW
      // 2) the ranges of frequency that each element has their radiance/
      // reflectivity/absorption defined at.

      // Float frequency_sample = sampler->next_1d(active);
      // Wavelength/time is not random. Do I do a mapping, or a sampling.
      // I guess I need to do a ray-weight as well.
      Float wavelength_sample = sampler->next_1d(active);

      // This is just a real (local) position.
      // Get the sample from 0-1, take away the
      // Old:
      // Vector2f adjusted_position =
      //   (position_sample - sensor->film()->crop_offset()) /
      //   sensor->film()->crop_size();
        // Perhaps sensor can have film & ADC?
      // New:
      Vector2f adjusted_position = position_sample;

      // ray_weight will be like the bsdf of the sensor. I can also use it to
      // apply my phase shift. Maybe not, this more relates to the wavelength.
      // Could make things more complicated, pass a frequency sample, then
      // convert to wavelength inside using medium and bandwidth conversion.
      // std::cout << time << wavelength_sample << adjusted_position << aperture_sample << std::cout;

      // std::cout << aperture_sample << std::endl;

      // In the context of time, this differential should be which t,r/f,d bin
      auto [ray, ray_weight] = receiver->sample_ray_differential(
          time, wavelength_sample, adjusted_position, aperture_sample);
      // this says that the incoming ray has a wavelength and time already.
      // These are the rays landing on the sensor, we propagate back into the
      // scene, lets interpret t and λ as tx values.

      // I can't know the wavelength/freq until I know the travel time.

      // This version guarentees a pixel hit, but time/freq is unconstrained

      // Here we make a copy of ray.

      ray.scale_differential(diff_scale_factor);

      // std::cout << ray << std::endl;

      const Medium *medium = receiver->medium();
      // std::tuple<Spectrum, Mask, Float> result =
      //   sample(scene, sampler, ray, medium, aovs+1, active);
      // When doing this, we only render the 'yellow'
      // std::tuple<Spectrum, Mask, Float> result =
      // sample(scene, sampler, ray, medium, aovs+5, active);
      // std::tuple<Spectrum, Mask, Float> result =
      // sample(scene, sampler, ray, medium, aovs+3, active);

      std::pair<Spectrum, Mask> result = sample(scene, sampler, ray, medium,
          aovs + 3, active);

      // Make a new sample routine: one which allows the ray to be changed.
      // In it, the ray's wavelength can change by hitting moving targets. It's
      // time changes with propagation time and length (ray.t) changes with
      // propagation distance. ray.t can be also updated to include phase
      // changes at boundaries, or simply be phase.
      // Maybe I don't even need...

      // ScalarVector2i rd = (std::get<2>(result)/
      //   (sensor->far_clip()-sensor->near_clip()), 0);

      // We'll have to get this from our signal/film prop.
      // Alternatively from our integrator.
      // ScalarVector2i rd = (std::get<2>(result)/
      //   (10-0.1), 0);
      Vector2f rd;
      // range / range interval * nbins
      // rd[0] = std::get<2>(result)/(10.0-0.1)*400;
      // rd[0] = (std::get<2>(result) - (receiver->adc()->centres().x() - receiver->adc()->bandwidth().x()/2)) * receiver->adc()->size().x()/receiver->adc()->bandwidth().x();
      rd[0] = (ray.time - (receiver->adc()->centres().x() - receiver->adc()->bandwidth().x()/2)) * receiver->adc()->size().x()/receiver->adc()->bandwidth().x();

      // Which r bin?
      // I probably want to capture the if signal

      // realrange/(maxrange-minrange)*binsresolution
      // realrange * bins/(rangebandwidth)
      // r * samples/(dr)
      // rd[0] = std::get<2>(result);
      // Put it in the middle
      rd[1] = 0.5f;

      // The point of this function is to get a real value and convert it to a
      // bin idx. This should be similar to the hardware and include things
      // like mixing, filtering etc. Note this is a single return value, the
      // effects of multiple returns comes later.
      //
      // Fuck this, I think it needs to be t/f
      // The internal representation is tf. My only concern is high res.
      // Or is it just always magically delay time?
      // Receiver n_bins
      // Can I skip ahead or do i lose important info? I need to put things
      // into the correct time bins because I need to account for interference.

      // receiver->process_signal

      // Change this back, but allow rays to be modified.
      // std::pair<Spectrum, Mask> result =
      // sample(scene, sampler, ray, medium, aovs+5, active);
      // Gives power, wavelength.

        // std::cout << std::get<0>(result) << std::endl;

        // We currently have a path tracer. ie ray time is at rx.

        // Unconstrained radar path tracer.
        // t_rx is receive time
        // w_rx is transmit/clock frequency now.
        // We don't know what the wavelength will be, because the path time is
        // unknown. Lets send the tx signal NOW, and do everything backward.

        // propagate through scene.
        // when it hits the transmitter, get the time, and tx freq.
        // apply 'mixnfilter'. just do incoherent atm.  to do coherent, we'd
        // need to append each ray landing in a block then do combination sum.

        // Each ray gets multiplied by a wigner layer wrt the path length..or t
        // Do we include bounces? Or does that already happen when we hit the
        // ground? It is simply a -1 after all. If we do path length, that
        // assumes that all rays left with 0 phase. Phase will need to be
        // calculated based on time, we can say 0 phase at time 0, it's all
        // relative anyway....maybe this will even give interference between
        // targets!!

        // I could put phase calc here, but that would force it to be always
        // be wigner. Maybe that's ok considering we're saying that signals are
        // tf anyway. How does this limit my ability to do diff rendering?
        // Let's build the code properly. Keep classic mods so that maybe we
        // can go back.
      result.first = ray_weight * result.first;
      UnpolarizedSpectrum spec_u = depolarize(result.first);

      // std::cout<<spec_u<<std::endl;
      // spec_u should be like radiances at wavelengths

      // This should return a simple power at 1 distance. The accumulated power
      // of the entire sim

      Color3f xyz;
      if constexpr (is_monochromatic_v<Spectrum>) {
          xyz = spec_u.x();
      } else if constexpr (is_rgb_v<Spectrum>) {
          xyz = srgb_to_xyz(spec_u, active);
      } else {
          static_assert(is_spectral_v<Spectrum>);
          // This is in receive sample, so should be only used with adc film
          // xyz = spectrum_to_xyz(spec_u, ray.wavelengths, active);
          xyz = spec_u.x();
          // std::cout<<xyz<<std::endl;
      }
      // if constexpr (is_wigner_v<Spectrum>){}

      // What should I output?
      // I have wavelength

      // float freq = 94E9;
      // float lambda = math::CVac/freq;

      // aovs[0] = xyz.x();

      // aovs[1] = lambda;
      // aovs[1] = select(std::get<1>(result), Float(1.f), Float(0.f));
      // aovs[1] = xyz.y();
      // aovs[2] = xyz.z();
      // aovs[3] = select(std::get<1>(result), Float(1.f), Float(0.f));
      // aovs[4] = 1.f;

      // aovs[0] = xyz.x();
      // aovs[1] = xyz.y();
      // aovs[2] = xyz.z();
      // aovs[3] = select(std::get<1>(result), Float(1.f), Float(0.f));
      // aovs[4] = 1.f;

      // aovs[0] = xyz.x();
      // aovs[1] = select(std::get<1>(result), Float(1.f), Float(0.f));
      // aovs[2] = 1.f;

      aovs[0] = xyz.x();
      aovs[1] = select(result.second, Float(1.f), Float(0.f));
      aovs[2] = 1.f;

      // For interference, what about a weighted sum of phase? Each ray
      // contributes its phase * power. Or If all rays were equal, phase
      // What about tracking phase variance?
      // Or is it another aov?
      // Normalised power*phase.

      // Modify sample so that ray is not const.The sample routine can now
      // change the ray.time and ray.wavelength.
      // aovs become result(basically power), time, wavelength....or can I
      // encode these factors into what they already are?

      // tf_sample = ray.time, ray.wavelengths/frequency

      // block->put(tf_sample, aovs, active);

      // Old:
      // block->put(position_sample, aovs, active);
      // New:
      block->put(rd, aovs, active);
      // block should have tx channel and rx channel. These are aovs...maybe
      // and nah, cause too fast.

      // block->put(position_sample, time_bin (as a fraction of the span), frequency_bin, aovs, active);

      sampler->advance();
}

// MTS_VARIANT std::pair<Spectrum, typename SamplingIntegrator<Float, Spectrum>::Mask>
// SamplingIntegrator<Float, Spectrum>::sample(const Scene * /* scene */,
//                                             Sampler * /* sampler */,
//                                             const RayDifferential3f & /* ray */,
//                                             const Medium * /* medium */,
//                                             Float * /* aovs */,
//                                             Mask /* active */) const {
//     NotImplementedError("sample");
// }

MTS_VARIANT std::pair<Spectrum, typename SamplingIntegrator<Float, Spectrum>::Mask>
SamplingIntegrator<Float, Spectrum>::sample(const Scene * /* scene */,
                                            Sampler * /* sampler */,
                                            RayDifferential3f & /* ray */,
                                            const Medium * /* medium */,
                                            Float * /* aovs */,
                                            Mask /* active */) const {
    NotImplementedError("sample");
}

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
