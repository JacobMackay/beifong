#include <random>
#include <enoki/stl.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/transmitter.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>

// I Guess different integrators have different needs, ie do not NEED transmitters,
// can just do transmitters if required. This may help.

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-path:

Path tracer (:monosp:`path`)
-------------------------------------------

.. pluginparameters::

 * - max_depth
   - |int|
   - Specifies the longest path depth in the generated output image (where -1 corresponds to
     :math:`\infty`). A value of 1 will only render directly visible light sources. 2 will lead
     to single-bounce (direct-only) illumination, and so on. (Default: -1)
 * - rr_depth
   - |int|
   - Specifies the minimum path depth, after which the implementation will start to use the
     *russian roulette* path termination criterion. (Default: 5)
 * - hide_transmitters
   - |bool|
   - Hide directly visible transmitters. (Default: no, i.e. |false|)

This integrator implements a basic path tracer and is a **good default choice**
when there is no strong reason to prefer another method.

To use the path tracer appropriately, it is instructive to know roughly how
it works: its main operation is to trace many light paths using *random walks*
starting from the sensor. A single random walk is shown below, which entails
casting a ray associated with a pixel in the output image and searching for
the first visible intersection. A new direction is then chosen at the intersection,
and the ray-casting step repeats over and over again (until one of several
stopping criteria applies).

.. image:: ../images/integrator_path_figure.png
    :width: 95%
    :align: center

At every intersection, the path tracer tries to create a connection to
the light source in an attempt to find a *complete* path along which
light can flow from the transmitter to the sensor. This of course only works
when there is no occluding object between the intersection and the transmitter.

This directly translates into a category of scenes where
a path tracer can be expected to produce reasonable results: this is the case
when the transmitters are easily "accessible" by the contents of the scene. For instance,
an interior scene that is lit by an area light will be considerably harder
to render when this area light is inside a glass enclosure (which
effectively counts as an occluder).

Like the :ref:`direct <integrator-direct>` plugin, the path tracer internally relies on multiple importance
sampling to combine BSDF and transmitter samples. The main difference in comparison
to the former plugin is that it considers light paths of arbitrary length to compute
both direct and indirect illumination.

.. _sec-path-strictnormals:

.. Commented out for now
.. Strict normals
   --------------

.. Triangle meshes often rely on interpolated shading normals
   to suppress the inherently faceted appearance of the underlying geometry. These
   "fake" normals are not without problems, however. They can lead to paradoxical
   situations where a light ray impinges on an object from a direction that is
   classified as "outside" according to the shading normal, and "inside" according
   to the true geometric normal.

.. The :paramtype:`strict_normals` parameter specifies the intended behavior when such cases arise. The
   default (|false|, i.e. "carry on") gives precedence to information given by the shading normal and
   considers such light paths to be valid. This can theoretically cause light "leaks" through
   boundaries, but it is not much of a problem in practice.

.. When set to |true|, the path tracer detects inconsistencies and ignores these paths. When objects
   are poorly tesselated, this latter option may cause them to lose a significant amount of the
   incident radiation (or, in other words, they will look dark).

.. note:: This integrator does not handle participating media

 */

template <typename Float, typename Spectrum>
class PathTimeFrequencyIntegrator : public MonteCarloIntegrator<Float, Spectrum> {
 public:
    MTS_IMPORT_BASE(MonteCarloIntegrator, m_max_depth, m_rr_depth)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium, Transmitter, TransmitterPtr, BSDF, BSDFPtr)

    explicit PathTimeFrequencyIntegrator(const Properties &props) : Base(props) { }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray_,
                                     const Medium * /* medium */,
                                     Float *aovs ,
                                     Mask active) const override {

        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        RayDifferential3f ray = ray_;

        // Tracks radiance scaling due to index of refraction changes
        Float eta(1.f);

        // MIS weight for intersected transmitters (set by prev. iteration)
        Float emission_weight(1.f);
        // emission_weight *= ray.maxt;
        // ray.maxt = math::Infinity<float>;
        Spectrum throughput(1.f), result(0.f);


        // ==============================================================
        // Took doppler out to test, don't forget to put back
        // ==============================================================

        // ---------------------- First intersection ----------------------

        // In the future this si should probably include doppler, phase and time
        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        Mask valid_ray = si.is_valid();
        TransmitterPtr transmitter = si.transmitter(scene);

        // ray.update_state(-si.t);
        // si.time = ray.time;
        // if(all(si.is_valid())){
        //     std::cout << "Distance: " << si.t << " Phase: " << ray.phase << std::endl;
        // }

        // Doppler
        // if(any_or<true>(neq(si.shape, nullptr))){
        //     const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, valid_ray), 0.f);
        // }



        // This statement sharpens the range
        if(all(si.is_valid())){
        // if(all(active)){
            ray.update_state(-si.t);
            si.time = ray.time;
        }
        for (int depth = 1;; ++depth) {
            // if(all(si.is_valid())){
            //     ray.update_state(-si.t);
            //     si.time = ray.time;
            // }
            // ray.update_state(-si.t);
            // si.time = ray.time;
            // ---------------- Intersection with transmitters ----------------
            // Cancelled for some reason////////double sampling
            if (any_or<true>(neq(transmitter, nullptr))) {

                // if(all(si.is_valid())){
                // // if(all(active)){
                //     ray.update_state(-si.t);
                //     si.time = ray.time;
                // }

                // Something weird is happening here
                // std::cout << "Hit??" << std::endl;

                // Advance the ray travel time ---------------------
                // Have I updated the ray twice here? or no bc int depth 1?
                // ray.update_state(-si.t);
                // si.time = ray.time;
                // =================================================

                // Apply doppler from tx hit -----------------------
                // if(any_or<true>(neq(si.shape, nullptr))){
                //     const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, active), 0.f);
                // }
                // =================================================

                // if(all(si.is_valid())){
                //     ray.update_state(-si.t);
                //     si.time = ray.time;
                // }

                // Evaluate the direct hit illumination -------------
                result[active] +=
                    emission_weight * throughput * transmitter->eval(si, active);

                // if(all(si.is_valid())){
                //     ray.update_state(-si.t);
                //     si.time = ray.time;
                // }

                // these are problems. Currently our only solution without a
                // path tracer, but have to piggyback. si.phase will be some
                // value from the transmitter.
                // const_cast<RayDifferential3f&>(ray_).wavelengths = si.wavelengths;
                // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase;
                // const_cast<RayDifferential3f&>(ray_).phase += ray.phase;
                // ==================================================
            }

            // const_cast<RayDifferential3f&>(ray_).wavelengths = si.wavelengths;
            // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase;
            // const_cast<RayDifferential3f&>(ray_).phase += ray.phase;

            active &= si.is_valid();

            /* Russian roulette: try to keep path weights equal to one,
               while accounting for the solid angle compression at refractive
               index boundaries. Stop with at least some probability to avoid
               getting stuck (e.g. due to total internal reflection) */
            if (depth > m_rr_depth) {
                Float q = min(hmax(depolarize(throughput)) * sqr(eta), .95f);
                active &= sampler->next_1d(active) < q;
                throughput *= rcp(q);
            }

            // Stop if we've exceeded the number of requested bounces, or
            // if there are no more active lanes. Only do this latter check
            // in GPU mode when the number of requested bounces is infinite
            // since it causes a costly synchronization.
            if ((uint32_t) depth >= (uint32_t) m_max_depth ||
                ((!is_cuda_array_v<Float> || m_max_depth < 0) && none(active)))
                break;

            // --------------------- Transmitter sampling ---------------------
            // Ray has hit a scene object. Find the transmitter illumination
            // at this point/time

            BSDFContext ctx;
            BSDFPtr bsdf = si.bsdf(ray);
            Mask active_e =
                active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            // signal_weight = 1.f;

            if (likely(any_or<true>(active_e))) {

                // Advance the ray travel time ---------------------
                // ray.update_state(-si.t);
                // si.time = ray.time;
                // =================================================

                // Apply doppler from bsdf->tx hit -----------------
                if(any_or<true>(neq(si.shape, nullptr))){
                    // // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, active), 0.f);
                    // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, active), 0.f);
                    //
                    // // If the transmitter val != 0, we have a tx hit
                    // // Intersect the BSDF ray against the scene geometry
                    // SurfaceInteraction3f si_tx = scene->ray_intersect(si.spawn_ray(si.to_world(wo)), active);
                    //
                    // if(any_or<true>(neq(si_tx.transmitter(scene, active), nullptr))){
                    //     // const_cast<RayDifferential3f&>(ray_).wavelengths
                    //     //     += select(neq(si_tx.transmitter(scene, active), nullptr),
                    //     //         si_tx.transmitter(scene, active)->doppler(si, active), 0.f);
                    //     const_cast<RayDifferential3f&>(ray_).wavelengths
                    //         += select(neq(si_tx.transmitter(scene, active), nullptr),
                    //             si_tx.transmitter(scene, active)->doppler(si_tx, active), 0.f);
                    //
                    //     std::cout << "Second" << std::endl;
                    //     signal_weight = si_tx.transmitter(scene, active)->eval_signal(ray.time, MTS_C/((ray_).wavelengths[0]*1e-9));
                    //     std::cout << "After" << std::endl;
                    //
                    // }

                }
                // =================================================


                // si is not updated here. I probably need to hit the tx first
                // SurfaceInteraction3f si_tx = scene->ray_intersect(si.spawn_ray(si.to_world(si.to_local(ds.d))), active);

                // Second bounce?
                // Also sharpens. Is this an effect of fmcw?
                // if(all(si.is_valid())){
                // // if(all(active)){
                //     ray.update_state(-si.t);
                //     si.time = ray.time;
                // }

                // So this samples the transmitter directly?
                auto [ds, transmitter_val] = scene->sample_transmitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);
                // ray should be updated now

                // // Being here splits the bounce returns
                // if(all(si.is_valid())){
                //     ray.update_state(-si.t);
                //     si.time = ray.time;
                // }

                // This is good for the main path
                // ....and now seems useless

                // ray.update_state(-ds.dist);
                // si.time = ray.time;

                // const_cast<RayDifferential3f&>(ray_).wavelengths = si.wavelengths;
                // // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase - MTS_P;
                // // // std::cout << ray.phase / math::TwoPi<float> * si.wavelengths[0]*1e-9 << std::endl;
                // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - MTS_P;

                // New
                // ray.update_state(-ds.dist);
                // si.time = ray.time;
                // --
                active_e &= neq(ds.pdf, 0.f);
                // const_cast<RayDifferential3f&>(ray_).wavelengths = si.wavelengths;
                // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase - MTS_P;
                // // std::cout << ray.phase / math::TwoPi<float> * si.wavelengths[0]*1e-9 << std::endl;
                // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - MTS_P;
                // const_cast<RayDifferential3f&>(ray_).time = ray.time;

                // Query the BSDF for that transmitter-sampled direction
                Vector3f wo = si.to_local(ds.d);
                Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // Determine density of sampling that same direction using BSDF
                // sampling
                Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);

                // Evaluate the bsdf->tx illumination -------------
                // When the bsdf IS the transmitter, its bsdfval = 0
                // Does this mis have a huge bias and therefore effect??
                Float mis = select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));
                result[active_e] += mis * throughput * bsdf_val * transmitter_val;
                // ================================================

            }

            // This may not be necessary here, but we'll wait
            // const_cast<RayDifferential3f&>(ray_).time = ray.time;

            // ----------------------- BSDF sampling ----------------------

            // Sample BSDF * cos(theta)
            auto [bs, bsdf_val] = bsdf->sample(ctx, si,
                sampler->next_1d(active), sampler->next_2d(active), active);
            bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

            throughput = throughput * bsdf_val;
            active &= any(neq(depolarize(throughput), 0.f));
            if (none_or<false>(active))
                break;

            eta *= bs.eta;

            // if(all(si.is_valid())){
            //     ray.update_state(-si.t);
            //     si.time = ray.time;
            // }

            // Intersect the BSDF ray against the scene geometry
            ray = si.spawn_ray(si.to_world(bs.wo));
            SurfaceInteraction3f si_bsdf = scene->ray_intersect(ray, active);

            // Ground and generally 1st bounce
            // if(all(si_bsdf.is_valid())){
            if(all(active)){
                ray.update_state(-si_bsdf.t);
                si_bsdf.time = ray.time;
            }

            // // This is good for the main path
            // The ghost signal disappers by changing this
            // is it because we spawn additional rays? That still doesn't add up
            // ray.update_state(-si_bsdf.t);
            // si_bsdf.time = ray.time;

            /* Determine probability of having sampled that same
               direction using transmitter sampling. */
            transmitter = si_bsdf.transmitter(scene, active);
            DirectionSample3f ds(si_bsdf, si);
            ds.object = transmitter;
            // if(all(si_bsdf.is_valid())){
            //     ray.update_state(-si_bsdf.t);
            //     si_bsdf.time = ray.time;
            // }

            // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - MTS_P;
            // const_cast<RayDifferential3f&>(ray_).time = ray.time;

            if (any_or<true>(neq(transmitter, nullptr))) {
                Float transmitter_pdf =
                    select(neq(transmitter, nullptr) &&
                        !has_flag(bs.sampled_type, BSDFFlags::Delta),
                           scene->pdf_transmitter_direction(si, ds),
                           0.f);

                emission_weight = mis_weight(bs.pdf, transmitter_pdf);

                // if(all(si_bsdf.is_valid())){
                //     ray.update_state(-si_bsdf.t);
                //     si_bsdf.time = ray.time;
                // }

                // ray.update_state(-ds.dist);
                // si_bsdf.time = ray.time;
                // if(all(si_bsdf.is_valid())){
                //     ray.update_state(-si_bsdf.t);
                //     si_bsdf.time = ray.time;
                // }

                // const_cast<RayDifferential3f&>(ray_).wavelengths = si_bsdf.wavelengths;
                // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase - MTS_P;
                // // std::cout << ray.phase / math::TwoPi<float> * si.wavelengths[0]*1e-9 << std::endl;
                // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - MTS_P;
                // const_cast<RayDifferential3f&>(ray_).time = ray.time;

                // // New hacks
                // // // If the transmitter val != 0, we have a tx hit
                // // // Intersect the BSDF ray against the scene geometry
                // SurfaceInteraction3f si_tx = scene->ray_intersect(si.spawn_ray(si.to_world(si.to_local(ds.d))), active);
                // // signal_weight = si_tx.transmitter(scene, active)->eval_signal(ray.time, MTS_C/((ray_).wavelengths[0]*1e-9));
                // Spectrum extra_weight = transmitter->eval(si_tx, active);
                // const_cast<RayDifferential3f&>(ray_).wavelengths = si_tx.wavelengths;
                // // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase;
                // const_cast<RayDifferential3f&>(ray_).phase += si_tx.phase;
                // // ray.update_state(-si_bsdf.t);
                // // si_bsdf.time = ray.time;

            }

            si = std::move(si_bsdf);

            // Really not sure about this section, but it definitely does something
            // ray.update_state(-si.t);

            // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - MTS_P;
            // const_cast<RayDifferential3f&>(ray_).time = ray.time;

            // // Another doppler?
            // // if(any_or<true>(neq(si.shape, nullptr))){
            // if(any_or<true>(neq(transmitter, nullptr))){
            //     // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(transmitter, nullptr), transmitter->shape()->doppler(si, active), 0.f);
            //     const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(transmitter, nullptr), transmitter->doppler(si, active), 0.f);
            // }
        }

        if(all(valid_ray)) {
            const_cast<RayDifferential3f&>(ray_).time = ray.time;
            const_cast<RayDifferential3f&>(ray_).wavelengths = si.wavelengths;
            // const_cast<RayDifferential3f&>(ray_).phase += ray.phase - si.phase;
            const_cast<RayDifferential3f&>(ray_).phase += ray.phase;
        }

        // This may not be right
        // const_cast<RayDifferential3f&>(ray_).time = ray.time;
        // result = 1.f;
        return {result, valid_ray};
    }

    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("PathTimeFrequencyIntegrator[\n"
            "  max_depth = %i,\n"
            "  rr_depth = %i\n"
            "]", m_max_depth, m_rr_depth);
    }

    Float mis_weight(Float pdf_a, Float pdf_b) const {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        return select(pdf_a > 0.f, pdf_a / (pdf_a + pdf_b), 0.f);
        // return pdf_a / (pdf_a + pdf_b);
    }

    MTS_DECLARE_CLASS()

// private:
    // ref<Base> m_integrator;
};

MTS_IMPLEMENT_CLASS_VARIANT(PathTimeFrequencyIntegrator, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PathTimeFrequencyIntegrator, "PathTracerTimeFrequencyintegrator");
NAMESPACE_END(mitsuba)
