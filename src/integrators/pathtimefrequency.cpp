#include <random>
#include <enoki/stl.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>

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
 * - hide_emitters
   - |bool|
   - Hide directly visible emitters. (Default: no, i.e. |false|)

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
light can flow from the emitter to the sensor. This of course only works
when there is no occluding object between the intersection and the emitter.

This directly translates into a category of scenes where
a path tracer can be expected to produce reasonable results: this is the case
when the emitters are easily "accessible" by the contents of the scene. For instance,
an interior scene that is lit by an area light will be considerably harder
to render when this area light is inside a glass enclosure (which
effectively counts as an occluder).

Like the :ref:`direct <integrator-direct>` plugin, the path tracer internally relies on multiple importance
sampling to combine BSDF and emitter samples. The main difference in comparison
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
    MTS_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr)

    explicit PathTimeFrequencyIntegrator(const Properties &props) : Base(props) { }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray_,
                                     const Medium * /* medium */,
                                     Float *aovs ,
                                     Mask active) const override {

    // std::pair<Spectrum, Mask> sample(const Scene *scene,
    //                                  Sampler *sampler,
    //                                  RayDifferential3f &ray_,
    //                                  const Medium * /* medium */,
    //                                  Float *aovs ,
    //                                  Mask active) const override {

    // std::pair<Spectrum, Mask> sample(const Scene *scene,
    //                                  Sampler *sampler,
    //                                  RayDifferential3f *ray_,
    //                                  const Medium * /* medium */,
    //                                  Float *aovs ,
    //                                  Mask active) const override {

    // std::tuple<Spectrum, Mask, Float> sample(const Scene *scene,
    //                                  Sampler *sampler,
    //                                  const RayDifferential3f &ray_,
    //                                  const Medium * /* medium */,
    //                                  Float * /* aovs*/ ,
    //                                  Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        // Possibility to change ray from const so that when function returns
        // we keep ray.

        // RayDifferential3f ray = *ray_;
        RayDifferential3f ray = ray_;

        // std::cout << ray.time << std:: endl;

        // Tracks radiance scaling due to index of refraction changes
        Float eta(1.f);

        // MIS weight for intersected emitters (set by prev. iteration)
        Float emission_weight(1.f);

        Spectrum throughput(1.f), result(0.f);
        // Float pathlength(0.f);

        // ---------------------- First intersection ----------------------

        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        Mask valid_ray = si.is_valid();
        EmitterPtr emitter = si.emitter(scene);

        // pathlength = select(si.is_valid(), si.t, 0.f);
        // pathlength = select(si.is_valid(), si.t, math::Infinity<Float>);
        // pathlength = select(valid_ray, si.t, math::Infinity<Float>);

        // pathlength += select(valid_ray, si.t, 0.f);
        // const_cast<RayDifferential3f&>(ray_).time -= select(valid_ray, si.t / math::CVac<float>, 0.f);
        ray.time -= select(valid_ray, si.t / math::CVac<float>, 0.f);
        // ray_time -= select(valid_ray, si.t/C_AIR, 0.f);

        // pathlength += select(si.is_valid(), si.t, 0.f);

        // pathlength += si.t;
        // pathlength[active] += si.t;

        for (int depth = 1;; ++depth) {
            // ---------------- Intersection with emitters ----------------

            if (any_or<true>(neq(emitter, nullptr))) {
                result[active] +=
                    emission_weight * throughput * emitter->eval(si, active);

                // pathlength += select(si.is_valid(), si.t, 0.f);
                // const_cast<RayDifferential3f&>(ray_).time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
                ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
                // ray_time -= select(valid_ray, si.t/C_AIR, 0.f);
                // Lets say the transmitter is moving at 10 m/s
                // ray.wavelengths += (100.f/math::CVac<float>)*ray.wavelengths;
                // delf = delv/c * f0
                // c/dell = delv /c * c/l0
                // del l = c / delv * l0
            }

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

            // --------------------- Emitter sampling ---------------------

            BSDFContext ctx;
            BSDFPtr bsdf = si.bsdf(ray);
            Mask active_e =
                active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            if (likely(any_or<true>(active_e))) {
                auto [ds, emitter_val] = scene->sample_emitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);
                active_e &= neq(ds.pdf, 0.f);

                // Query the BSDF for that emitter-sampled direction
                Vector3f wo = si.to_local(ds.d);
                Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // Determine density of sampling that same direction using BSDF
                // sampling
                Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);

                Float mis = select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));
                result[active_e] += mis * throughput * bsdf_val * emitter_val;

                // pathlength += select(si.is_valid(), si.t, 0.f);
                // const_cast<RayDifferential3f&>(ray_).time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
                ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
                // ray_time -= select(valid_ray, si.t/C_AIR, 0.f);
                // output_wavelength = ray_wavelength - tx.get_wavelength
                // put into mixed equation. get result.

                // New class, like film But for tx.
                // Transmitter has signal.
                // Different types of signals with diff parameters
                // But basically, given a time, return a wavelength/amplitude.
                // Ideally it'd be a wdf.
                // Can have a member, is delta. If it is, get exact result.

                // The returned ray should have its time as time, and
                // wavelength as beat wavelength.
                // We do the AA filter later in rx.
                // For a given time and frequency, the tx should have an amp.
                // The signal should be analytic or texture.

                // There's signal bandwidth, but also tx/rx bandwidth.

                // Hardcode for now.
                // Given a t and f, return the power; alternatively, given an t

                Wavelength f_tx(0.f);
                // f0 = fc - bandwidth/2;
                // f1 = fc + bandwidth/2;
                //
                // t1 = rise;
                // t2 = t1 + hold;
                // t3 = t2 + fall;
                // t4 = t3 + wait;
                // tn = t % t4;

                // Given a time, return the transmit frequency.
                // This means that the transmitter needs to have a signal.

                Wavelength f0 = 94e9 - 6e9/2;
                Wavelength f1 = 94e9 + 6e9/2;

                Float t1 = 240e-6;
                Float t2 = t1 + 10e-6;
                Float t3 = t2 + 240e-6;
                Float t4 = t3 + 10e-6;

                Float tn = math::modulo(ray.time, t4);
                // Float tn = ray.time;

                // std::cout << tn << ", " << ray.time << std::endl;

                if (all(tn < t1)) {
                    // f = 2*((f1 - f0)/(t2 - t1))*tn + f0;
                    // f = 2*((6e9)/(240e-6))*tn + f0;
                    f_tx = ((6e9)/(240e-6))*tn + f0;
                } else if (all(tn < t2)) {
                    f_tx = f1;
                } else if (all(tn < t3)){
                    // f = 2*((f0 - f1)/(t3 - t2))*tn + f1;
                    // f = 2*((-6e9)/(240e-6))*tn + f1;
                    f_tx = ((-6e9)/(240e-6))*(tn-t2) + f1;
                } else {
                    f_tx = f0;
                }

                // const_cast<RayDifferential3f&>(ray_).wavelengths = math::CVac<double>/f_tx *1e9;

                // ray.wavelengths += (100.f/math::CVac<float>)*ray.wavelengths;
                // Imagine the tx is moving at 20 m/s
                f_tx += 20.f/math::CVac<float> * f_tx;

                // std::cout << ray.time << " " << tn << std::endl;

                // std::cout << "tx_t: " << tn * 1e6 << " rx_t: " << ray_.time * 1e6
                //     << " tx_f: " << f_tx[0]*1e-9 << " rx_f: " << math::CVac<double>/(ray.wavelengths[0]*1e-9)*1e-9
                //     << " d_f: " << (f_tx[0] - math::CVac<double>/(ray.wavelengths[0]*1e-9))
                //     << " tx_λ: " << ray_.wavelengths[0]*1e-9*1e3 << " rx_λ: " << ray.wavelengths[0]*1e-9*1e3 << std::endl;

                // const_cast<RayDifferential3f&>(ray_).wavelengths = 0.5*(abs(f - math::CVac<double>/(ray.wavelengths*1e-9)))*math::InvTwoPi<double>;
                // const_cast<RayDifferential3f&>(ray_).wavelengths = (abs(f - math::CVac<double>/(ray.wavelengths*1e-9)))*math::InvTwoPi<double>;
                // const_cast<RayDifferential3f&>(ray_).wavelengths = (abs(f_tx - math::CVac<double>/(ray.wavelengths*1e-9)));
                // const_cast<RayDifferential3f&>(ray_).wavelengths -= ray.wavelengths;

                // const_cast<RayDifferential3f&>(ray_).wavelengths = (math::CVac<double>/ray.wavelengths - math::CVac<double>/ray_.wavelengths)*1e9;
                const_cast<RayDifferential3f&>(ray_).wavelengths = abs(f_tx[0] - math::CVac<double>/(ray.wavelengths[0]*1e-9));

                // std::cout << ray_.wavelengths[0] << std::endl;

            }

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

            // Intersect the BSDF ray against the scene geometry
            ray = si.spawn_ray(si.to_world(bs.wo));
            SurfaceInteraction3f si_bsdf = scene->ray_intersect(ray, active);

            /* Determine probability of having sampled that same
               direction using emitter sampling. */
            emitter = si_bsdf.emitter(scene, active);
            DirectionSample3f ds(si_bsdf, si);
            ds.object = emitter;

            if (any_or<true>(neq(emitter, nullptr))) {
                Float emitter_pdf =
                    select(neq(emitter, nullptr) &&
                        !has_flag(bs.sampled_type, BSDFFlags::Delta),
                           scene->pdf_emitter_direction(si, ds),
                           0.f);

                emission_weight = mis_weight(bs.pdf, emitter_pdf);
            }

            si = std::move(si_bsdf);
            // pathlength += select(si.is_valid(), si.t, 0.f);
            // const_cast<ray&>(ray_).time -= select(si.is_valid(), si.t/ math::CVac<float>, 0.f);
            // const_cast<RayDifferential3f&>(ray_).time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
            ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
            // ray_time -= select(valid_ray, si.t/C_AIR, 0.f);

            // pathlength += select(si.is_valid(), si.t, math::Infinity<Float>);
            // pathlength += select(active, si.t, math::Infinity<Float>);

            // pathlength += select(active, si.t, 0.f);
            // pathlength += select(si.is_valid(), si.t, 0.f);

            // pathlength += select(active_e, si.t, math::Infinity<Float>);
            // pathlength += select(active_e, si.t, 0.f);
            // pathlength += si.t;
            // pathlength[active] += si.t;
        }

        // std::cout << pathlength << std::endl;

        // return { result, valid_ray, pathlength};

        // ray_.time = pathlength / math::CVac<float>;

        // Float f_tx(0.f);
        // // f0 = fc - bandwidth/2;
        // // f1 = fc + bandwidth/2;
        // //
        // // t1 = rise;
        // // t2 = t1 + hold;
        // // t3 = t2 + fall;
        // // t4 = t3 + wait;
        // // tn = t % t4;
        //
        // // Given a time, return the transmit frequency.
        // // This means that the transmitter needs to have a signal.
        //
        // Float f0 = 94e9 - 6e9/2;
        // Float f1 = 94e9 + 6e9/2;
        //
        // Float t1 = 240e-6;
        // Float t2 = t1 + 10e-6;
        // Float t3 = t2 + 240e-6;
        // Float t4 = t3 + 10e-6;
        // // Float tn = math::modulo(ray_.time, t4);
        //
        // // Float tn = math::modulo(t4, ray.time);
        // // Float tn = math::modulo(ray.time, t4);
        // Float tn = ray.time;
        //
        // // Float tn = ray.time % t4;
        //
        // // std::cout << tn << ", " << ray.time << std::endl;
        //
        // if (all(tn < t1)) {
        //     // f = 2*((f1 - f0)/(t2 - t1))*tn + f0;
        //     // f = 2*((6e9)/(240e-6))*tn + f0;
        //     f_tx = ((6e9)/(240e-6))*tn + f0;
        // } else if (all(tn < t2)) {
        //     f_tx = f1;
        // } else if (all(tn < t3)){
        //     // f = 2*((f0 - f1)/(t3 - t2))*tn + f1;
        //     // f = 2*((-6e9)/(240e-6))*tn + f1;
        //     f_tx = ((-6e9)/(240e-6))*tn + f1;
        // } else {
        //     f_tx = f0;
        // }
        //
        // std::cout << tn * 1e6 << " " << f_tx*1e-9 << " " << math::CVac<double>/(ray.wavelengths[0]*1e-9)*1e-9  << " " << ray.wavelengths[0]*1e-9 << std::endl;
        //
        // // const_cast<RayDifferential3f&>(ray_).wavelengths = 0.5*(abs(f - math::CVac<double>/(ray.wavelengths*1e-9)))*math::InvTwoPi<double>;
        // // const_cast<RayDifferential3f&>(ray_).wavelengths = (abs(f - math::CVac<double>/(ray.wavelengths*1e-9)))*math::InvTwoPi<double>;
        // const_cast<RayDifferential3f&>(ray_).wavelengths = (abs(f_tx - math::CVac<double>/(ray.wavelengths*1e-9)));

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
