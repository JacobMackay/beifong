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

        // Tracks radiance scaling due to index of refraction changes
        Float eta(1.f);

        // MIS weight for intersected transmitters (set by prev. iteration)
        Float emission_weight(1.f);

        Spectrum throughput(1.f), result(0.f);

        // ==============================================================
        // Took doppler out to test, don't forget to put back
        // ==============================================================

        // ---------------------- First intersection ----------------------

        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        Mask valid_ray = si.is_valid();
        TransmitterPtr transmitter = si.transmitter(scene);

        ray.time -= select(valid_ray, si.t / math::CVac<float>, 0.f);
        if(any_or<true>(neq(si.shape, nullptr))){
            // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, valid_ray), 0.f);
        }

        for (int depth = 1;; ++depth) {
            // ---------------- Intersection with transmitters ----------------

            if (any_or<true>(neq(transmitter, nullptr))) {
                // I can sample from the transmitter here: sample_ray with random numbers.
                // If it has a delta distribution, we can clamp it to what it should be.
                // usually takes 0-1 argument, what should i do when i have a real number for time?
                result[active] +=
                    emission_weight * throughput * transmitter->eval(si, active);

                // We have the time t, and beat frequency f.
                // We can find the tx freq when we have a time.
                // We can now find the frequency which gives us the required beat frequency.
                // With the tx freq we can find a corresponding time
                // With the difference in time we can have a range/delay.
                // Also, beat freq corresponds to range

                // transmitter eval should call the signal.si will have time and wavelength.
                // the receiver
                // to this function came in a ray with time t and wavelengths λ
                //
                // Also, the antenna would have a 'filter' which smears the signal,
                // creates delay and amplitude change
                // So call this function, I need a si.t, si.λ
                // Alternatively: time is real arrival t, λ is tx λ-> this implies a travel time
                // the antennas should have a function which converts wavelength to frequency...but later

                // pathlength += select(si.is_valid(), si.t, 0.f);
                // const_cast<RayDifferential3f&>(ray_).time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);

                // std::cout << "tx" << std::endl;

                ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
                if(any_or<true>(neq(si.shape, nullptr))){
                    // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(si.shape, nullptr), si.shape->doppler(si, active), 0.f);
                }

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

            // --------------------- Transmitter sampling ---------------------

            BSDFContext ctx;
            BSDFPtr bsdf = si.bsdf(ray);
            Mask active_e =
                active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            if (likely(any_or<true>(active_e))) {
                auto [ds, transmitter_val] = scene->sample_transmitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);
                active_e &= neq(ds.pdf, 0.f);

                // Query the BSDF for that transmitter-sampled direction
                Vector3f wo = si.to_local(ds.d);
                Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // Determine density of sampling that same direction using BSDF
                // sampling
                Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);

                Float mis = select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));
                result[active_e] += mis * throughput * bsdf_val * transmitter_val;

                ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
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
                    // }

                }


                // std::cout << si.shape->doppler(si) << std::endl;

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

                // Wavelength f_tx(0.f);
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
                // Wavelength f0 = 94e9 - 6e9/2;
                // Wavelength f1 = 94e9 + 6e9/2;
                //
                // Float t1 = 240e-6;
                // Float t2 = t1 + 10e-6;
                // Float t3 = t2 + 240e-6;
                // Float t4 = t3 + 10e-6;
                //
                // // std::cout << "pt in" << std::endl;
                // Float tn = math::fmodulo(ray.time, t4);
                // // Float tn = math::modulo(ray.time, t4);
                // // std::cout << "pt out" << std::endl;
                // // Float tn = ray.time;
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
                //     f_tx = ((-6e9)/(240e-6))*(tn-t2) + f1;
                // } else {
                //     f_tx = f0;
                // }

                // Wavelength doppler() const {
                //
                //     // if (shape) {
                //
                //         // Transform4f vel = shape->velocity;
                //         // ScalarTransform4f vel = shape->velocity();
                //         ScalarTransform4f vel = ScalarTransform4f(shape->velocity());
                //         // DirectionSample3f ws = m_shape->sample_wigner(ds, wavelength, active);
                //
                //         // return dot(wi, shape.velocity*Point3f(to_local(p)))/math::CVac<float>*wavelengths;
                //         return dot(wi, vel*Point3f(to_local(p)))/math::CVac<float>*wavelengths;
                //         // return dot(wi, shape.velocity*Point3f(to_local(p)))/math::CVac<float>*wavelengths;
                //
                //
                //         // return dot(si.wi, m_velocity*Point3f(si.to_local(si.p)))/math::CVac<float>*si.wavelengths;
                //
                //     // } else {
                //     //     return 0;
                //     // }
                //
                //
                //
                // }

                // ScalarTransform4f vel = si.instance(scene)->to_world();

                // Wavelength dop = si.instance->doppler(si);
                // Should I do active checks or non-nullptr checks?
                // Wavelength dop = si.shape->doppler(si);
                // const_cast<RayDifferential3f&>(ray_).wavelengths += si.shape->doppler(si);

                // std::cout << dop << std::endl;
                //
                // f_tx += math::CVac<float>/dop;


                // ScalarTransform4f vel = ScalarTransform4f(si.shape->velocity());
                // ScalarTransform4f vel = si.instance->to_world();
                // ScalarTransform4f vel = si.instance.m_to_world;

                // ShapePtr *SP = si.instance;
                // ScalarFloat SA = si.instance->surface_area();
                // ScalarTransform4f vel = si.instance->velocity();

                // const_cast<RayDifferential3f&>(ray_).wavelengths = math::CVac<double>/f_tx *1e9;

                // ray.wavelengths += (100.f/math::CVac<float>)*ray.wavelengths;
                // Imagine the tx is moving at 20 m/s
                // f_tx += -100.f/math::CVac<float> * f_tx;
                // f_tx += math::CVac<float>/si.shape->eval_doppler(si);
                // f_tx += math::CVac<float>/si.doppler();

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
                // const_cast<RayDifferential3f&>(ray_).wavelengths = abs(f_tx[0] - math::CVac<double>/(ray.wavelengths[0]*1e-9));

                // const_cast<RayDifferential3f&>(ray_).wavelengths = abs(f_tx[0] - math::CVac<double>/(ray.wavelengths[0]*1e-9));

                // const_cast<RayDifferential3f&>(ray_).wavelengths = abs(math::CVac<double>/f_tx*1e9 - ray.wavelengths);
                // const_cast<RayDifferential3f&>(ray_).wavelengths = ray.wavelengths;
                const_cast<RayDifferential3f&>(ray_).time = ray.time;


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
               direction using transmitter sampling. */
            transmitter = si_bsdf.transmitter(scene, active);
            DirectionSample3f ds(si_bsdf, si);
            ds.object = transmitter;

            if (any_or<true>(neq(transmitter, nullptr))) {
                Float transmitter_pdf =
                    select(neq(transmitter, nullptr) &&
                        !has_flag(bs.sampled_type, BSDFFlags::Delta),
                           scene->pdf_transmitter_direction(si, ds),
                           0.f);

                emission_weight = mis_weight(bs.pdf, transmitter_pdf);
            }

            si = std::move(si_bsdf);

            // Really not sure about this section, but it definitely does something
            ray.time -= select(si.is_valid(), si.t / math::CVac<float>, 0.f);
            // if(any_or<true>(neq(si.shape, nullptr))){
            if(any_or<true>(neq(transmitter, nullptr))){
                // // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(transmitter, nullptr), transmitter->shape()->doppler(si, active), 0.f);
                // const_cast<RayDifferential3f&>(ray_).wavelengths += select(neq(transmitter, nullptr), transmitter->doppler(si, active), 0.f);



            }
        }

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
