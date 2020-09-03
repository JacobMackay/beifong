#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-stokes:

Stokes vector integrator (:monosp:`stokes`)
-----------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - :paramtype:`integrator`
   - Sub-integrator (only one can be specified) which will be sampled along the Stokes
     integrator. In polarized rendering modes, its output Stokes vector is written
     into distinct images.

This integrator returns a multi-channel image describing the complete measured
polarization state at the sensor, represented as a Stokes vector :math:`\mathbf{s}`.

Here we show an example monochrome output in a scene with two dielectric and one
conductive sphere that all affect the polarization state of the
(initially unpolarized) light.

The first entry corresponds to usual radiance, whereas the remaining three entries
describe the polarization of light shown as false color images (green: positive, red: negative).

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/integrator_stokes_cbox.jpg
   :caption: ":math:`\mathbf{s}_0`": radiance
.. subfigure:: ../../resources/data/docs/images/render/integrator_stokes_cbox_s1.jpg
   :caption: ":math:`\mathbf{s}_1`": horizontal vs. vertical polarization
.. subfigure:: ../../resources/data/docs/images/render/integrator_stokes_cbox_s2.jpg
   :caption: ":math:`\mathbf{s}_2`": positive vs. negative diagonal polarization
.. subfigure:: ../../resources/data/docs/images/render/integrator_stokes_cbox_s3.jpg
   :caption: ":math:`\mathbf{s}_3`": right vs. left circular polarization
.. subfigend::
   :label: fig-stokes

In the following example, a normal path tracer is nested inside the Stokes vector
integrator:

.. code-block:: xml

    <integrator type="stokes">
        <integrator type="path">
            <!-- path tracer parameters -->
        </integrator>
    </integrator>


 */

template <typename Float, typename Spectrum>
class RangeIntegrator final : public SamplingIntegrator<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(SamplingIntegrator)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium)

    RangeIntegrator(const Properties &props) : Base(props) {
        // if constexpr (!is_polarized_v<Spectrum>)
        //     Throw("This integrator should only be used in polarized mode!");
        for (auto &kv : props.objects()) {
            Base *integrator = dynamic_cast<Base *>(kv.second.get());
            if (!integrator)
                Throw("Child objects must be of type 'SamplingIntegrator'!");
            if (m_integrator)
                Throw("More than one sub-integrator specified!");
            m_integrator = integrator;

            m_dr = props.float_("dr", -1.f);
            m_bins = props.int_("bins", -1);
        }

        if (!m_integrator)
            Throw("Must specify a sub-integrator!");
    }

    // std::pair<Spectrum, Mask> sample(const Scene *scene,
    //                                  Sampler * sampler,
    //                                  const RayDifferential3f &ray,
    //                                  const Medium *medium,
    //                                  Float *aovs,
    //                                  Mask active) const override {
    //     MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);
    std::tuple<Spectrum, Mask, Float> sample(const Scene *scene,
                                     Sampler * sampler,
                                     const RayDifferential3f &ray,
                                     const Medium *medium,
                                     Float *aovs,
                                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        // auto result = m_integrator->sample(scene, sampler, ray, medium, aovs + 12, active);
        auto result = m_integrator->sample(scene, sampler, ray, medium, aovs, active);
        // Don't know if/why we need to specify number of aovs.
        // Lets say 10 bins for now, get it from input later

        // ------------- from image block ----------------------------------
        // Convert to pixel coordinates within the image block
        // pos_ in fractional pixel coordintes? Bins already defined by pixels
        // time should be scaled into bins
        // Point2f pos = pos_ - (m_offset - m_border_size + .5f);
        // Point2u lo = ceil2int<Point2i>(pos - .5f);
        // UInt32 offset = m_channel_count * (lo.y() * size.x() + lo.x());
        //
        // Mask enabled = active && all(lo >= 0u && lo < size);
        // ENOKI_NOUNROLL for (uint32_t k = 0; k < m_channel_count; ++k)
        //     scatter_add(m_data, value[k], offset + k, enabled);
        //---------------------------------------------------------------


        //
        // using Int32 = int32_array_t<Value>;
        // Int32 index = arange<Int32>(n * size.y());  // 1:len

        // Float dt = 200.0e-6;
        // Float dt = 0.5e-9;
        // Float dr = 0.1;

        // std::cout << m_bins <<std::endl;

        // Color3f xyz;
        // if constexpr (is_monochromatic_v<Spectrum>) {
        //     xyz = spec_u.x();
        // } else if constexpr (is_rgb_v<Spectrum>) {
        //     xyz = srgb_to_xyz(spec_u, active);
        // } else {
        //     static_assert(is_spectral_v<Spectrum>);
        //     xyz = spectrum_to_xyz(spec_u, ray.wavelengths, active);
        // }

        auto const &ranges = std::get<2>(result);
        auto const &spec_u = depolarize(std::get<0>(result));
        // auto const &rads = depolarize(std::get<0>(result));
        // auto const &rads = std::get<0>(result);
        for (int i = 0; i < m_bins; ++i){

            Point1f lo = (Float)i *m_dr;
            Point1f hi = (Float)i *m_dr + m_dr;

            Color3f xyz;
            if constexpr (is_monochromatic_v<Spectrum>) {
                // xyz = spec_u.x();
                xyz = select(all(ranges>=lo && ranges<hi), spec_u.x(), 0.f);
                *aovs++ = xyz.x(); *aovs++ = xyz.y(); *aovs++ = xyz.z();
            } else if constexpr (is_rgb_v<Spectrum>) {
                xyz = select(all(ranges>=lo && ranges<hi), srgb_to_xyz(spec_u, active), 0.f);
                *aovs++ = xyz.x(); *aovs++ = xyz.y(); *aovs++ = xyz.z();
            } else {
                static_assert(is_spectral_v<Spectrum>);
                // xyz = select(all(times>=lo && times<hi), spectrum_to_xyz(spec_u, ray.wavelengths, active), 0.f);
                // *aovs++ = xyz.x(); *aovs++ = xyz.y(); *aovs++ = xyz.z();
                xyz = select(all(ranges>=lo && ranges<hi), spec_u.x(), 0.f);
                *aovs++ = xyz.x(); *aovs++ = xyz.y(); *aovs++ = xyz.z();
            }

            // Color3f rgb;
            // Point3f rgb;
            // Point1f lo = Float i * dt;
            // Point1f hi = Float i * dt + dt;

            // rgb = rads[active && all(ray.time>=lo && ray.time<hi)];
            // rgb = rads[all(ray.time>=lo && ray.time<hi)];
            // rgb = select(all(ray.time>=lo && ray.time<hi), rads, 0.f);
            // rgb = select(all(times>=lo && times<hi), rads, 0.f)
            // std::cout<<ray.time<<std::endl;
            // rgb = rads
            // *aovs++ = rgb.r(); *aovs++ = rgb.g(); *aovs++ = rgb.b();
        }

        // if constexpr (is_polarized_v<Spectrum>) {
        //     auto const &stokes = result.first.coeff(0);
        //     for (int i = 0; i < 4; ++i) {
        //         Color3f rgb;
        //         if constexpr (is_monochromatic_v<Spectrum>) {
        //             rgb = stokes[i].x();
        //         } else if constexpr (is_rgb_v<Spectrum>) {
        //             rgb = stokes[i];
        //         } else {
        //             static_assert(is_spectral_v<Spectrum>);
        //             /// Note: this assumes that sensor used sample_rgb_spectrum() to generate 'ray.wavelengths'
        //             auto pdf = pdf_rgb_spectrum(ray.wavelengths);
        //             UnpolarizedSpectrum spec = stokes[i] * select(neq(pdf, 0.f), rcp(pdf), 0.f);
        //             rgb = xyz_to_srgb(spectrum_to_xyz(spec, ray.wavelengths, active));
        //         }
        //
        //         *aovs++ = rgb.r(); *aovs++ = rgb.g(); *aovs++ = rgb.b();
        //     }
        // }

        return result;
    }

    // std::vector<std::string> aov_names() const override {
    //     std::vector<std::string> result = m_integrator->aov_names();
    //     for (int i = 0; i < 4; ++i)
    //         for (int j = 0; j < 3; ++j)
    //             result.insert(result.begin() + 3*i + j, "S" + std::to_string(i) + "." + ("RGB"[j]));
    //     return result;
    // }

    std::vector<std::string> aov_names() const override {
        std::vector<std::string> result = m_integrator->aov_names();
        for (int i = 0; i < 50; ++i)
            for (int j = 0; j < 3; ++j)
                result.insert(result.begin() + 3*i + j, "S" + std::to_string(i) + "." + ("RGB"[j]));
        return result;
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("integrator", m_integrator.get());
    }

    MTS_DECLARE_CLASS()
private:
    ref<Base> m_integrator;
    /// Jacobs mess
    float m_dr;
    float m_bins;
};

MTS_IMPLEMENT_CLASS_VARIANT(RangeIntegrator, SamplingIntegrator)
MTS_EXPORT_PLUGIN(RangeIntegrator, "Range integrator");
NAMESPACE_END(mitsuba)
