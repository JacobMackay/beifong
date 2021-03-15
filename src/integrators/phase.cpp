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
class PhaseIntegrator final : public SamplingIntegrator<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(SamplingIntegrator)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium)

    PhaseIntegrator(const Properties &props) : Base(props) {
        // if constexpr (!is_polarized_v<Spectrum>)
        //     Throw("This integrator should only be used in polarized mode!");
        for (auto &kv : props.objects()) {
            Base *integrator = dynamic_cast<Base *>(kv.second.get());
            if (!integrator)
                Throw("Child objects must be of type 'SamplingIntegrator'!");
            if (m_integrator)
                Throw("More than one sub-integrator specified!");
            m_integrator = integrator;
        }

        if (!m_integrator)
            Throw("Must specify a sub-integrator!");

        m_bins = props.int_("bins", 1);
        m_bin_width = math::TwoPi<float>/m_bins;
        for (int k = 0; k < m_bins; ++k) {
            m_bin_centres.push_back(m_bin_width*(k + 0.5));
        }
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler * sampler,
                                     const RayDifferential3f &ray,
                                     const Medium *medium,
                                     Float *aovs,
                                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        auto result = m_integrator->sample(scene, sampler, ray, medium, aovs, active);

        UnpolarizedSpectrum spec_u = depolarize(result.first);
        Color3f xyz;
        if constexpr (is_monochromatic_v<Spectrum>) {
            xyz = spec_u.x();
        } else if constexpr (is_rgb_v<Spectrum>) {
            xyz = srgb_to_xyz(spec_u, active);
        } else {
            static_assert(is_spectral_v<Spectrum>);
            xyz = spec_u.x();
        }

// std::cout << math::fmodulo(ray.phase, math::TwoPi<float>) << " " << std::fmod(hsum(ray.phase), math::TwoPi<float>) << std::endl;

        // Float phase = math::fmodulo(ray.phase, math::TwoPi<float>);
        Float phase = std::fmod(hsum(ray.phase), math::TwoPi<float>);
        phase += select(all(phase<0.f), math::TwoPi<float>, 0.f);
        // if (phase != 0){
        for (int k = 0; k < m_bins; ++k){
            *aovs++ = select(
                math::rect((phase - m_bin_centres[k])/m_bin_width) > 0.f,
                xyz.x(),
                0.f);
        }
        // }

        return result;
    }

    std::vector<std::string> aov_names() const override {
        std::vector<std::string> result = m_integrator->aov_names();
        for (int i = 0; i < m_bins; ++i)
            result.insert(result.begin() + i, "S" + std::to_string(i) + "." + "Y");
        return result;
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("integrator", m_integrator.get());
    }

    MTS_DECLARE_CLASS()
private:
    ref<Base> m_integrator;
    int m_bins;
    float m_bin_width;
    std::vector<float> m_bin_centres;
};

MTS_IMPLEMENT_CLASS_VARIANT(PhaseIntegrator, SamplingIntegrator)
MTS_EXPORT_PLUGIN(PhaseIntegrator, "Phase integrator");
NAMESPACE_END(mitsuba)
