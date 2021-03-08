#pragma once

#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/render/endpoint.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/adc.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/records.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER Receiver : public Endpoint<Float, Spectrum> {
public:
    MTS_IMPORT_TYPES(ADC, Sampler)
    MTS_IMPORT_BASE(Endpoint, sample_ray, m_needs_sample_3)

    // =============================================================
    //! @{ \name Sensor-specific sampling functions
    // =============================================================

    /**
     * \brief Importance sample a ray differential proportional to the sensor's
     * sensitivity profile.
     *
     * The sensor profile is a six-dimensional quantity that depends on time,
     * wavelength, surface position, and direction. This function takes a given
     * time value and five uniformly distributed samples on the interval [0, 1]
     * and warps them so that the returned ray the profile. Any
     * discrepancies between ideal and actual sampled profile are absorbed into
     * a spectral importance weight that is returned along with the ray.
     *
     * In contrast to \ref Endpoint::sample_ray(), this function returns
     * differentials with respect to the X and Y axis in screen space.
     *
     * \param time
     *    The scene time associated with the ray_differential to be sampled
     *
     * \param sample1
     *     A uniformly distributed 1D value that is used to sample the spectral
     *     dimension of the sensitivity profile.
     *
     * \param sample2
     *    This argument corresponds to the sample position in fractional pixel
     *    coordinates relative to the crop window of the underlying film.
     *
     * \param sample3
     *    A uniformly distributed sample on the domain <tt>[0,1]^2</tt>. This
     *    argument determines the position on the aperture of the sensor. This
     *    argument is ignored if <tt>needs_sample_3() == false</tt>.
     *
     * \return
     *    The sampled ray differential and (potentially spectrally varying)
     *    importance weights. The latter account for the difference between the
     *    sensor profile and the actual used sampling density function.
     */
    virtual std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float sample1,
                            const Point2f &sample2, const Point2f &sample3,
                            Mask active = true) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Additional query functions
    // =============================================================

    /// Return the time value of the shutter opening event
    ScalarFloat adc_sampling_start() const { return m_adc_sampling_start; }

    /// Return the length, for which the shutter remains open
    ScalarFloat adc_sampling_time() const { return m_adc_sampling_time; }

    /// Does the sampling technique require a sample for the aperture position?
    bool needs_aperture_sample() const { return m_needs_sample_3; }

    /// Return the \ref ADC instance associated with this sensor
    ADC *adc() { return m_adc; }

    /// Return the \ref ADC instance associated with this sensor (const)
    const ADC *adc() const { return m_adc.get(); }

    // virtual Wavelength doppler(SurfaceInteraction3f si, Mask active) const;

    /**
     * \brief Return the sensor's sample generator
     *
     * This is the \a root sampler, which will later be cloned a
     * number of times to provide each participating worker thread
     * with its own instance (see \ref Scene::sampler()).
     * Therefore, this sampler should never be used for anything
     * except creating clones.
     */
    Sampler *sampler() { return m_sampler; }

    /**
     * \brief Return the sensor's sampler (const version).
     *
     * This is the \a root sampler, which will later be cloned a
     * number of times to provide each participating worker thread
     * with its own instance (see \ref Scene::sampler()).
     * Therefore, this sampler should never be used for anything
     * except creating clones.
     */
    const Sampler *sampler() const { return m_sampler.get(); }

    std::string receive_type() const { return m_receive_type; }

    //! @}
    // =============================================================

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("adc_sampling_start", m_adc_sampling_start);
        callback->put_parameter("adc_sampling_time", m_adc_sampling_time);
        callback->put_object("adc", m_adc.get());
        callback->put_object("sampler", m_sampler.get());
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {})
        override {
        m_resolution = ScalarVector2f(m_adc->window_size());
    }

    std::pair<wavelength_t<Spectrum>, Spectrum> sample_frequency(Float time, Float sample) const;
    Spectrum eval_signal(Float time, wavelength_t<Spectrum> frequency) const;

    ENOKI_CALL_SUPPORT_FRIEND()
    MTS_DECLARE_CLASS()
protected:
    Receiver(const Properties &props);

    virtual ~Receiver();

protected:
    ref<ADC> m_adc;
    ref<Sampler> m_sampler;
    ScalarVector2f m_resolution;
    ScalarFloat m_adc_sampling_start;
    // Ignoring the tx, when does the receiver start measuring and for how long
    // This should correspond to the adc collection window. So if we're doing
    // high res, long range it will be the whole thing. If we're doing high res
    // , particular region we have this as well.
    ScalarFloat m_adc_sampling_time;
    std::string m_receive_type;
};

MTS_EXTERN_CLASS_RENDER(Receiver)
NAMESPACE_END(mitsuba)
