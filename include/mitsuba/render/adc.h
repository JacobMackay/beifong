#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/** \brief Abstract ADC base class - used to store samples
 * generated by \ref Integrator implementations.
 *
 * To avoid lock-related bottlenecks when rendering with many cores,
 * rendering threads first store results in a "signal block", which
 * is then committed to the ADC using the \ref put() method.
 */
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER ADC : public Object {
public:
    MTS_IMPORT_TYPES(SignalBlock, ReconstructionFilter)

    /// Configure the ADC for rendering a specified set of channels
    virtual void prepare(const std::vector<std::string> &channels) = 0;

    /// Merge a signal block into the ADC. This methods should be thread-safe.
    virtual void put(const SignalBlock *block) = 0;

    /// Develop the ADC and write the result to the previously specified
    // filename
    virtual void develop() = 0;

    /**
     * \brief Develop the contents of a subregion of the ADC and store
     * it inside the given bitmap
     *
     * This may fail when the ADC does not have an explicit representation
     * of the bitmap in question (e.g. when it is writing to a tiled EXR image)
     *
     * \return \c true upon success
     */
    virtual bool develop(
        const ScalarPoint2i  &offset,
        const ScalarVector2i &size,
        const ScalarPoint2i  &target_offset,
        Bitmap *target) const = 0;

    /// Return a bitmap object storing the developed contents of the ADC
    virtual ref<Bitmap> bitmap(bool raw = false) = 0;

    /// Set the target filename (with or without extension)
    virtual void set_destination_file(const fs::path &filename) = 0;

    /// Does the destination file already exist?
    virtual bool destination_exists(const fs::path &basename) const = 0;

    /**
     * Should regions slightly outside the signal range be sampled to improve
     * the quality of the reconstruction at the edges? This only makes
     * sense when reconstruction filters other than the box filter are used.
     */
    bool has_high_quality_edges() const { return m_high_quality_edges; }

    // =============================================================
    //! @{ \name Accessor functions
    // =============================================================

    /// Ignoring the crop window, return the resolution (number of range and
    // doppler bins) of the ADC
    const ScalarVector2i &size() const { return m_size; }

    /// Ignoring the crop window, return the resolution (number of range and
    // doppler bins) of the ADC
    // const ScalarFloat &sample_rate() const { return m_sample_rate; }

    /// Return the size of the crop window
    const ScalarVector2i &window_size() const { return m_window_size; }

    /// Return the offset of the crop window
    const ScalarPoint2i &window_offset() const { return m_window_offset; }

    /// Return the bandwidth of the ADC
    const ScalarVector2f &bandwidth() const { return m_bandwidth; }
    //
    // /// Return the centres of the ADC
    // const ScalarVector2f &centres() const { return m_centres; }

    /// Set the size and offset of the crop window.
    void set_window(const ScalarPoint2i &window_offset,
                         const ScalarVector2i &window_size);

    /// Return the image reconstruction filter (const version)
    const ReconstructionFilter *reconstruction_filter() const {
        return m_filter.get();
    }

    //! @}
    // =============================================================

    virtual std::string to_string() const override;

    MTS_DECLARE_CLASS()
protected:
    /// Create an ADC
    ADC(const Properties &props);

    /// Virtual destructor
    virtual ~ADC();

protected:
    ScalarVector2i m_size;
    // ScalarFloat m_sample_rate;
    ScalarVector2i m_window_size;
    ScalarPoint2i m_window_offset;
    ScalarVector2f m_bandwidth;
    // ScalarVector2f m_centres;
    bool m_high_quality_edges;
    ref<ReconstructionFilter> m_filter;
};

MTS_EXTERN_CLASS_RENDER(ADC)
NAMESPACE_END(mitsuba)