#include <mitsuba/render/adc.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT ADC<Float, Spectrum>::ADC(const Properties &props) : Object() {
    bool is_m_adc = string::to_lower(props.plugin_name()) == "madc";

    // Time/range and frequency/doppler ADC resolution in bins
    m_size = ScalarVector2i(
        props.int_("tr_bins", is_m_adc ? 1 : 2e8),
        props.int_("fd_bins", is_m_adc ? 1 : 2e8)
    );

    // Window specified in bins - by default, this matches the full ADC range.
    ScalarVector2i window_size = ScalarVector2i(
        props.int_("window_tr_bins", m_size.x()),
        props.int_("window_fd_bins", m_size.y())
    );

    ScalarPoint2i window_offset = ScalarPoint2i(
        props.int_("window_offset_tr", 0),
        props.int_("window_offset_fd", 0)
    );

    set_window(window_offset, window_size);

    // Time/range and frequency/doppler ADC bandwidth in s/m and Hz/(s/m)
    m_bandwidth = ScalarVector2f(
        (props.float_("tr_max", 1.f) - props.float_("tr_min", 1.f)),
        props.float_("fd_bandwidth", 1.f)
    );

    // Time/range and frequency/doppler ADC centres
    m_centres = ScalarVector2f(
        (props.float_("tr_max", 0.f) + props.float_("tr_min", 0.f))/2,
        props.float_("fd_centre", 0.f)
    );

    /* If set to true, regions slightly outside of the ADC range will also be
       sampled, which improves the signal quality at the edges especially with
       large reconstruction filters. */
    m_high_quality_edges = props.bool_("high_quality_edges", false);

    // Use the provided reconstruction filter, if any.
    for (auto &[name, obj] : props.objects(false)) {
        auto *rfilter = dynamic_cast<ReconstructionFilter *>(obj.get());
        if (rfilter) {
            if (m_filter)
                Throw("An ADC can only have one reconstruction filter.");
            m_filter = rfilter;
            props.mark_queried(name);
        }
    }

    if (!m_filter) {
        // No reconstruction filter has been selected. Load a Gaussian filter
        // by default
        m_filter = PluginManager::instance()
            ->create_object<ReconstructionFilter>(Properties("gaussian"));
    }
}

MTS_VARIANT ADC<Float, Spectrum>::~ADC() {}

MTS_VARIANT void ADC<Float, Spectrum>::
    set_window(const ScalarPoint2i &window_offset,
               const ScalarVector2i &window_size) {
        if (any(window_offset < 0 || window_size <= 0 ||
                window_offset + window_size > m_size))
            Throw("Invalid window specification!\n"
                  "offset %s + window size %s vs full size %s",
                  window_offset, window_size, m_size);

        m_window_size   = window_size;
        m_window_offset = window_offset;
}

MTS_VARIANT std::string ADC<Float, Spectrum>::to_string() const {
    std::ostringstream oss;
    oss << "ADC[" << std::endl
        << "  size = "        << m_size        << "," << std::endl
        << "  window_size = "   << m_window_size   << "," << std::endl
        << "  window_offset = " << m_window_offset << "," << std::endl
        << "  high_quality_edges = " << m_high_quality_edges << "," << std::endl
        << "  m_filter = " << m_filter << std::endl
        << "]";
    return oss.str();
}


MTS_IMPLEMENT_CLASS_VARIANT(ADC, Object, "adc")
MTS_INSTANTIATE_CLASS(ADC)
NAMESPACE_END(mitsuba)
