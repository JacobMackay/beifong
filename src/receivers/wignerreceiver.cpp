#include <mitsuba/core/fwd.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/receiver.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _sensor-irradiancemeter:

Irradiance meter (:monosp:`irradiancemeter`)
--------------------------------------------

.. pluginparameters::

 * - none

This sensor plugin implements an irradiance meter, which measures
the incident power per unit area over a shape which it is attached to.
This sensor is used with films of 1 by 1 pixels.

If the irradiance meter is attached to a mesh-type shape, it will measure the
irradiance over all triangles in the mesh.

This sensor is not instantiated on its own but must be defined as a child
object to a shape in a scene. To create an irradiance meter,
simply instantiate the desired sensor shape and specify an
:monosp:`irradiancemeter` instance as its child:

.. code-block:: xml
    :name: sphere-meter

    <shape type="sphere">
        <sensor type="irradiancemeter">
            <!-- film -->
        </sensor>
    </shape>
*/

MTS_VARIANT class Wignerreceiver final : public Receiver<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Receiver, m_adc, m_world_transform, m_shape, m_receive_type)
    MTS_IMPORT_TYPES(Shape)

    Wignerreceiver(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The wigner receiver inherits this transformation "
                  "from its parent shape.");

        // if (m_adc->size() != ScalarPoint2i(1, 1))
        //     Throw("This sensor only supports films of size 1x1 Pixels!");

        if (m_adc->reconstruction_filter()->radius() >
            0.5f + math::RayEpsilon<Float>)
            Log(Warn, "This sensor should only be used with a reconstruction filter"
               "of radius 0.5 or lower(e.g. default box)");

               // m_receive_type = props.string("signaltype", "raw");

               if (m_receive_type == "raw") {
                   m_sig_f_centre = props.float_("freq_centre", 1.f);
                   m_sig_f_ext = props.float_("freq_ext", 1.f);
                   m_gain = props.float_("gain", 1.f);
               } else {
                   m_signal = props.string("signaltype", "cw");
                   if (m_signal == "linfmcw"){
                       m_sig_amplitude = props.float_("amplitude", 1.f);
                       m_sig_repfreq = props.float_("crf", 1.f);
                       m_sig_t_ext = props.float_("chirp_len", 1.f);
                       m_sig_f_centre = props.float_("freq_centre", 1.f);
                       m_sig_f_ext = props.float_("freq_sweep", 1.f);
                       m_sig_phi0 = props.float_("phase", 0.f);
                       m_sig_is_delta = props.bool_("sig_is_delta", true);
                       m_gain = props.float_("gain", 1.f);
                   } else if (m_signal == "pulse") {
                       m_sig_amplitude = props.float_("amplitude", 1.f);
                       m_sig_repfreq = props.float_("prf", 1.f);
                       m_sig_t_ext = props.float_("pulse_len", 1.f);
                       m_sig_f_centre = props.float_("freq_centre", 1.f);
                       m_sig_f_ext = props.float_("freq_ext", 1.f);
                       m_sig_phi0 = props.float_("phase", 0.f);
                       m_sig_is_delta = props.bool_("sig_is_delta", false);
                       m_gain = props.float_("gain", 1.f);
                   } else {
                       m_sig_amplitude = props.float_("amplitude", 1.f);
                       m_sig_repfreq = props.float_("prf", 1.f);
                       m_sig_t_ext = props.float_("pulse_len", 1.f);
                       m_sig_f_centre = props.float_("freq_centre", 1.f);
                       m_sig_f_ext = static_cast<Float>(rcp(m_sig_t_ext));
                       m_sig_phi0 = props.float_("phase", 0.f);
                       m_sig_is_delta = props.bool_("sig_is_delta", true);
                       m_gain = props.float_("gain", 1.f);
                   }
               }
    }

    // Return the spectral flux/instantaneous signal power spectral density in
    // units V^2/Hz
    // ------------------------------------------------------------------------
    Spectrum eval_signal(Float time, Wavelength frequency) const {

        Spectrum result(0.f);
        Float t_norm;
        Float t_hat;
        Wavelength f_hat;

        if (m_signal == "linfmcw") {
            t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
            t_hat = t_norm/m_sig_t_ext;
            f_hat = frequency - (m_sig_f_centre - m_sig_f_ext/2  + 0.5*m_sig_f_ext/m_sig_t_ext*t_norm);

            result = select(math::rect(t_hat) > 0.f,
                2*m_sig_amplitude*m_sig_amplitude * m_sig_t_ext*math::tri(t_hat) *
                math::sinc(math::TwoPi<Float>*f_hat[0]*m_sig_t_ext*math::tri(t_hat)),
                0.f);

        } else if (m_signal == "pulse") {
            t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
            t_hat = t_norm/m_sig_t_ext;
            f_hat = frequency - m_sig_f_centre;

            result = select(math::rect(t_hat) > 0.f,
                2*m_sig_amplitude*m_sig_amplitude * m_sig_t_ext*math::tri(t_hat) *
                math::sinc(math::TwoPi<Float>*f_hat[0]*m_sig_t_ext*math::tri(t_hat)),
                0.f);
        } else if (m_signal == "cw"){
            result = m_sig_amplitude*m_sig_amplitude;
        } else {
            result = m_sig_amplitude*m_sig_amplitude;
        }

        return result;
    }
    // ========================================================================

    // Return a frequency sample drawn from a delta distribution in units V^2/Hz
    // ------------------------------------------------------------------------
    std::pair<Wavelength, Spectrum> sample_delta_frequency(Float time) const {
        Wavelength frequencies;
        if (m_signal == "linfmcw") {
            // Sample a frequency from 1st order chirp -------
            Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
            frequencies = (m_sig_f_centre - m_sig_f_ext/2)
                + 0.5*m_sig_f_ext/m_sig_t_ext*t_norm;
        } else if (m_signal == "cw") {
            frequencies = m_sig_f_centre;
        }

        return {frequencies, eval_signal(time, frequencies)};
        // ===============================================
    }
    // ========================================================================

    // Return a frequency sample from signal at time t
    // in addition to the spectral/signal weight/power in V^2/Hz
    // ------------------------------------------------------------------------
    std::pair<Wavelength, Spectrum>
    sample_frequency(Float time, Float sample) const {
        // Make a pair, frequencies and frequency weight
        std::pair<Wavelength, Spectrum> result;

        if (m_receive_type == "raw"){
            auto freq_sample = math::sample_shifted<Wavelength>(sample);
            result.first =
                freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
            result.second = 1.f;
        } else {
            if (m_sig_is_delta == true) {
                result = sample_delta_frequency(time);
            } else {
                auto freq_sample = math::sample_shifted<Wavelength>(sample);
                result.first =
                    freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
                result.second = eval_signal(time, result.first);
            }
        }
        return result;
    }
    // ========================================================================


    // Sample a ray and return the power
    // emanating at a position, in a direction and with a wavelength
    // ------------------------------------------------------------------------
    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float frequency_sample,
                            const Point2f & position_sample,
                            const Point2f & direction_sample,
                            Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Evaluate the signal power in [V^2/Hz] ------------
        auto [frequencies, signal_power] =
            sample_frequency(time, frequency_sample);
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum reception_gain = m_gain;
        // =====================================================

        // 3. Convert to wavelength ----------------------------
        Wavelength wavelength = MTS_C*rcp(frequencies)*1e9;
        // =====================================================

        // 4. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 4a. Positional part from surface in [1/m^2] --
        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;
        Float pdf = 1.f;
        // 'm_antenna_texture' not a thing in RX for now
        PositionSample3f ps =
            m_shape->sample_position(time, position_sample, active);

        // Radiance not spatially varying, use area-based sampling of shape
        si = SurfaceInteraction3f(ps, zero<Wavelength>());
        pdf = ps.pdf;
        // 4b. Directional part from WDF in [1/sr] --
        DirectionSample3f ds(si);
        // ds.d = si.to_world(warp::square_to_cosine_hemisphere(direction_sample));
        ds.d = warp::square_to_cosine_hemisphere(direction_sample);
        DirectionSample3f ws =
            m_shape->sample_wigner(ds, wavelength, active);
        Spectrum geom_gain = ws.pdf * pdf;
        // ===========================================
        // ====================================================

        // 5. Evaluate various extents
        // to find individual ray power -----------------------
        Spectrum extents = m_shape->surface_area() * math::Pi<Float>;
        if (!m_sig_is_delta) {
            extents *= MTS_C*rcp(m_sig_f_ext)*1e9;
        }
        // ====================================================

        // 6. Return the ray, and ray power -------------------
        return std::make_pair(
            Ray3f(si.p, si.to_world(ds.d), time, wavelength),
            unpolarized<Spectrum>(signal_power * reception_gain * geom_gain * extents)
        );
        // ====================================================
    }
    // ============================================================================

    // Lazy, these aren't touched atm.
    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        std::cout << "Touched RX sample_direction" << std::endl;
        return std::make_pair(m_shape->sample_direction(it, sample, active), math::Pi<ScalarFloat>);
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        std::cout << "Touched RX pdf_direction" << std::endl;
        return m_shape->pdf_direction(it, ds, active);
    }

    Spectrum eval(const SurfaceInteraction3f &/*si*/, Mask /*active*/) const override {
        std::cout << "Touched RX eval" << std::endl;
        return math::Pi<ScalarFloat> / m_shape->surface_area();
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Wignerreceiver[" << std::endl
            << "  shape = " << m_shape << "," << std::endl
            << "  ADC = " << m_adc << "," << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    // std::string m_receive_type;
    std::string m_signal;
    Float m_gain;
    Float m_sig_amplitude;
    Float m_sig_repfreq;
    Float m_sig_t_ext;
    Float m_sig_f_centre;
    Float m_sig_f_ext;
    Float m_sig_phi0;
    bool m_sig_is_delta;
};

MTS_IMPLEMENT_CLASS_VARIANT(Wignerreceiver, Receiver)
MTS_EXPORT_PLUGIN(Wignerreceiver, "Wignerreceiver");
NAMESPACE_END(mitsuba)
