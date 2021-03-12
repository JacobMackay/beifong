#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/transmitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

#include <unistd.h>

NAMESPACE_BEGIN(mitsuba)
/**!

.. _transmitter-area:

Area light (:monosp:`area`)
---------------------------

.. pluginparameters::

 * - radiance
   - |spectrum|
   - Specifies the transmitted radiance in units of power per unit area per unit steradian.

This plugin implements an area light, i.e. a light source that transmits
diffuse illumination from the exterior of an arbitrary shape.
Since the transmission profile of an area light is completely diffuse, it
has the same apparent brightness regardless of the observer's viewing
direction. Furthermore, since it occupies a nonzero amount of space, an
area light generally causes scene objects to cast soft shadows.

To create an area light source, simply instantiate the desired
transmitter shape and specify an :monosp:`area` instance as its child:

.. code-block:: xml
    :name: sphere-light

    <shape type="sphere">
        <transmitter type="area">
            <spectrum name="radiance" value="1.0"/>
        </transmitter>
    </shape>

 */

template <typename Float, typename Spectrum>
class WignerTransmitter final : public Transmitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Transmitter, m_flags, m_shape, m_medium)//,
        // m_sig_amplitude, m_sig_repfreq, m_sig_t_ext, m_sig_f0, m_sig_is_delta)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    WignerTransmitter(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        m_antenna_texture = props.texture<Texture>("antenna_texture", Texture::D65(1.f));

        m_flags = +TransmitterFlags::Surface;
        if (m_antenna_texture->is_spatially_varying())
            m_flags |= +TransmitterFlags::SpatiallyVarying;

        m_signal = props.string("signaltype", "cw");
        m_resample_freq = props.bool_("resample_freq", false);

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
        } else if (m_signal == "cw"){
            m_sig_amplitude = props.float_("amplitude", 1.f);
            m_sig_f_centre = props.float_("freq_centre", 1.f);
            m_sig_f_ext = props.float_("freq_ext", 0.f);
            m_sig_phi0 = props.float_("phase", 0.f);
            m_sig_is_delta = props.bool_("sig_is_delta", true);
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

    // Return the spectral flux/instantaneous signal power spectral density in
    // units V^2/Hz
    // ------------------------------------------------------------------------
    Spectrum eval_signal(Float time, Wavelength frequency) const {

        Spectrum result(0.f);
        Float t;
        // Float t_hat;
        Float ti;
        // Wavelength f_hat;
        Wavelength f_i;

        if (m_signal == "linfmcw") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            f_i = (m_sig_f_centre) + (m_sig_f_ext/(m_sig_t_ext))*(t - ti);
            result = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - f_i[0], m_sig_t_ext, m_sig_amplitude),
                0.f);
        } else if (m_signal == "pulse") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            f_i = m_sig_f_centre;
            result = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - f_i[0], m_sig_t_ext, m_sig_amplitude),
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
            // Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
            // Float f_i = (m_sig_f_centre - m_sig_t_ext/2) + (m_sig_f_ext/m_sig_t_ext*t_norm);
            Float t_i = 0 + m_sig_t_ext/2;
            Float t_norm = math::fmodulo(time - t_i, rcp(m_sig_repfreq));
            Wavelength f_i = (m_sig_f_centre - m_sig_f_ext/2) + (m_sig_f_ext/m_sig_t_ext*t_norm);
            frequencies = f_i;
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
        if (m_sig_is_delta == true) {
            auto [frequencies, freq_weight] = sample_delta_frequency(time);
            return {frequencies, freq_weight};
        } else {
            auto freq_sample = math::sample_shifted<Wavelength>(sample);
            Wavelength frequencies =
                freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
            Spectrum freq_weight = eval_signal(time, frequencies);
            return {frequencies, freq_weight};
        }
    }
    // ========================================================================

    // Return the spectral radiance of an impacting ray in [W/(sr*m^2*sm)]
    // ------------------------------------------------------------------------
    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

        // 1. Evaluate the signal power in [V^2/Hz] -----------------
        Spectrum signal_power;
        if (m_resample_freq == true) {
            // Resample frequency ------------------
            // This is used for forcing results, usually for cw transmission
            auto [frequencies, freq_weight] = sample_delta_frequency(si.time);
            signal_power = freq_weight;
            const_cast<SurfaceInteraction3f&>(si).wavelengths = MTS_C*rcp(frequencies)*1e9;
            std::cout << "resample eval" << signal_power << std::endl;
            // =====================================
        } else {
            signal_power = eval_signal(si.time, MTS_C*rcp(si.wavelengths*1e-9));
        }
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 3a. Positional part from surface in [1/m^2] -
        Spectrum geom_gain = m_antenna_texture->eval(si, active) *
            rcp(m_shape->surface_area());
        // =============================
        // 3b. Directional part from WDF in [1/sr] --
        DirectionSample3f ds(si);
        ds.d *= -1.f;
        DirectionSample3f ws = m_shape->sample_wigner(ds, si.wavelengths, active);
        geom_gain *= ws.pdf;
        // =============================

        // Return the radiance in W/(sr*m^2*sm) --------------
        // Return the radiance in W/(sr*m^2*Hz) --------------
        return select(
            Frame3f::cos_theta(si.wi) > 0.f,
            signal_power * transmission_gain * geom_gain,
            0.f
        );
        // ======================================================
    }
    // ========================================================================

    // Sample a ray and return the power
    // emanating at a position, in a direction and with a wavelength
    // ------------------------------------------------------------------------
    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float frequency_sample,
                                          const Point2f &position_sample,
                                          const Point2f &direction_sample,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Evaluate the signal power in [V^2/Hz] ------------
        // std::pair<Wavelength, Spectrum> res =
        //     sample_frequency(time, frequency_sample);
        // Wavelength frequencies = res.first;
        // Spectrum signal_power= res.second;
        std::cout << "sample ray" << std::endl;
        auto [frequencies, signal_power] =
            sample_frequency(time, frequency_sample);
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Convert to wavelength ----------------------------
        Wavelength wavelength = MTS_C*rcp(frequencies)*1e9;
        // =====================================================

        // 4. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 4a. Positional part from surface in [1/m^2] --
        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;
        Float pdf = 1.f;
        // Two strategies to sample spatial component based on
        // 'm_antenna_texture'
        if (!m_antenna_texture->is_spatially_varying()) {
            PositionSample3f ps =
                m_shape->sample_position(time, position_sample, active);

            // Radiance not spatially varying, use area-based sampling of shape
            si = SurfaceInteraction3f(ps, zero<Wavelength>());
            pdf = ps.pdf;
        } else {
            // Ipmortance sample texture
            std::tie(si.uv, pdf) =
                m_antenna_texture->sample_position(position_sample, active);
            active &= neq(pdf, 0.f);

            si = m_shape->eval_parameterization(Point2f(si.uv), active);
            active &= si.is_valid();

            pdf /= norm(cross(si.dp_du, si.dp_dv));
        }
        // 4b. Directional part from WDF in [1/sr] --
        DirectionSample3f ds(si);
        // ds.d = si.to_world(warp::square_to_cosine_hemisphere(direction_sample));
        ds.d = warp::square_to_cosine_hemisphere(direction_sample);
        DirectionSample3f ws =
            m_shape->sample_wigner(ds, si.wavelengths, active);
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
        // return std::make_pair(
        //     Ray3f(si.p, ds.d, time, wavelength),
        //     unpolarized<Spectrum>(signal_power * transmission_gain * geom_gain * extents)
        // );
        return std::make_pair(
            Ray3f(si.p, si.to_world(ds.d), time, wavelength),
            unpolarized<Spectrum>(signal_power * transmission_gain * geom_gain * extents)
        );
        // ====================================================
    }
// ============================================================================

    // Return the spectral radiant intensity in [W/(sr*m^2*sm)] passing through
    // a point it.
    // This includes the radiance of the transmitter, the cosine of the
    // direction and r^2 falloff.
    // ------------------------------------------------------------------------
    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area transmitter without an associated Shape.");

        // 1. Evaluate the signal power in [V^2/Hz] -----------------
        Spectrum signal_power;
        if (m_resample_freq == true) {
            // Resample frequency ------------------
            // This is used for forcing results, usually for cw transmission
            auto [frequencies, freq_weight] = sample_delta_frequency(it.time);
            signal_power = freq_weight;
            const_cast<Interaction3f&>(it).wavelengths = MTS_C*rcp(frequencies)*1e9;
            std::cout << "resample sample direction" << signal_power << std::endl;
            // =====================================
        } else {
            signal_power = eval_signal(it.time, MTS_C*rcp(it.wavelengths*1e-9));
        }
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        DirectionSample3f ds;
        Spectrum geom_gain;
        // 3a. Positional part from surface in [m^2/r^2] -
        // One of two very different strategies is used depending on
        // 'm_antenna_texture'
        if (!m_antenna_texture->is_spatially_varying()) {
            // Texture is uniform, try to importance sample the shape wrt.
            // solid angle at 'it'. Returns pdf = r^2/(Area*cosθ.)
            ds = m_shape->sample_direction(it, sample, active);
            active &= dot(ds.d, ds.n) < 0.f && neq(ds.pdf, 0.f);

            SurfaceInteraction3f si(ds, it.wavelengths);
            geom_gain = m_antenna_texture->eval(si, active) / ds.pdf;
        } else {
            // Importance sample the texture, then map onto the shape
            auto [uv, pdf] = m_antenna_texture->sample_position(sample, active);
            active &= neq(pdf, 0.f);

            SurfaceInteraction3f si = m_shape->eval_parameterization(uv, active);
            si.wavelengths = it.wavelengths;
            active &= si.is_valid();

            ds.p = si.p;
            ds.n = si.n;
            ds.uv = si.uv;
            ds.time = it.time;
            ds.delta = false;
            ds.d = ds.p - it.p;

            Float dist_squared = squared_norm(ds.d);
            ds.dist = sqrt(dist_squared);
            ds.d /= ds.dist;

            Float dp = dot(ds.d, ds.n);
            active &= dp < 0;
            ds.pdf = select(active, pdf / norm(cross(si.dp_du, si.dp_dv)) *
                                        dist_squared / -dp, 0.f);

            geom_gain = m_antenna_texture->eval(si, active) / ds.pdf;
        }
        // =============================
        // 3b. Directional part from WDF in [1/sr] --
        ds.d *= -1.f;
        DirectionSample3f ws = m_shape->sample_wigner(ds, it.wavelengths, active);
        ds.d *= -1.f;
        geom_gain *= ws.pdf;
        ds.pdf *= ws.pdf;
        // =============================

        // 4. Evaluate various extents to find ds radiance -----
        Float extents = math::Pi<Float>;
        // ====================================================

        // std::cout << (signal_power*transmission_gain*geom_gain* extents) << std::endl;

        // 6. Return the ds, and ds radiant intensity ---------
        ds.object = this;
        return { ds,
            unpolarized<Spectrum>((signal_power*transmission_gain*geom_gain)
         * extents) & active };
        // ====================================================
    }
    // ========================================================================

    // Return the probability per steradian
    // Ignore signal for now?
    // ------------------------------------------------------------------------
    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);
        Float dp = dot(ds.d, ds.n);
        active &= dp < 0.f;

        Float value;
        // 1. Evaluate the geometric probability in [1/(sr)] ------
        // 1a. Positional part from surface in [m^2/r^2] -
        if (!m_antenna_texture->is_spatially_varying()) {
            value = m_shape->pdf_direction(it, ds, active);
        } else {
            // This surface intersection would be nice to avoid..
            SurfaceInteraction3f si = m_shape->eval_parameterization(ds.uv, active);
            active &= si.is_valid();

            value = m_antenna_texture->pdf_position(ds.uv, active) * sqr(ds.dist) /
                    (norm(cross(si.dp_du, si.dp_dv)) * -dp);
        }
        // =============================
        // 1b. Directional part from WDF in [1/sr] --
        DirectionSample3f ds_ = ds;
        ds_.d *= -1.f;
        DirectionSample3f ws = m_shape->sample_wigner(ds_, it.wavelengths, active);
        value *= ws.pdf;
        // =============================

        // 2. Evaluate various extents to find ds probability -----
        Float extents = math::Pi<Float>;
        // ====================================================

        return select(active, value*extents, 0.f);
    }
    // ========================================================================

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("antenna_texture", m_antenna_texture.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "WignerTransmitter[" << std::endl
            << "  antenna_texture = " << string::indent(m_antenna_texture) << "," << std::endl
            << "  surface_area = ";
        if (m_shape) oss << m_shape->surface_area();
        else         oss << "  <no shape attached!>";
        oss << "," << std::endl;
        if (m_medium) oss << string::indent(m_medium);
        else         oss << "  <no medium attached!>";
        oss << std::endl << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Float m_gain;
    ref<Texture> m_antenna_texture;
    bool m_resample_freq;
    std::string m_signal;
    Float m_sig_amplitude;
    Float m_sig_repfreq;
    Float m_sig_t_ext;
    Float m_sig_f_centre;
    Float m_sig_f_ext;
    Float m_sig_phi0;
    bool m_sig_is_delta;
};

MTS_IMPLEMENT_CLASS_VARIANT(WignerTransmitter, Transmitter)
MTS_EXPORT_PLUGIN(WignerTransmitter, "Wigner transmitter")
NAMESPACE_END(mitsuba)
