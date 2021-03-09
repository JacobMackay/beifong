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
        m_change_freq = props.bool_("change_freq", false);

        if (m_signal == "linfmcw"){
            m_sig_amplitude = props.float_("amplitude", 1.f);
            m_sig_repfreq = props.float_("crf", 1.f);
            m_sig_t_ext = props.float_("chirp_len", 1.f);
            m_sig_f_centre = props.float_("freq_centre", MTS_C/((MTS_WAVELENGTH_MAX + MTS_WAVELENGTH_MIN)/2));
            m_sig_f_ext = props.float_("freq_sweep", MTS_C/(MTS_WAVELENGTH_MAX - MTS_WAVELENGTH_MIN));
            m_sig_phi0 = props.float_("phase", 0.f);
            m_sig_is_delta = props.bool_("sig_is_delta", true);
        } else if (m_signal == "pulse") {
            m_sig_amplitude = props.float_("amplitude", 1.f);
            m_sig_repfreq = props.float_("prf", 1.f);
            m_sig_t_ext = props.float_("pulse_len", 1.f);
            m_sig_f_centre = props.float_("freq_centre", MTS_C/((MTS_WAVELENGTH_MAX + MTS_WAVELENGTH_MIN)/2));
            m_sig_f_ext = rcp(m_sig_t_ext);
            m_sig_phi0 = props.float_("phase", 0.f);
            m_sig_is_delta = props.bool_("sig_is_delta", false);
        } else {
            m_sig_amplitude = props.float_("amplitude", 1.f);
            m_sig_repfreq = props.float_("prf", 1.f);
            m_sig_t_ext = props.float_("pulse_len", 1.f);
            m_sig_f_centre = props.float_("freq_centre", MTS_C/((MTS_WAVELENGTH_MAX + MTS_WAVELENGTH_MIN)/2));
            m_sig_f_ext = rcp(m_sig_t_ext);
            m_sig_phi0 = props.float_("phase", 0.f);
            m_sig_is_delta = props.bool_("sig_is_delta", true);
        }


    }

    // Return the spectral flux/instantaneous signal power spectral density in
    // units V^2/Hz
    // ------------------------------------------------------------------------
    Float eval_signal(Float time, Float frequency) const {

        Float result(0.f);
        Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
        Float t_hat;
        Float f_hat;

        if (m_signal == "linfmcw") {
            t_hat = t_norm/m_sig_t_ext;
            f_hat = frequency - (m_sig_f_centre - m_sig_f_ext/2  + 0.5*m_sig_f_ext/m_sig_t_ext*t_norm);

            result = select(math::rect(t_hat) > 0.f,
                2*m_sig_amplitude*m_sig_amplitude * m_sig_t_ext*math::tri(t_hat) *
                math::sinc(math::TwoPi<Float>*f_hat*m_sig_t_ext*math::tri(t_hat)),
                0.f);
        } else if (m_signal == "pulse") {
            t_hat = t_norm/m_sig_t_ext;
            f_hat = frequency - m_sig_f_centre;

            result = select(math::rect(t_hat) > 0.f,
                2*m_sig_amplitude*m_sig_amplitude * m_sig_t_ext*math::tri(t_hat) *
                math::sinc(math::TwoPi<Float>*f_hat*m_sig_t_ext*math::tri(t_hat)),
                0.f);
        } else {
            result = m_sig_amplitude*m_sig_amplitude;
        }

        return result;
    }
    // ========================================================================

    // Return a frequency sample drawn from a delta distribution in units V^2/Hz
    // ------------------------------------------------------------------------
    std::pair<wavelength_t<Spectrum>, Spectrum> sample_delta_frequency(Float time) const {
        // Sample a frequency from 1st order chirp -------
        Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
        Float frequencies = (m_sig_f_centre - m_sig_f_ext/2)
            + 0.5*m_sig_f_ext/m_sig_t_ext*t_norm;
        return {frequencies, eval_signal(time, frequencies)};
        // ===============================================
    }
    // ========================================================================

    // Return a frequency sample from signal at time t
    // in addition to the spectral/signal weight/power in V^2/Hz
    // ------------------------------------------------------------------------
    std::pair<wavelength_t<Spectrum>, Spectrum>
    sample_frequency(Float time, Float sample) const {
        if (m_sig_is_delta == true) {
            auto [frequencies, freq_weight] = sample_delta_frequency(time);
            return {frequencies, freq_weight};
        } else {
            auto freq_sample = math::sample_shifted<wavelength_t<Spectrum>>(sample);
            wavelength_t<Spectrum> frequencies =
                freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
            Spectrum freq_weight = eval_signal(time, frequencies);
            return {frequencies, freq_weight};
        }
    }
    // ========================================================================

    // Return the spectral radiance of an impacting ray in W/(sr * m^2 * sm)
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
            const_cast<SurfaceInteraction3f&>(si).wavelengths = to_wavelength(frequencies);
            // =====================================
        } else {
            signal_power = eval_signal(time, to_frequency(si.wavelengths[0]));
        }
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 3a. Directional part from WDF in [1/sr] --
        DirectionSample3f ds(si);
        ds.d *= -1.f;
        DirectionSample3f ws = m_shape->sample_wigner(ds, si.wavelengths, active);
        Float geom_gain = ws.pdf;
        // =============================
        // 3b. Areal part from surface in [1/m^2] ---
        geom_gain *= m_antenna_texture->eval(si, active) *
            rcp(m_shape->surface_area());
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
        auto [frequencies, signal_power] =
            sample_frequency(time, frequency_sample);
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Convert to wavelength ----------------------------
        Wavelength wavelength = to_wavelength(frequencies);
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
        ds.d = si.to_world(warp::square_to_cosine_hemisphere(direction_sample));
        DirectionSample3f ws =
            m_shape->sample_wigner(ds, si.wavelengths, active);
        Float geom_gain = ws.pdf * pdf;
        // ===========================================
        // ====================================================

        // 5. Evaluate various extents
        // to find individual ray power -----------------------
        Float extents = m_shape->surface_area() * math::Pi<Float>;
        if (!m_sig_is_delta) {
            extents *= to_wavelength(m_sig_f_ext);
        }
        // ====================================================

        return (signal_power * transmission_gain * geom_gain) * extents
        // return radiance * area * pi


        return std::make_pair(
            Ray3f(ws.p, ws.d, ws.time, wavelength),
            select(abs(ws.pdf)>math::Epsilon<Float>, unpolarized<Spectrum>(spec_weight) / ws.pdf, 0.f)
        );
        return std::make_pair(
            Ray3f(si.p, ds.d, time, wavelength),
            unpolarized<Spectrum>(spec_weight) * (math::Pi<Float> / pdf)
        );

        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;

        Float pdf = 1.f;

        // 1. Two strategies to sample spatial component based on 'm_antenna_texture'
        if (!m_antenna_texture->is_spatially_varying()) {
            PositionSample3f ps = m_shape->sample_position(time, sample2, active);

            // Radiance not spatially varying, use area-based sampling of shape
            si = SurfaceInteraction3f(ps, zero<Wavelength>());
            pdf = ps.pdf;
        } else {
            // Ipmortance sample texture
            std::tie(si.uv, pdf) = m_radiance->sample_position(sample2, active);
            active &= neq(pdf, 0.f);

            si = m_shape->eval_parameterization(Point2f(si.uv), active);
            active &= si.is_valid();

            pdf /= norm(cross(si.dp_du, si.dp_dv));
        }

        // 2. Sample directional component
        Vector3f local = warp::square_to_cosine_hemisphere(sample3);

        Wavelength wavelength;
        Spectrum spec_weight;

        if constexpr (is_spectral_v<Spectrum>) {
            std::tie(wavelength, spec_weight) = m_radiance->sample_spectrum(
                si, math::sample_shifted<Wavelength>(wavelength_sample), active);
        } else {
            wavelength = zero<Wavelength>();
            spec_weight = m_radiance->eval(si, active);
        }

        DirectionSample3f ds;

        // Somehow in whats happening now, perhaps the integrator? We're never
        // entering this function

        ds.p = si.p;
        ds.n = si.n;
        ds.uv = si.uv;
        ds.time = time;
        ds.delta = false;
        ds.d = si.to_world(warp::square_to_cosine_hemisphere(sample3));
        // ds.d = si.to_world(warp::square_to_uniform_hemisphere(sample3));
        Float dist_squared = squared_norm(ds.d);
        ds.dist = sqrt(dist_squared);
        ds.d /= ds.dist;

        ds.pdf = pdf;
        // Assuming wigner doesn't take look angle into account
        Float dp = abs_dot(ds.d, ds.n);
        ds.pdf *= select(neq(dp, 0.f), dist_squared / dp, 0.f);

        ds.object = this;

        // ds.d *= -1.f;

        DirectionSample3f ws = m_shape->sample_wigner(ds, wavelength, active);

        // find spectral radiance at pos/dir/λ, then upgrade to power

        // ws.d *= -1.f;

        // ray weight = spec_weight / m_inv_surface_area

        // return std::make_pair(
        //     Ray3f(si.p, si.to_world(local), time, wavelength),
        //     unpolarized<Spectrum>(spec_weight) * (math::Pi<Float> / pdf)
        // );
        // return std::make_pair(
        //     Ray3f(ws.p, ws.d, ws.time, wavelength),
        //     unpolarized<Spectrum>(spec_weight) * (math::Pi<Float> / ws.pdf)
        // );
        // return std::make_pair(
        //     Ray3f(ws.p, ws.d, ws.time, wavelength),
        //     unpolarized<Spectrum>(spec_weight) / ws.pdf
        // );

        // spec_weight *= 1/(4*math::Pi<Float>);
        // spec_weight *= m_shape->surface_area();

        // In the original area light, the weight was (spectral band)*(pi)*(Area)
        // This implies possibly two things:
        // Functions we sampled from are in units W/(sr m^2 Hz)
        // Well actually W/(sr m^2 nu)
        // This means the function returns the power not rdiance.
        // Or that we multiply by extent because integral?

        // original has radiance * spectral width * pi /(1/A)
        // ie returns power

        return std::make_pair(
            Ray3f(ws.p, ws.d, ws.time, wavelength),
            select(abs(ws.pdf)>math::Epsilon<Float>, unpolarized<Spectrum>(spec_weight) / ws.pdf, 0.f)
        );
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area transmitter without an associated Shape.");
        DirectionSample3f ds;
        Spectrum spec;

        DirectionSample3f ws;

        // One of two very different strategies is used depending on 'm_radiance'
        if (!m_radiance->is_spatially_varying()) {
            // Texture is uniform, try to importance sample the shape wrt. solid angle at 'it'
            ds = m_shape->sample_direction(it, sample, active);
            active &= dot(ds.d, ds.n) < 0.f && neq(ds.pdf, 0.f);

            SurfaceInteraction3f si(ds, it.wavelengths);
            // spec = m_radiance->eval(si, active) / ds.pdf;
            // Convert to wigner space. We want to take in the 'outgoing' ray.
            ds.d *= -1.f;
            ws = m_shape->sample_wigner(ds, it.wavelengths, active);
            ws.d *= -1.f;
            // There is a possibility to actually return the correct weight per wlen sample.
            spec = m_radiance->eval(si, active) / ws.pdf;
        } else {
            // Importance sample the texture, then map onto the shape
            auto [uv, pdf] = m_radiance->sample_position(sample, active);
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

            // spec = m_radiance->eval(si, active) / ds.pdf;

            ds.d *= -1.f;
            ws = m_shape->sample_wigner(ds, it.wavelengths, active);
            ws.d *= -1.f;

            // After/before taking the wigner sample, use the it.time to find a
            // phase component. Then it's a simple multiplication. If we were
            // going from transmitter to receiver, how would this be different?
            // It would be nice to have a symmetric calculation, but rays have
            // diff start times, yes they would leave with correct phase, but
            // that's not useful as some could be 0.
            // There is a possibility to actually return the correct weight per wlen sample.
            spec = m_radiance->eval(si, active) / ws.pdf;
        }

        ds.object = this;

        spec /= (4*math::Pi<Float>)*(math::Pi<Float>);

        spec *= eval_signal(it.time, MTS_C/(it.wavelengths[0]*1e-9));

        // std::cout << ds.pdf << spec << std::endl;


        // The original gives radiance / dspdf
        // radiance / (1/A * r^2/cos(θ))
        // radiance * A * cos(θ) / r^2
        // So, is this power in W, or is this radiance normalised for sampling?


        // spec *= 1/(4*math::Pi<Float>);
        // spec *= m_shape->surface_area();

        // return { ds, unpolarized<Spectrum>(spec) & active };
        // return { ws, unpolarized<Spectrum>(spec) & active };
        return { ws, select(abs(ws.pdf)>math::Epsilon<Float>, unpolarized<Spectrum>(spec) & active, 0.f) };
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);
        Float dp = dot(ds.d, ds.n);
        active &= dp < 0.f;

        Float value;
        if (!m_radiance->is_spatially_varying()) {
            value = m_shape->pdf_direction(it, ds, active);
        } else {
            // This surface intersection would be nice to avoid..
            SurfaceInteraction3f si = m_shape->eval_parameterization(ds.uv, active);
            active &= si.is_valid();

            value = m_radiance->pdf_position(ds.uv, active) * sqr(ds.dist) /
                    (norm(cross(si.dp_du, si.dp_dv)) * -dp);
        }

        DirectionSample3f ds2 = ds;
        ds2.pdf = value;
        ds2.d *= -1.f;
        // This looks after the sampling
        DirectionSample3f ws = m_shape->sample_wigner(ds2, it.wavelengths, active);

        // value *= ws.pdf;
        // value = 1;
        // value = abs(ws.pdf);
        value = ws.pdf;
        // value = ws.pdf*ws.pdf;
        // value = ws.pdf*(math::TwoPi<Float>*math::TwoPi<Float>)*m_shape->surface_area()*it.wavelengths[0]*1e-9*it.wavelengths[0]*1e-9;
        value *= value;

        // value = 1.f;

        // value *= 1/(4*math::Pi<Float>);

        active &= abs(ws.pdf) > math::Epsilon<Float>;

        // original pdf is r^2 / (A*cos(θ))
        // value is A*cos(θ)/r^2

        return select(active, value, 0.f);
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("radiance", m_radiance.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "WignerTransmitter[" << std::endl
            << "  radiance = " << string::indent(m_radiance) << "," << std::endl
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
    Spectrum m_gain;
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
