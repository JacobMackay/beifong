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
class PhasedTransmitter final : public Transmitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Transmitter, m_flags, m_shape, m_medium)//,
        // m_sig_amplitude, m_sig_repfreq, m_sig_t_ext, m_sig_f0, m_sig_is_delta)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    PhasedTransmitter(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        m_antenna_texture = props.texture<Texture>("antenna_texture", Texture::D65(1.f));

        m_flags = +TransmitterFlags::Surface;
        if (m_antenna_texture->is_spatially_varying())
            m_flags |= +TransmitterFlags::SpatiallyVarying;

        // Signal Components --------------------------------------------------
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
        // ====================================================================

        // Phased Components --------------------------------------------------
        // Array Details
        m_n_elems = props.int_("n_elems", 1);
        ScalarVector3f steer_vec = props.vector3f("steering_vector", ScalarVector3f());
        steer_vec = sin(steer_vec);
        ScalarTransform4f array_to_world = props.transform("array_loc", ScalarTransform4f());
        m_wid = props.vector3f("elem_dims", ScalarVector3f());

        // These are local
        ScalarVector3f elem_spacing = props.vector3f("elem_spacing", ScalarVector3f());
        ScalarVector3f elem_axis = props.vector3f("elem_axis", ScalarVector3f());
        ScalarPoint3f array_centre(0.f, 0.f, 0.f);

        std::vector<ScalarPoint3f> elem_locs;

        for (int i = 0; i < m_n_elems; i++) {
            if (math::modulo(m_n_elems,2)==0) {
                elem_locs.push_back(array_centre - elem_spacing*elem_axis*(i-(m_n_elems/2.f)+0.5));
            } else {
                elem_locs.push_back(array_centre - elem_spacing*elem_axis*(i-(m_n_elems-1.f)/2.f));
            }
        }

        ScalarVector3f dp_du, dp_dv;
        ScalarNormal3f normal;
        ScalarVector3f r_v;

        std::complex<ScalarFloat> myI(0,1);
        for (int i = 0; i < m_n_elems; i++) {
            for (int j = 0; j < m_n_elems; j++) {
                r_v = (elem_locs[i] + elem_locs[j])/2;
                m_r_dash.push_back(elem_locs[i] - elem_locs[j]);

                m_velem_to_world.push_back(
                    array_to_world *
                    ScalarTransform4f::translate(r_v) *
                    ScalarTransform4f::scale(ScalarVector3f(m_wid.x()/2, m_wid.y()/2, m_wid.z()))
                );
                m_velem_to_object.push_back(m_velem_to_world.back().inverse());

                dp_du = m_velem_to_world.back() * ScalarVector3f(2.f, 0.f, 0.f);
                dp_dv = m_velem_to_world.back() * ScalarVector3f(0.f, 2.f, 0.f);
                normal = normalize(m_velem_to_world.back() * ScalarNormal3f(0.f, 0.f, 1.f));
                m_velem_frame.push_back(ScalarFrame3f(dp_du, dp_dv, normal));

                Transform4f trafto1;
                m_dir_to_local_velem.push_back( trafto1.from_frame(
                    Frame3f(normalize(m_velem_frame.back().s),
                            normalize(m_velem_frame.back().t),
                            normalize(m_velem_frame.back().n))) );

                m_psi_dash.push_back(
                    // exp( math::TwoPi<ScalarFloat>*myI*
                    exp( myI*
                    static_cast<ScalarFloat>(rcp((MTS_WAVELENGTH_MAX - MTS_WAVELENGTH_MIN)*1e-9/2))*
                    dot(array_centre - m_r_dash.back(), steer_vec) ));
            }
        }
        // ====================================================================


    }

    // Signal Methods ---------------------------------------------------------
    // ------------------------------------------------------------------------

    // Return the spectral flux/instantaneous signal power spectral density in
    // units V^2/Hz
    // ------------------------------------------------------------------------
    // Spectrum eval_signal(Float time, Wavelength frequency) const {
    std::pair<Spectrum, Float> eval_signal(Float time, Wavelength frequency) const {

        std::pair<Spectrum, Float> result;
        Float t;
        Float ti;
        Wavelength fi;

        if (m_signal == "linfmcw") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            fi = m_sig_f_centre + (m_sig_f_ext/m_sig_t_ext)*(t - ti);
            result.first = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - fi[0], m_sig_t_ext, m_sig_amplitude),
                0.f);
            result.second = m_sig_phi0 + math::TwoPi<Float>*(time - ti)*(m_sig_f_centre +
                0.5*(m_sig_f_ext/m_sig_t_ext)*(time - ti));
        } else if (m_signal == "pulse") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            fi = m_sig_f_centre;
            result.first = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - fi[0], m_sig_t_ext, m_sig_amplitude),
                0.f);
            result.second = m_sig_phi0 + math::TwoPi<Float>*(t - ti)*m_sig_f_centre;
        } else if (m_signal == "cw"){
            result.first = m_sig_amplitude*m_sig_amplitude;
            result.second = m_sig_phi0 + math::TwoPi<Float>*(t - ti)*m_sig_f_centre;
        } else {
            result.first = m_sig_amplitude*m_sig_amplitude;
            result.second = m_sig_phi0 + math::TwoPi<Float>*(t - ti)*m_sig_f_centre;
        }

        result.second = 0.f;

        return result;
    }
    // ========================================================================

    // Return a frequency sample drawn from a delta distribution in units V^2/Hz
    // ------------------------------------------------------------------------
    // std::pair<Wavelength, Spectrum> sample_delta_frequency(Float time) const {
    std::tuple<Wavelength, Spectrum, Float> sample_delta_frequency(Float time) const {
        // return tuple with phase
        Wavelength frequencies;
        if (m_signal == "linfmcw") {
            // Sample a frequency from 1st order chirp -------
            Float t = math::fmodulo(time, rcp(m_sig_repfreq));
            Float ti = 0 + m_sig_t_ext/2;
            Wavelength fi = m_sig_f_centre + (m_sig_f_ext/m_sig_t_ext)*(t - ti);
            frequencies = fi;
        } else if (m_signal == "cw") {
            frequencies = m_sig_f_centre;
        }
        auto [weight, phase] = eval_signal(time, frequencies);
        weight = 1.f;
        return {frequencies, weight, phase};
        // ===============================================
    }
    // ========================================================================

    // Return a frequency sample from signal at time t
    // in addition to the spectral/signal weight/power in V^2/Hz
    // ------------------------------------------------------------------------
    // std::pair<Wavelength, Spectrum>
    std::tuple<Wavelength, Spectrum, Float>
    sample_frequency(Float time, Float frequency_sample) const {
        // return tuple with phase
        Wavelength frequency;
        if (m_sig_is_delta == true) {
            return sample_delta_frequency(time);
        } else {
            auto freq_sample = math::sample_shifted<Wavelength>(frequency_sample);
            frequency =
                freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
            auto [weight, phase] = eval_signal(time, frequency);
            return {frequency, weight, phase};
        }
    }
    // ========================================================================

    // ========================================================================
    // ========================================================================

    // Phased Methods ---------------------------------------------------------
    // ------------------------------------------------------------------------
    Spectrum W_rect_2D(const Point3f &r_hat,
                        const Normal3f &nu_hat,
                        const ScalarVector3f &wid) const {

        return 4*wid.x()*wid.y() * math::tri(r_hat.x())*math::tri(r_hat.y()) *
                math::sinc(math::TwoPi<ScalarFloat>*nu_hat.x()*wid.x()*math::tri(r_hat.x())) *
                math::sinc(math::TwoPi<ScalarFloat>*nu_hat.y()*wid.y()*math::tri(r_hat.y()));
        // return 4* math::tri(r_hat.x())*math::tri(r_hat.y()) *
        //         math::sinc(math::TwoPi<ScalarFloat>*nu_hat.x()*wid.x()*math::tri(r_hat.x())) *
        //         math::sinc(math::TwoPi<ScalarFloat>*nu_hat.y()*wid.y()*math::tri(r_hat.y()));
    }

    Spectrum sample_wigner(const DirectionSample3f &ds, Wavelength wavelength) const {
        Point3f r_hat;
        Normal3f nu_hat;
        std::complex<Spectrum> W(0, 0);
        std::complex<ScalarFloat> myI(0,1);

        for (int i = 0; i < m_n_elems*m_n_elems; i++) {
            r_hat = m_velem_to_object[i] * ds.p/2;

            if (all(math::jabs(r_hat.x())<=0.5 & math::jabs(r_hat.y())<=0.5)) {
                nu_hat = m_dir_to_local_velem[i].transform_affine(ds.d)*rcp(wavelength[0]*1e-9);

                W += hsum(W_rect_2D(r_hat, nu_hat, m_wid)[0])
                 * exp(math::TwoPi<ScalarFloat>*myI*hsum(dot(nu_hat, m_r_dash[i])))*m_psi_dash[i];
            }
        }
        return real(W);
    }
    // ========================================================================
    // ========================================================================

    // Return the spectral radiance of an impacting ray in [W/(sr*m^2*sm)]
    // ------------------------------------------------------------------------
    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

        // It seems that the direct path has flat power, no 1/r^2 dropoff
        // Makes sense kinda, the 1/r^2 happens by sampling. Lets try with it in.

        // If I was going to go crazy, eval would return phase.
        // Maybe it should, that's more in spirit of returning a phasor
        // BUT I don't think mitsuba does complex numbers atm.
        // Also, endpoint and python prototypes will hate it.
        // So no. Just do  hacks

        // std::cout << "Here eval" << std::endl;

        // 1. Evaluate the signal power in [V^2/Hz] -----------------
        Wavelength frequency;
        Spectrum signal_power;
        Float phase;
        if (m_resample_freq == true) {
            // Resample frequency ------------------
            // This is used for forcing results, usually for cw transmission
            // auto [frequencies, freq_weight] = sample_delta_frequency(si.time);
            // auto [frequency, signal_power, phase] = sample_delta_frequency(si.time);
            auto [f, s, p] = sample_delta_frequency(si.time);
            frequency = f;
            signal_power = s;
            phase = p;
            const_cast<SurfaceInteraction3f&>(si).wavelengths = MTS_C*rcp(frequency)*1e9;
            const_cast<SurfaceInteraction3f&>(si).phase = phase;
            // std::cout << "resample eval" << signal_power << std::endl;
            // =====================================
        } else {
            // signal_power = eval_signal(si.time, MTS_C*rcp(si.wavelengths*1e-9));
            // auto [signal_power, phase] = eval_signal(si.time, MTS_C*rcp(si.wavelengths*1e-9));
            auto [s, p] = eval_signal(si.time, MTS_C*rcp(si.wavelengths*1e-9));
            signal_power = s;
            phase = p;
            const_cast<SurfaceInteraction3f&>(si).phase = phase;
        }
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/??] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 3a. Positional part from surface in [1/m^2] -
        // Original::
        Spectrum geom_gain = m_antenna_texture->eval(si, active) *
            rcp(m_shape->surface_area());
        // // Precaution::
        // Spectrum geom_gain = m_antenna_texture->eval(si, active);
        // =============================
        // 3b. Directional part from WDF in [1/sr] --
        DirectionSample3f ds(si);
        ds.d *= -1.f;
        // DirectionSample3f ws = m_shape->sample_wigner(ds, si.wavelengths, active);
        // geom_gain *= ws.pdf;
        geom_gain *= sample_wigner(ds, si.wavelengths);
        // =============================

        // Extent because mitsuba is inconsistent
        // Spectrum extent = rcp(m_shape->surface_area())/math::Pi<float>;
        // I feel like something catastrophic happens when things dip below 1.
        // Spectrum extent = rcp(m_shape->surface_area());
        // Currently a sample
        // Spectrum extent = rcp(m_shape->surface_area())*math::Pi<float>;
        // Spectrum extent = rcp(m_shape->surface_area())*math::Pi<float>*rcp(si.t*si.t);
        // ratio of intensities = d1^2/d2^2

        // Return the radiance in W/(sr*m^2*sm) --------------
        // Return the radiance in W/(sr*m^2*Hz) --------------
        // The last statement shouldn't be necessary, but its fucked atm
        return select(
            Frame3f::cos_theta(si.wi) > 0.f,
            // Dodgy
            // signal_power * transmission_gain * geom_gain * rcp(m_shape->surface_area()),
            // // Precaution:
            // signal_power * transmission_gain * geom_gain * math::TwoPi<Float>,
            // Original:
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
        // std::cout << "sample ray" << std::endl;
        Wavelength frequency;
        Spectrum signal_power;
        Float phase;
        // auto [frequencies, signal_power] =
        auto [f, s, p] =
            sample_frequency(time, frequency_sample);
        frequency = f;
        signal_power = s;
        phase = p;
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/??] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Convert to wavelength ----------------------------
        Wavelength wavelength = MTS_C*rcp(frequency)*1e9;
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
        // // DirectionSample3f ws =
        // //     m_shape->sample_wigner(ds, si.wavelengths, active);
        // DirectionSample3f ws =
        //     m_shape->sample_wigner(ds, wavelength, active);
        // Spectrum geom_gain = ws.pdf * pdf;
        Spectrum geom_gain = sample_wigner(ds, wavelength) * pdf;
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
        std::cout << "Sample here??" << std::endl;
        return std::make_pair(
            Ray3f(si.p, si.to_world(ds.d), time, phase, wavelength),
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

        // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // The wavelength is not exact, but skip for now
        DirectionSample3f ds;
        Spectrum geom_gain;
        // 3a. Positional part from surface in [m^2/r^2] -
        // One of two very different strategies is used depending on
        // 'm_antenna_texture'
        if (!m_antenna_texture->is_spatially_varying()) {
            // Texture is uniform, try to importance sample the shape wrt.
            // solid angle at 'it'. Returns pdf = r^2/(Area*cos??.)
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
        // Update the time when sampled by the environment
        // Don't if we're ON the transmitter...actually this isn't too bad
        if (all(ds.dist>5e-7)) {
            ds.time += -ds.dist / MTS_C;
            // const_cast<Interaction3f&>(it).time = ds.time;
        }

        // 1. Evaluate the signal power in [V^2/Hz] -----------------
        Wavelength frequency;
        Spectrum signal_power;
        Float phase;
        if (m_resample_freq == true) {
            // Resample frequency ------------------
            // This is used for forcing results, usually for cw transmission
            // auto [f, s, p] = sample_delta_frequency(it.time);
            auto [f, s, p] = sample_delta_frequency(ds.time);
            frequency = f;
            signal_power = s;
            phase = p;
            const_cast<Interaction3f&>(it).wavelengths = MTS_C*rcp(frequency)*1e9;
            const_cast<Interaction3f&>(it).phase = phase;
            // std::cout << "resample sample direction" << signal_power << std::endl;
            // =====================================
        } else {
            // auto [s, p] = eval_signal(it.time, MTS_C*rcp(it.wavelengths*1e-9));
            auto [s, p] = eval_signal(ds.time, MTS_C*rcp(it.wavelengths*1e-9));
            signal_power = s;
            phase = p;
            const_cast<Interaction3f&>(it).phase = phase;
            // it.phase = phase
        }
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/??] -------------------
        Spectrum transmission_gain = m_gain;
        // Spectrum transmission_gain = 1.f;
        // =====================================================

        // // 3. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // DirectionSample3f ds;
        // Spectrum geom_gain;
        // // 3a. Positional part from surface in [m^2/r^2] -
        // // One of two very different strategies is used depending on
        // // 'm_antenna_texture'
        // if (!m_antenna_texture->is_spatially_varying()) {
        //     // Texture is uniform, try to importance sample the shape wrt.
        //     // solid angle at 'it'. Returns pdf = r^2/(Area*cos??.)
        //     ds = m_shape->sample_direction(it, sample, active);
        //     active &= dot(ds.d, ds.n) < 0.f && neq(ds.pdf, 0.f);
        //
        //     SurfaceInteraction3f si(ds, it.wavelengths);
        //     geom_gain = m_antenna_texture->eval(si, active) / ds.pdf;
        // } else {
        //     // Importance sample the texture, then map onto the shape
        //     auto [uv, pdf] = m_antenna_texture->sample_position(sample, active);
        //     active &= neq(pdf, 0.f);
        //
        //     SurfaceInteraction3f si = m_shape->eval_parameterization(uv, active);
        //     si.wavelengths = it.wavelengths;
        //     active &= si.is_valid();
        //
        //     ds.p = si.p;
        //     ds.n = si.n;
        //     ds.uv = si.uv;
        //     ds.time = it.time;
        //     ds.delta = false;
        //     ds.d = ds.p - it.p;
        //
        //     Float dist_squared = squared_norm(ds.d);
        //     ds.dist = sqrt(dist_squared);
        //     ds.d /= ds.dist;
        //
        //     Float dp = dot(ds.d, ds.n);
        //     active &= dp < 0;
        //     ds.pdf = select(active, pdf / norm(cross(si.dp_du, si.dp_dv)) *
        //                                 dist_squared / -dp, 0.f);
        //
        //     geom_gain = m_antenna_texture->eval(si, active) / ds.pdf;
        // }
        // =============================
        // 3b. Directional part from WDF in [1/sr] --
        // ds.d *= -1.f;
        // DirectionSample3f ws = m_shape->sample_wigner(ds, it.wavelengths, active);
        // ds.d *= -1.f;
        // geom_gain *= ws.pdf;
        // ds.pdf *= ws.pdf;
        ds.d *= -1.f;
        Spectrum Wgain = sample_wigner(ds, it.wavelengths);
        ds.d *= -1.f;
        geom_gain *= Wgain;
        ds.pdf *= Wgain[0];
        ds.pdf = sqrt(ds.pdf * ds.pdf);
        // =============================

        // 4. Evaluate various extents to find ds radiance -----
        // Float extents = math::Pi<Float>;
        // Extent because mitsuba is inconsistent
        // Spectrum extents = rcp(m_shape->surface_area())*rcp(m_shape->surface_area())*math::Pi<float>;
        // Spectrum extents = rcp(m_shape->surface_area())*rcp(m_shape->surface_area())*rcp(m_shape->surface_area())*math::Pi<float>;
        // Spectrum extents = rcp(m_shape->surface_area())*rcp(m_shape->surface_area())*rcp(m_shape->surface_area()) * (ds.dist*ds.dist);
        // Maybe keep dsdist because this magically goes to tx, rather than naturally
        // Spectrum extents = rcp(m_shape->surface_area())*rcp(m_shape->surface_area())*rcp(m_shape->surface_area());
        // Spectrum extents = rcp(m_shape->surface_area())*rcp(m_shape->surface_area());
        // // Precaution:
        // Spectrum extents = rcp(m_shape->surface_area())*math::TwoPi<Float>;
        // Original:
        Float extents = 1.f;
        // extents *= math::Pi<float>;
        // ====================================================

        // std::cout << (signal_power*transmission_gain*geom_gain* extents) << std::endl;

        // std::cout << "Here direc" << std::endl;

        // 6. Return the ds, and ds radiant intensity ---------
        ds.object = this;
        // return { ds,
        //     unpolarized<Spectrum>(signal_power*transmission_gain*geom_gain*rcp(m_shape->surface_area())*rcp(m_shape->surface_area())) & active };
        // std::cout << " " << signal_power << " " << transmission_gain << " " << geom_gain << " " << signal_power*transmission_gain*geom_gain << std::endl;
        return { ds,
            unpolarized<Spectrum>(signal_power*transmission_gain*geom_gain*extents) & active };
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

        // std::cout << "tx pdf" << std::endl;

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
        // DirectionSample3f ws = m_shape->sample_wigner(ds_, it.wavelengths, active);
        // value *= ws.pdf;
        value *= sample_wigner(ds_, it.wavelengths)[0];
        value = sqrt(value*value);
        // =============================

        // 2. Evaluate various extents to find ds probability -----
        // // Precaution:
        // Float extents = math::Pi<Float>;
        // Original:
        Float extents = 1.f;
        // ====================================================

        // value = 1.f;

        return select(active, value*extents, 0.f);
        // return select(active, value, 0.f);
    }
    // ========================================================================

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("antenna_texture", m_antenna_texture.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "PhasedTransmitter[" << std::endl
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

    int m_n_elems;
    ScalarVector3f m_wid;
    ScalarFrame3f m_array_frame;

    std::vector<ScalarTransform4f> m_velem_to_world;
    std::vector<Transform4f> m_dir_to_local_velem;
    std::vector<ScalarTransform4f> m_velem_to_object;
    std::vector<ScalarFrame3f> m_velem_frame;
    std::vector<std::complex<ScalarFloat>> m_psi_dash;
    std::vector<ScalarVector3f> m_r_dash;
};

MTS_IMPLEMENT_CLASS_VARIANT(PhasedTransmitter, Transmitter)
MTS_EXPORT_PLUGIN(PhasedTransmitter, "Phased transmitter")
NAMESPACE_END(mitsuba)
