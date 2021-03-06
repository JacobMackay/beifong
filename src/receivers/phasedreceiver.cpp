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

// template <typename Float, typename Spectrum>
// template <typename Float, typename Spectrum>
MTS_VARIANT class Phasedreceiver final : public Receiver<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Receiver, m_adc, m_world_transform, m_shape, m_receive_type)
    MTS_IMPORT_TYPES(Shape)

    Phasedreceiver(const Properties &props) : Base(props) {
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

        // Signal Components --------------------------------------------------
       if (m_receive_type == "raw") {
           m_sig_f_centre = props.float_("freq_centre", 1.f);
           m_sig_f_ext = props.float_("freq_ext", 1.f);
           m_gain = props.float_("gain", 1.f);
       } else if (m_receive_type == "raw_resample") {
           m_sig_f_centre = props.float_("freq_centre", 1.f);
           m_sig_f_ext = props.float_("freq_ext", 1.f);
           m_gain = props.float_("gain", 1.f);
       } else if (m_receive_type == "mix_resample") {
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
    Spectrum eval_signal(Float time, Wavelength frequency) const {

        Spectrum result(0.f);
        Float t;
        Float ti;
        Wavelength fi;

        if (m_signal == "linfmcw") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            fi = m_sig_f_centre + (m_sig_f_ext/m_sig_t_ext)*(t - ti);
            result = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - fi[0], m_sig_t_ext, m_sig_amplitude),
                0.f);
        } else if (m_signal == "pulse") {
            t = math::fmodulo(time, rcp(m_sig_repfreq));
            ti = 0 + m_sig_t_ext/2;
            fi = m_sig_f_centre;
            result = select(math::rect((t - ti)/m_sig_t_ext) > 0.f,
                math::wchirp(t - ti, frequency[0] - fi[0], m_sig_t_ext, m_sig_amplitude),
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
            Float t = math::fmodulo(time, rcp(m_sig_repfreq));
            Float ti = 0 + m_sig_t_ext/2;
            Wavelength fi = m_sig_f_centre + (m_sig_f_ext/m_sig_t_ext)*(t - ti);
            frequencies = fi;
        } else if (m_signal == "cw") {
            frequencies = m_sig_f_centre;
        }

        // return {frequencies, eval_signal(time, frequencies)};
        return {frequencies, 1.f};
        // ===============================================
    }
    // ========================================================================

    // Return a frequency sample from signal at time t
    // in addition to the spectral/signal weight/power in V^2/Hz
    // ------------------------------------------------------------------------
    std::pair<Wavelength, Spectrum>
    sample_frequency(Float time, Float frequency_sample) const {
        std::pair<Wavelength, Spectrum> result;

        if (m_receive_type == "raw" | m_receive_type == "raw_resample"){
            auto freq_sample = math::sample_shifted<Wavelength>(frequency_sample);
            result.first =
                freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
            result.second = 1.f;
        } else if (m_receive_type == "mix_resample") {
            // Probably include phase here, but ignore for now. Eg tellurometer
            if (m_sig_is_delta == true) {
                result = sample_delta_frequency(time);
            } else {
                auto freq_sample = math::sample_shifted<Wavelength>(frequency_sample);
                result.first =
                    freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
                result.second = eval_signal(time, result.first);
                // result.second = 1.f;
            }
        } else {
            if (m_sig_is_delta == true) {
                result = sample_delta_frequency(time);
            } else {
                auto freq_sample = math::sample_shifted<Wavelength>(frequency_sample);
                result.first =
                    freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2);
                result.second = eval_signal(time, result.first);
            }
        }
        return result;
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
        // and radiation resistance in [1/??] -------------------
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
        // DirectionSample3f ws =
        //     m_shape->sample_wigner(ds, wavelength, active);
        // Spectrum geom_gain = ws.pdf * pdf;
        Spectrum geom_gain = sample_wigner(ds, wavelength) * pdf * (1 - dot(ds.d, ds.n)* dot(ds.d, ds.n)* dot(ds.d, ds.n)* dot(ds.d, ds.n));
        // ===========================================
        // ====================================================

        // 5. Evaluate various extents
        // to find individual ray power -----------------------
        // Spectrum extents = m_shape->surface_area() * math::Pi<Float>;
        Spectrum extents = m_shape->surface_area() * math::Pi<Float>;
        if (!m_sig_is_delta) {
            extents *= MTS_C*rcp(m_sig_f_ext)*1e9;
        }
        // ====================================================

        Float phase = 0;

        // 6. Return the ray, and ray power -------------------
        return std::make_pair(
            Ray3f(si.p, si.to_world(ds.d), time, phase, wavelength),
            unpolarized<Spectrum>(signal_power * reception_gain * geom_gain * extents)
            // unpolarized<Spectrum>(signal_power * reception_gain * geom_gain)
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
        oss << "Phasedreceiver[" << std::endl
            << "  shape = " << m_shape << "," << std::endl
            << "  ADC = " << m_adc << "," << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    Float m_gain;
    bool m_resample_freq;
    std::string m_signal;
    Float m_sig_amplitude;
    Float m_sig_repfreq;
    Float m_sig_t_ext;
    Float m_sig_f_centre;
    Float m_sig_f_ext;
    Float m_sig_phi0;
    bool m_sig_is_delta;
    // ref<Texture> m_antenna_texture;

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

MTS_IMPLEMENT_CLASS_VARIANT(Phasedreceiver, Receiver)
MTS_EXPORT_PLUGIN(Phasedreceiver, "Phasedreceiver");
NAMESPACE_END(mitsuba)
