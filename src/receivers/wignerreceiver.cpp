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

               m_receive_type = props.string("signaltype", "raw");

               if (m_receive_type == "raw") {
                   m_sig_f_centre = props.float_("freq_centre", MTS_C/((MTS_WAVELENGTH_MAX + MTS_WAVELENGTH_MIN)/2));
                   m_sig_f_ext = props.float_("freq_ext", MTS_C/(MTS_WAVELENGTH_MAX - MTS_WAVELENGTH_MIN));
               } else if (m_receive_type == "mixer_dodgy") {
                   m_signal = props.string("signaltype", "cw");
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
    }

    // Spectrum eval_signal(Float time, Spectrum frequency) const {
    // Spectrum eval_signal(Float time, wavelength_t<Spectrum> frequency) const {
    Float eval_signal(Float time, Float frequency) const {

        Float result(0.f);
        // Spectrum result(0.f);
        Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
        Float t_hat;
        Float f_hat;
        // wavelength_t<Spectrum> f_hat;

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

    std::pair<wavelength_t<Spectrum>, Spectrum> sample_frequency(Float time, Float sample) const {
        auto freq_sample = math::sample_shifted<wavelength_t<Spectrum>>(sample);

        if (m_receive_type == "raw") {
            // Randomly select a freqency in bounds
            return {freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2),
                m_sig_f_ext};
        } else if (m_receive_type == "mixer_dodgy") {
            Float t_norm = math::fmodulo(time, rcp(m_sig_repfreq));
            Float frequencies = (m_sig_f_centre - m_sig_f_ext/2) + 0.5*m_sig_f_ext/m_sig_t_ext*t_norm;
            return {frequencies, eval_signal(time, frequencies)};
        } else {
            return {freq_sample * m_sig_f_ext + (m_sig_f_centre - m_sig_f_ext/2),
                m_sig_f_ext};
        }
    }



    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float wavelength_sample,
                            const Point2f & position_sample,
                            const Point2f & direction_sample,
                            Mask active) const override {

        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;
        Float pdf = 1.f;

        // 1. Sample spatial component ---------------
        PositionSample3f ps =
            m_shape->sample_position(time, position_sample, active);
        si = SurfaceInteraction3f(ps, zero<Wavelength>());
        pdf = ps.pdf;
        // ===========================================

        // 2. Sample directional component -----------
        Vector3f local = warp::square_to_cosine_hemisphere(direction_sample);
        // ===========================================

        // 3.a Sample frequency spectrum -------------
        auto [frequencies, freq_weight] = sample_frequency(time, wavelength_sample);

        wavelength_t<Spectrum> wavelengths = MTS_C/frequencies * 1e9;
        Spectrum wav_weight = MTS_C/freq_weight * 1e9; // Seems to be a sampling thing, but not sure exactly what
        // ===========================================

        // 3. Sample spectrum ------------------------
        // auto [wavelengths, wav_weight] =
        //     sample_wavelength<Float, Spectrum>(wavelength_sample);
        // ===========================================

        // Create direction sample for wigner --------
        DirectionSample3f ds;
        ds.p = si.p;
        ds.n = si.n;
        ds.uv = si.uv;
        ds.time = time;
        ds.delta = false;
        ds.d = local;
        Float dist_squared = squared_norm(ds.d);
        ds.dist = sqrt(dist_squared);
        ds.d /= ds.dist;

        // pdf here is the inverse surface area.
        ds.pdf = pdf;
        // // Assuming wigner doesn't take look angle into account
        // Float dp = abs_dot(ds.d, ds.n);
        // // pdf is now 1/A * r^2/cos(θ)cos(φ)
        // ds.pdf *= select(neq(dp, 0.f), dist_squared / dp, 0.f);

        ds.object = this;
        // ===========================================

        // also is done by area alrteady
        // but probably not got the norm from radiance?
        DirectionSample3f ws = m_shape->sample_wigner(ds, wavelengths, active);
        // The function should be divided by λ^2, pdf multiplied.
        // ws.pdf *= (wavelengths[0]*1e-9)*(wavelengths[0]*1e-9);
        // ws.pdf *= 2;

        // return std::make_pair(
        //     RayDifferential3f(ws.p, Frame3f(ps.n).to_world(local), time,
        //         wavelengths),
        //         unpolarized<Spectrum>(wav_weight)*math::Pi<ScalarFloat>*(rcp(ws.pdf))
        // );

        // return std::make_pair(
        //     RayDifferential3f(ws.p, Frame3f(ps.n).to_world(local), time,
        //         wavelengths),
        //         unpolarized<Spectrum>(wav_weight)
        //         * m_shape->surface_area()/(4*math::Pi<ScalarFloat>)*(rcp(ws.pdf))
        // );
        // return std::make_pair(
        //     RayDifferential3f(ws.p, Frame3f(ps.n).to_world(local), time,
        //         wavelengths),
        //         unpolarized<Spectrum>(wav_weight)
        //         * 1/(4*math::Pi<ScalarFloat>)*(rcp(ws.pdf)*math::Pi<ScalarFloat>)
        // );

        // std::cout << "samp_rdiff: " << unpolarized<Spectrum>(wav_weight)
        // *(rcp(ws.pdf)) << std::endl;

        return std::make_pair(
            RayDifferential3f(ws.p, Frame3f(ps.n).to_world(local), time,
                wavelengths),
                unpolarized<Spectrum>(wav_weight)
                *(rcp(ws.pdf))
        );
        // return std::make_pair(
        //     RayDifferential3f(ws.p, Frame3f(ps.n).to_world(local), time,
        //         wavelengths),
        //         unpolarized<Spectrum>(wav_weight)
        //         * m_shape->surface_area()/(4*math::Pi<ScalarFloat>))
        // );

        // return std::make_pair(
        //     RayDifferential3f(ps.p, Frame3f(ps.n).to_world(local), time,
        //         wavelengths),
        //     unpolarized<Spectrum>(wav_weight)
        //         * math::Pi<ScalarFloat> / m_shape->surface_area()
        // );
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample,
        Mask active) const override {

        DirectionSample3f ds = m_shape->sample_direction(it, sample, active);
        ds.d *= -1.f;
        DirectionSample3f ws = m_shape->sample_wigner(ds, it.wavelengths, active);
        ws.d *= -1.f;
        Spectrum spec = rcp(ws.pdf);

        // std::cout << "sample_direction" << std::endl;

        // return std::make_pair(ws, math::Pi<ScalarFloat>*spec);

        // return std::make_pair(ws, 1.f/(4*math::Pi<ScalarFloat>));
        // return std::make_pair(ws, spec*1.f/(4*math::Pi<ScalarFloat>)*math::Pi<ScalarFloat>);
        return std::make_pair(ws, spec*math::Pi<ScalarFloat>);

        // return std::make_pair(m_shape->sample_direction(it, sample, active),
        //     math::Pi<ScalarFloat>);
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        // std::cout << "pdf_direction" << std::endl;
        return m_shape->pdf_direction(it, ds, active);
    }

    // I'm guessing pi comes from integrating over directions, cosine hemisphere
    Spectrum eval(const SurfaceInteraction3f &/*si*/,
        Mask /*active*/) const override {
        // std::cout << "eval" << std::endl;
        // return math::Pi<ScalarFloat> / m_shape->surface_area();
        // return math::Pi<ScalarFloat>;
        // return m_shape->surface_area()/(4*math::Pi<ScalarFloat>)*math::Pi<ScalarFloat>;
        // std::cout << "eval" << std:: endl;
        return m_shape->surface_area()*math::Pi<ScalarFloat>;
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
