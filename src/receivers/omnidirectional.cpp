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

// ==============================
// Note, this receiver reports back in W/m^2
// Changed to give back W
// Hypothesis: This says it measures the power coming in over area
// I want the power in total, but what is coming in is in w/m^2/sr
// I need to multiply by A
// ==============================

MTS_VARIANT class Omnidirectional final : public Receiver<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Receiver, m_adc, m_world_transform, m_shape)
    MTS_IMPORT_TYPES(Shape)

    Omnidirectional(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The omnidirectional receiver inherits this transformation "
                  "from its parent shape.");

        // if (m_adc->size() != ScalarPoint2i(1, 1))
        //     Throw("This sensor only supports films of size 1x1 Pixels!");

        if (m_adc->reconstruction_filter()->radius() >
            0.5f + math::RayEpsilon<Float>)
            Log(Warn, "This sensor should only be used with a reconstruction filter"
               "of radius 0.5 or lower(e.g. default box)");
    }

    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float wavelength_sample,
                            const Point2f & position_sample,
                            const Point2f & direction_sample,
                            Mask active) const override {

        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Sample spatial component
        PositionSample3f ps =
            m_shape->sample_position(time, position_sample, active);

        // 2. Sample directional component
        Vector3f local = warp::square_to_cosine_hemisphere(direction_sample);

        // 3. Sample spectrum
        auto [wavelengths, wav_weight] =
            sample_wavelength<Float, Spectrum>(wavelength_sample);

        // std::cout << "samp_rdiff: " << unpolarized<Spectrum>(wav_weight)
        //     * m_shape->surface_area() << std::endl;

        return std::make_pair(
            RayDifferential3f(ps.p, Frame3f(ps.n).to_world(local), time,
                wavelengths),
            // unpolarized<Spectrum>(wav_weight)
            //     * math::Pi<ScalarFloat> / m_shape->surface_area()
            // unpolarized<Spectrum>(wav_weight)
            //     * math::Pi<ScalarFloat>
            // unpolarized<Spectrum>(wav_weight)
            //     * m_shape->surface_area()/(math::Pi<ScalarFloat>/2)
            // unpolarized<Spectrum>(wav_weight)
            //     * m_shape->surface_area()/(4*math::Pi<ScalarFloat>)
            unpolarized<Spectrum>(wav_weight)
                * m_shape->surface_area()
        );
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample,
        Mask active) const override {
        // return std::make_pair(m_shape->sample_direction(it, sample, active),
        //     math::Pi<ScalarFloat>);
        // return std::make_pair(m_shape->sample_direction(it, sample, active),
        //     1.f/(math::Pi<ScalarFloat>/2));
        // return std::make_pair(m_shape->sample_direction(it, sample, active),
        //     1.f/(4*math::Pi<ScalarFloat>));

        // std::cout << "samp_direc: " << std::endl;

        return std::make_pair(m_shape->sample_direction(it, sample, active),
            1.f);
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        return m_shape->pdf_direction(it, ds, active);
    }

    Spectrum eval(const SurfaceInteraction3f &/*si*/,
        Mask /*active*/) const override {
        // return math::Pi<ScalarFloat> / m_shape->surface_area();
        // return math::Pi<ScalarFloat>;
        // return m_shape->surface_area()/(math::Pi<ScalarFloat>/2);
        // return m_shape->surface_area()/(4*math::Pi<ScalarFloat>);
        return m_shape->surface_area();
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Omnidirectional[" << std::endl
            << "  shape = " << m_shape << "," << std::endl
            << "  ADC = " << m_adc << "," << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_VARIANT(Omnidirectional, Receiver)
MTS_EXPORT_PLUGIN(Omnidirectional, "Omnidirectional");
NAMESPACE_END(mitsuba)
