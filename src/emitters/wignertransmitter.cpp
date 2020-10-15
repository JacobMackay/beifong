#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

#include <unistd.h>

NAMESPACE_BEGIN(mitsuba)
/**!

.. _emitter-area:

Area light (:monosp:`area`)
---------------------------

.. pluginparameters::

 * - radiance
   - |spectrum|
   - Specifies the emitted radiance in units of power per unit area per unit steradian.

This plugin implements an area light, i.e. a light source that emits
diffuse illumination from the exterior of an arbitrary shape.
Since the emission profile of an area light is completely diffuse, it
has the same apparent brightness regardless of the observer's viewing
direction. Furthermore, since it occupies a nonzero amount of space, an
area light generally causes scene objects to cast soft shadows.

To create an area light source, simply instantiate the desired
emitter shape and specify an :monosp:`area` instance as its child:

.. code-block:: xml
    :name: sphere-light

    <shape type="sphere">
        <emitter type="area">
            <spectrum name="radiance" value="1.0"/>
        </emitter>
    </shape>

 */

template <typename Float, typename Spectrum>
class WignerTransmitter final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    WignerTransmitter(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        m_radiance = props.texture<Texture>("radiance", Texture::D65(1.f));

        m_flags = +EmitterFlags::Surface;
        if (m_radiance->is_spatially_varying())
            m_flags |= +EmitterFlags::SpatiallyVarying;
    }

    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

        // Ensures only forward emission
        return select(
            Frame3f::cos_theta(si.wi) > 0.f,
            unpolarized<Spectrum>(m_radiance->eval(si, active)),
            0.f
        );

        // return unpolarized<Spectrum>(m_radiance->eval(si, active));
        // return 0.f;
    }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &sample2, const Point2f &sample3,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;

        Float pdf = 1.f;

        // 1. Two strategies to sample spatial component based on 'm_radiance'
        if (!m_radiance->is_spatially_varying()) {
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
        return std::make_pair(
            Ray3f(ws.p, ws.d, ws.time, wavelength),
            select(abs(ws.pdf)>math::Epsilon<Float>, unpolarized<Spectrum>(spec_weight) / ws.pdf, 0.f)
        );
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area emitter without an associated Shape.");
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

        active &= abs(ws.pdf) > math::Epsilon<Float>;

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
    ref<Texture> m_radiance;
};

MTS_IMPLEMENT_CLASS_VARIANT(WignerTransmitter, Emitter)
MTS_EXPORT_PLUGIN(WignerTransmitter, "Wigner transmitter")
NAMESPACE_END(mitsuba)
