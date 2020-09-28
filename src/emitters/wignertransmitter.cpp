#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

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
class AreaLight final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    AreaLight(const Properties &props) : Base(props) {
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

        return select(
            Frame3f::cos_theta(si.wi) > 0.f,
            unpolarized<Spectrum>(m_radiance->eval(si, active)),
            0.f
        );
    }

    // DirectionSample3f sample_wigner(Float time, Float wavelength,
    //                                       const Point2f &sample2,
    //                                       const Point2f &sample3,
    //                                       Mask active) const override {
    //     MTS_MASK_ARGUMENT(active);
    //
    //     SurfaceInteraction3f si;
    //
    //     DirectionSample3f ds;
    //     ds.p = m_to_world.transform_affine(
    //         Point3f(sample2.x() * 2.f - 1.f, sample2.y() * 2.f - 1.f, 0.f));
    //     ps.n    = m_frame.n;
    //     ps.uv   = sample2;
    //     ps.time = time;
    //     ps.delta = false;
    //     ds.d = warp::square_to_cosine_hemisphere(sample3);
    //     // I forgot what dp actually was. Dot product! Ie how aligned is the
    //     // normal and query direction.
    //     ds.dist = 0;
    //     ds.obj = (const Object *) this;
    //
    //     ps.pdf  = m_inv_surface_area;m
    //     // pdf is abs wigner. pdf integrates to 1. Need to normalise.
    //
    //     // Static gain of rectangular transmitter
    //     Float g_static = 4*math::pi*m_surface_area/(wavelength*1e9)^2;
    //
    //     // For single rectangular element we can use sample space to cut some
    //     // corners. For more complex geometry, probably not.
    //
    //     // convert ds.d to
    //
    //     Float w_val = tri(sample2.x - 0.5) * sinc(2*nu*tri(sample2.x - 0.5))
    //
    //     return ds;
    // }

    DirectionSample3f sample_wigner(SurfaceInteraction3f si,
                                          const Point2f &sample2,
                                          Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // si is either a point in the scene, or
        // sample2 is used as position sample

        DirectionSample3f ds;
        // Convert sample2 to a world position
        ds.p = m_to_world.transform_affine(
            Point3f(sample2.x() * 2.f - 1.f, sample2.y() * 2.f - 1.f, 0.f));
        // The normal is the frame/object normal. May only be valid for flat
        // objects
        ps.n    = m_frame.n;
        // Optional,
        ps.uv   = sample2;
        ps.time = si.time;
        ps.delta = false;

        // ds.d = warp::square_to_cosine_hemisphere(sample3);
        // Should have a check if distance = or is close to 0
        ds.d = ds.p - it.p;
        Float dist_squared = squared_norm(ds.d);
        ds.dist = sqrt(dist_squared);
        ds.d /= ds.dist;
        // I forgot what dp actually was. Dot product! Ie how aligned is the
        // normal and query direction.
        // ds.dist = 0;
        ds.obj = (const Object *) this;

        // ps.pdf  = m_inv_surface_area;
        // pdf is abs wigner. pdf integrates to 1. Need to normalise.
        // No, the 'pdf' could be a qdf here. If it needs to be pdf later,
        // whatever


        // Static gain of rectangular transmitter

        float nu_0 = 1/si.wavelengths*1e9;



        Float g_static = 4*math::pi*m_surface_area*nu_0^2;

        // For single rectangular element we can use sample space to cut some
        // corners. For more complex geometry, probably not.

        // convert ds.d to
        // Float w_val = tri(m_to_local.transform_affine(ds.p))
            // * sinc(2*nu*tri(sample2.x - 0.5))

        // nuxz is the z component of the projection of the wavevector onto
        // plane zy plane

        // this should be the dot product of u axis and wavevector

        Float w_val = tri(sample2.x - 0.5) * sinc(2*ds.d**tri(sample2.x - 0.5))

        return ds;
    }

    // MTS_VARIANT typename Shape<Float, Spectrum>::DirectionSample3f
    // Shape<Float, Spectrum>::sample_direction(const Interaction3f &it,
    //                                          const Point2f &sample,
    //                                          Mask active) const {
    //     MTS_MASK_ARGUMENT(active);
    //
    //     DirectionSample3f ds(sample_position(it.time, sample, active));
    //     ds.d = ds.p - it.p;
    //
    //     Float dist_squared = squared_norm(ds.d);
    //     ds.dist = sqrt(dist_squared);
    //     ds.d /= ds.dist;
    //
    //     Float dp = abs_dot(ds.d, ds.n);
    //     ds.pdf *= select(neq(dp, 0.f), dist_squared / dp, 0.f);
    //     ds.object = (const Object *) this;
    //
    //     return ds;
    // }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &sample2, const Point2f &sample3,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // given rand pos/dir, sample the world loc, and return weight
        // Take in a shape (must be rectangle), this is the template shape.
        // Take in a world transform, this is the 'centre' point
        // Take in list of locations, these are the sub parts of the tx
        // Then expand a bbox, will prob need to instantiate a shape based on
        // this so we can get tinteractions.

        // BUT for now. Just have rectangular single one and get tx working.

        // jbo.D = sqrt(4*pi)*(jbo.l.zbw)*sqrt(4*pi)*(jbo.l.ybw)/lambda^2;

        // jbo.K = @(z,nuxz, y,nuxy) ...
        //     jbo.D .* ...
        //     tri((z-jbo.l.z0)/jbo.l.zbw,2).*...
        //     sinc(2*nuxz.*tri((z-jbo.l.z0 )/jbo.l.zbw,2).*jbo.l.zbw) .*...
        //     tri((y-jbo.l.y0)/jbo.l.ybw,2).*...
        //     sinc(2*nuxy.*tri((y-jbo.l.y0 )/jbo.l.ybw,2).*jbo.l.ybw); % [W/W]


        // return a position sample on the object, and its weight, ie
        // outgoing wigner vb

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

        return std::make_pair(
            Ray3f(si.p, si.to_world(local), time, wavelength),
            unpolarized<Spectrum>(spec_weight) * (math::Pi<Float> / pdf)
        );
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area emitter without an associated Shape.");
        DirectionSample3f ds;
        Spectrum spec;

        // Given a position in the scene and a random 2d sample, return a ds
        // and a spectrum.

        // Would like to specify m_wigner
        // Make it a function for now.
        // m_wigner->eval()
        // Given sample, active and si, return ds.

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

        spec = m_radiance->eval(si, active) / ds.pdf;
    }

    ds.object = this;
    return { ds, unpolarized<Spectrum>(spec) & active };

        // One of two very different strategies is used depending on 'm_radiance'
        if (!m_radiance->is_spatially_varying()) {
            // Texture is uniform, try to importance sample the shape wrt. solid angle at 'it'
            ds = m_shape->sample_direction(it, sample, active);
            active &= dot(ds.d, ds.n) < 0.f && neq(ds.pdf, 0.f);

            SurfaceInteraction3f si(ds, it.wavelengths);
            spec = m_radiance->eval(si, active) / ds.pdf;
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

            spec = m_radiance->eval(si, active) / ds.pdf;
        }

        ds.object = this;
        return { ds, unpolarized<Spectrum>(spec) & active };
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

        return select(active, value, 0.f);
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("radiance", m_radiance.get());
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "AreaLight[" << std::endl
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

MTS_IMPLEMENT_CLASS_VARIANT(AreaLight, Emitter)
MTS_EXPORT_PLUGIN(AreaLight, "Area emitter")
NAMESPACE_END(mitsuba)
