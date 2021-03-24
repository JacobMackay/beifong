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
class Coherent final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    Coherent(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        // Antenna characteristics
        m_antenna_texture = props.texture<Texture>("antenna_texture", Texture::D65(1.f));
        // Signal characteristics
        m_power = props.float_("power", 1.f);
        m_gain = props.float_("gain", 1.f);

        m_flags = +EmitterFlags::Surface;
        if (m_antenna_texture->is_spatially_varying())
            m_flags |= +EmitterFlags::SpatiallyVarying;
    }

    // Return the spectral radiance of an impacting ray in [W/(sr*m^2*sm)]
    // ------------------------------------------------------------------------
    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

        // 1. Evaluate the signal power in [V^2/Hz] -----------------
        Spectrum signal_power = m_power;

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
        return select(
            Frame3f::cos_theta(si.wi) > 0.f,
            signal_power * transmission_gain * geom_gain,
            0.f
        );
        // ====================================================
    }

    // Sample a ray and return the power
    // emanating at a position, in a direction and with a wavelength
    // W/sr * π
    // ------------------------------------------------------------------------
    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &position_sample,
                                          const Point2f &direction_sample,
                                          Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Evaluate the signal power in [V^2/Hz] ------------
        Spectrum signal_power = m_power;
        // =====================================================

        // 2. Evaluate the line loss, amplifier gain
        // and radiation resistance in [1/Ω] -------------------
        Spectrum transmission_gain = m_gain;
        // =====================================================

        // 3. Sample wavelengths
        // -----------------------------------------------------
        Wavelength wavelength;
        Spectrum spec_weight;

        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = math::Infinity<Float>;
        Float pdf = 1.f;

        if constexpr (is_spectral_v<Spectrum>) {
            std::tie(wavelength, spec_weight) = m_antenna_texture->sample_spectrum(
                si, math::sample_shifted<Wavelength>(wavelength_sample), active);
        } else {
            wavelength = zero<Wavelength>();
            spec_weight = m_antenna_texture->eval(si, active);
        }
        // =====================================================

        // 4. Evaluate the geometric gain in [1/(sr*m^2)] ------
        // 4a. Positional part from surface in [1/m^2] --
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
            m_shape->sample_wigner(ds, wavelength, active);
        Spectrum geom_gain = ws.pdf * pdf;
        // ===========================================
        // ====================================================

        // 5. Evaluate various extents
        // to find individual ray power -----------------------
        Spectrum extents = m_shape->surface_area() * math::Pi<Float>;
        // ====================================================

        // 6. Return the ray, and ray power -------------------
        return std::make_pair(
            Ray3f(si.p, si.to_world(ds.d), time, wavelength),
            unpolarized<Spectrum>(signal_power * transmission_gain * geom_gain * extents * spec_weight)
        );
        // ====================================================
    }

    // W/sr *cos(θ)/r^2
    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area emitter without an associated Shape.");

        // 1. Evaluate the signal power in [V^2/Hz] ------------
        Spectrum signal_power = m_power;
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
        // Float extents = math::Pi<Float>;
        Float extents = 1.f;
        // ====================================================

        // 5. Return the ds, and ds radiant intensity ---------
        ds.object = this;
        return { ds,
            unpolarized<Spectrum>((signal_power*transmission_gain*geom_gain)
         * extents) & active };
        // ====================================================
    }

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
        DirectionSample3f ws = m_shape->sample_wigner(ds_, it.wavelengths, active);
        value *= ws.pdf;
        // =============================

        // 2. Evaluate various extents to find ds probability -----
        // Float extents = math::Pi<Float>;
        Float extents = 1.f;
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
        oss << "Coherent[" << std::endl
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
    Float m_power;
    ref<Texture> m_antenna_texture;
};

MTS_IMPLEMENT_CLASS_VARIANT(Coherent, Emitter)
MTS_EXPORT_PLUGIN(Coherent, "Coherent")
NAMESPACE_END(mitsuba)
