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
class Phased final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
    MTS_IMPORT_TYPES(Scene, Shape, Texture)

    Phased(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        // Antenna characteristics
        m_antenna_texture = props.texture<Texture>("antenna_texture", Texture::D65(1.f));
        // Signal characteristics
        m_power = props.float_("power", 1.f);
        m_gain = props.float_("gain", 1.f);

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

        m_flags = +EmitterFlags::Surface;
        if (m_antenna_texture->is_spatially_varying())
            m_flags |= +EmitterFlags::SpatiallyVarying;
    }

    Spectrum W_rect_2D(const Point3f &r_hat,
                        const Normal3f &nu_hat,
                        const ScalarVector3f &wid) const {

        return 4*wid.x()*wid.y() * math::tri(r_hat.x())*math::tri(r_hat.y()) *
                math::sinc(math::TwoPi<ScalarFloat>*nu_hat.x()*wid.x()*math::tri(r_hat.x())) *
                math::sinc(math::TwoPi<ScalarFloat>*nu_hat.y()*wid.y()*math::tri(r_hat.y()));
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
        // DirectionSample3f ws = m_shape->sample_wigner(ds, si.wavelengths, active);
        // geom_gain *= ws.pdf;
        geom_gain *= sample_wigner(ds, si.wavelengths);



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
        // DirectionSample3f ws =
        //     m_shape->sample_wigner(ds, wavelength, active);
        // Spectrum geom_gain = ws.pdf * pdf;
        Spectrum geom_gain = sample_wigner(ds, wavelength);
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

            geom_gain = m_antenna_texture->eval(si, active) / ds.pdf * dp * dp;
        }
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
        // DirectionSample3f ws = m_shape->sample_wigner(ds_, it.wavelengths, active);
        // value *= ws.pdf;
        value *= sample_wigner(ds_, it.wavelengths)[0];
        value = sqrt(value*value);
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
        oss << "Phased[" << std::endl
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
    ScalarFloat m_gain;
    ScalarFloat m_power;
    ref<Texture> m_antenna_texture;

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

MTS_IMPLEMENT_CLASS_VARIANT(Phased, Emitter)
MTS_EXPORT_PLUGIN(Phased, "Phased")
NAMESPACE_END(mitsuba)
