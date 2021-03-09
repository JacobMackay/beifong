#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/util.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/shape.h>

#if defined(MTS_ENABLE_OPTIX)
    #include "optix/rectangle.cuh"
#endif

NAMESPACE_BEGIN(mitsuba)

/**!

.. _shape-rectangle:

Rectangle (:monosp:`rectangle`)
-------------------------------------------------

.. pluginparameters::

 * - flip_normals
   - |bool|
   - Is the rectangle inverted, i.e. should the normal vectors be flipped? (Default: |false|)
 * - to_world
   - |transform|
   - Specifies a linear object-to-world transformation. (Default: none (i.e. object space = world space))

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/shape_rectangle.jpg
   :caption: Basic example
.. subfigure:: ../../resources/data/docs/images/render/shape_rectangle_parameterization.jpg
   :caption: A textured rectangle with the default parameterization
.. subfigend::
   :label: fig-rectangle

This shape plugin describes a simple rectangular shape primitive.
It is mainly provided as a convenience for those cases when creating and
loading an external mesh with two triangles is simply too tedious, e.g.
when an area light source or a simple ground plane are needed.
By default, the rectangle covers the XY-range :math:`[-1,1]\times[-1,1]`
and has a surface normal that points into the positive Z-direction.
To change the rectangle scale, rotation, or translation, use the
:monosp:`to_world` parameter.


The following XML snippet showcases a simple example of a textured rectangle:

.. code-block:: xml

    <shape type="rectangle">
        <bsdf type="diffuse">
            <texture name="reflectance" type="checkerboard">
                <transform name="to_uv">
                    <scale x="5" y="5" />
                </transform>
            </texture>
        </bsdf>
    </shape>
 */

template <typename Float, typename Spectrum>
class Rectangle final : public Shape<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Shape, m_to_world, m_velocity, m_to_object, set_children,
                    get_children_string, parameters_grad_enabled)
    MTS_IMPORT_TYPES()

    using typename Base::ScalarSize;

    Rectangle(const Properties &props) : Base(props) {
        if (props.bool_("flip_normals", false))
            m_to_world = m_to_world * ScalarTransform4f::scale(ScalarVector3f(1.f, 1.f, -1.f));

        update();
        set_children();
    }

    void update() {
        m_to_object = m_to_world.inverse();

        ScalarVector3f dp_du = m_to_world * ScalarVector3f(2.f, 0.f, 0.f);
        ScalarVector3f dp_dv = m_to_world * ScalarVector3f(0.f, 2.f, 0.f);
        ScalarNormal3f normal = normalize(m_to_world * ScalarNormal3f(0.f, 0.f, 1.f));
        m_frame = ScalarFrame3f(dp_du, dp_dv, normal);

        m_inv_surface_area = rcp(surface_area());
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarBoundingBox3f bbox;
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-1.f, -1.f, 0.f)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f( 1.f, -1.f, 0.f)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f( 1.f,  1.f, 0.f)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-1.f,  1.f, 0.f)));
        return bbox;
    }

    ScalarFloat surface_area() const override {
        return norm(cross(m_frame.s, m_frame.t));
    }

    // UnpolarizedSpectrum fourier_weight(const Vector3f &local_vect,
    //     Wavelength wavelengths) const override {
    //         // Really should be wavevectors
    //         // 2.0f / (I * q_parr^2) * q_x
    //         //
    //         // std::cout << props << std::endl;
    //
    //         // float width_x = m_frame.s;
    //         // float width_y = m_frame.t;
    //
    //         //
    //         // std::cout << m_frame.s << m_frame.t << m_frame.n << std::endl;
    //
    //         // Coord system in world is like x, y, z
    //
    //         // It would be good to be able to provide array of vectors and
    //         // evaluate them in parallel
    //         // Return should be a list of n wavelengths wide, m localvectors
    //         // long weights
    //
    //         // Wavevector q
    //         // if is constexpr spectral
    //         Array<Vector3f, 4> q = local_vect;
    //         for (int i = 0; i < 4; ++i) {
    //             q[i] /= wavelengths[i];
    //         }
    //
    //         Array<Vector3f, 5> vertices;
    //         vertices[0] = m_frame.to_local(
    //             m_to_world.transform_affine(ScalarPoint3f(-1.f, -1.f, 0.f)));
    //         vertices[1] = m_frame.to_local(
    //             m_to_world.transform_affine(ScalarPoint3f(1.f, -1.f, 0.f)));
    //         vertices[2] = m_frame.to_local(
    //             m_to_world.transform_affine(ScalarPoint3f(-1.f, 1.f, 0.f)));
    //         vertices[3] = m_frame.to_local(
    //             m_to_world.transform_affine(ScalarPoint3f(-1.f, 1.f, 0.f)));
    //         vertices[4] = vertices[0];
    //
    //         Array<Vector3f, 4> E_vert;
    //         Array<Vector3f, 4> R_vert;
    //
    //         // Create the simple polygonal vertex chain
    //         for (int i = 1; i < 5; ++i) {
    //             E_vert[0] = (vertices[i] - vertices[i-1])/2.f;  // diff
    //             R_vert[0] = (vertices[i] + vertices[i-1])/2.f;  // mean
    //         }
    //
    //         // Find the components of the wavevector compared to the surface
    //         Vector3f q_perp = dot(q, m_frame.n)*Vector3f(m_frame.n);
    //         Vector3f q_parr = q - q_perp;
    //         Vector3f q_crss = cross(m_frame.n, q_parr);
    //
    //         // Array<Vector3f, 4> E_vert;
    //         //
    //         // UnpolarizedSpectrum result;
    //
    //
    //
    //
    //         // std::cout << q_perp << std::endl;
    //
    //         // Vector3f vertices =
    //
    //         return 1.f;
    // }

    // =============================================================
    //! @{ \name Sampling routines
    // =============================================================

    PositionSample3f sample_position(Float time, const Point2f &sample,
                                     Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        PositionSample3f ps;
        ps.p = m_to_world.transform_affine(
            Point3f(sample.x() * 2.f - 1.f, sample.y() * 2.f - 1.f, 0.f));
        ps.n    = m_frame.n;
        ps.pdf  = m_inv_surface_area;
        ps.uv   = sample;
        ps.time = time;
        ps.delta = false;

        return ps;
    }

    Float pdf_position(const PositionSample3f & /*ps*/, Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        return m_inv_surface_area;
    }

    DirectionSample3f sample_wigner(const DirectionSample3f &ds,
                                    Wavelength wavelength,
                                    Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // ---------------------------------
        // Get the size of the target
        Float wid_x = norm(m_frame.s);
        Float wid_y = norm(m_frame.t);
        // =================================

        // ---------------------------------
        // Get the local position vector, normalised by width in x and y
        Point3f r_hat = m_to_object * ds.p/2;
        // =================================

        // ---------------------------------
        // Get the local wavevector
        Transform4f trafto1;
        trafto1 = trafto1.from_frame(
            Frame3f(normalize(m_frame.s),
                    normalize(m_frame.t),
                    normalize(m_frame.n)));
        // Normal3f nu_hat = trafto1 * ws.d;
        Normal3f nu_hat = trafto1.transform_affine(ds.d)*rcp(wavelength[0]*1e-9);
        // =================================

        // The wigner provides a function of 1/(sr*sm)
        // Intuitively: If the object gets bigger, the directional gain at 0
        // gets higher, as does the wdf. If the object gets smaller, the gain
        // at 0 is smaller, ie there is less power flowing through 1 sr at 0.
        // The WDF provides the directional gain as a function of position.
        // ---------------------------------
        // Find the wigner function value
        // Float gain = (math::TwoPi<Float>*rcp(wavelength[0]*1e-9))*(math::TwoPi<Float>*rcp(wavelength[0]*1e-9))*
        //         4*wid_x*wid_y * math::tri(r_hat.x())*math::tri(r_hat.y()) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.x()*wid_x*math::tri(r_hat.x())) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.y()*wid_y*math::tri(r_hat.y()));
        // Float gain = 4*wid_x*wid_y * math::tri(r_hat.x())*math::tri(r_hat.y()) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.x()*wid_x*math::tri(r_hat.x())) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.y()*wid_y*math::tri(r_hat.y()));

        // Float gain = 4*math::TwoPi<Float>*math::TwoPi<Float>* math::tri(r_hat.x())*math::tri(r_hat.y()) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.x()*wid_x*math::tri(r_hat.x())) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.y()*wid_y*math::tri(r_hat.y()));

        Float gain = 4*wid_x*wid_y * math::tri(r_hat.x())*math::tri(r_hat.y()) *
                math::sinc(math::TwoPi<Float>*nu_hat.x()*wid_x*math::tri(r_hat.x())) *
                math::sinc(math::TwoPi<Float>*nu_hat.y()*wid_y*math::tri(r_hat.y()));

        // Float gain = math::TwoPi<Float>* 4*math::tri(r_hat.x())*math::tri(r_hat.y()) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.x()*wid_x*math::tri(r_hat.x())) *
        //         math::sinc(math::TwoPi<Float>*nu_hat.y()*wid_y*math::tri(r_hat.y()));
        // =================================

        // gain *= math::InvTwoPi<Float>*2;

        // ---------------------------------
        // Add the wavelength normalisation that should be with radiance
        // gain *= (wavelength[0]*1e-9)*(wavelength[0]*1e-9);
        // It seems that they already divide by area, so let's undo our area normalisation
        // gain *= wid_x*wid_y;
        // =================================

        // This probably doesn't need to be a ds. A spectrum is probably more
        // approipriate

        // Maybe already multiplied by area
        DirectionSample3f ws = ds;
        // Will eventually remove this, but leave for now.
        // ws.pdf *= rcp(gain);
        ws.pdf = gain; // This follows the shape->sample position where pdf is 1/Area. wdf.pdf is 1/sr
        // ws.pdf = rcp(gain);

        return ws;
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

    PreliminaryIntersection3f ray_intersect_preliminary(const Ray3f &ray_,
                                                        Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        Ray3f ray     = m_to_object.transform_affine(ray_);
        Float t       = -ray.o.z() * ray.d_rcp.z();
        Point3f local = ray(t);

        // Is intersection within ray segment and rectangle?
        active = active && t >= ray.mint
                        && t <= ray.maxt
                        && abs(local.x()) <= 1.f
                        && abs(local.y()) <= 1.f;

        PreliminaryIntersection3f pi = zero<PreliminaryIntersection3f>();
        pi.t = select(active, t, math::Infinity<Float>);
        pi.prim_uv = Point2f(local.x(), local.y());
        pi.shape = this;

        return pi;
    }

    Mask ray_test(const Ray3f &ray_, Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        Ray3f ray     = m_to_object.transform_affine(ray_);
        Float t       = -ray.o.z() * ray.d_rcp.z();
        Point3f local = ray(t);

        // Is intersection within ray segment and rectangle?
        return active && t >= ray.mint
                      && t <= ray.maxt
                      && abs(local.x()) <= 1.f
                      && abs(local.y()) <= 1.f;
    }

    SurfaceInteraction3f compute_surface_interaction(const Ray3f &ray,
                                                     PreliminaryIntersection3f pi,
                                                     HitComputeFlags flags,
                                                     Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        bool differentiable = false;
        if constexpr (is_diff_array_v<Float>)
            differentiable = requires_gradient(ray.o) ||
                             requires_gradient(ray.d) ||
                             parameters_grad_enabled();

        // Recompute ray intersection to get differentiable prim_uv and t
        if (differentiable && !has_flag(flags, HitComputeFlags::NonDifferentiable))
            pi = ray_intersect_preliminary(ray, active);

        active &= pi.is_valid();

        SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
        si.t = select(active, pi.t, math::Infinity<Float>);

        si.p = ray(pi.t);

        si.n          = m_frame.n;
        si.sh_frame.n = m_frame.n;
        si.dp_du      = m_frame.s;
        si.dp_dv      = m_frame.t;
        si.uv         = Point2f(fmadd(pi.prim_uv.x(), .5f, .5f),
                                fmadd(pi.prim_uv.y(), .5f, .5f));

        si.dn_du = si.dn_dv = zero<Vector3f>();

        return si;
    }

    // Transform4f velocity() const override {
    //     Transform4f vel = m_velocity;
    //     // Transform4f vel = m_to_world;
    //     return vel;
    // }

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/) override {
        update();
        Base::parameters_changed();
#if defined(MTS_ENABLE_OPTIX)
        optix_prepare_geometry();
#endif
    }

#if defined(MTS_ENABLE_OPTIX)
    using Base::m_optix_data_ptr;

    void optix_prepare_geometry() override {
        if constexpr (is_cuda_array_v<Float>) {
            if (!m_optix_data_ptr)
                m_optix_data_ptr = cuda_malloc(sizeof(OptixRectangleData));

            OptixRectangleData data = { bbox(), m_to_object, m_frame.n,
                                        m_frame.s, m_frame.t };

            cuda_memcpy_to_device(m_optix_data_ptr, &data, sizeof(OptixRectangleData));
        }
    }
#endif

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Rectangle[" << std::endl
            << "  to_world = " << string::indent(m_to_world, 13) << "," << std::endl
            << "  velocity = " << string::indent(m_velocity, 13) << "," << std::endl
            << "  frame = " << string::indent(m_frame) << "," << std::endl
            << "  surface_area = " << surface_area() << "," << std::endl
            << "  " << string::indent(get_children_string()) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ScalarFrame3f m_frame;
    ScalarFloat m_inv_surface_area;
};

MTS_IMPLEMENT_CLASS_VARIANT(Rectangle, Shape)
MTS_EXPORT_PLUGIN(Rectangle, "Rectangle intersection primitive");
NAMESPACE_END(mitsuba)
