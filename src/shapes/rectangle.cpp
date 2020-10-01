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
    MTS_IMPORT_BASE(Shape, m_to_world, m_to_object, set_children,
                    get_children_string, parameters_grad_enabled)
    MTS_IMPORT_TYPES()

    using typename Base::ScalarSize;

    Rectangle(const Properties &props) : Base(props) {
        if (props.bool_("flip_normals", false))
            m_to_world = m_to_world * ScalarTransform4f::scale(ScalarVector3f(1.f, 1.f, -1.f));

        update();
        set_children();
        // std::cout << props << std::endl;
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

    // /// Sinc function basic. Sin(x)/x
    // template <typename T, typename Value = expr_t<T>>
    // Value sinc(const T &x) {
    //     return select(x <= math::Epsilon<T>, sin(x)/x, 1.f);
    // }
    //
    // /// Triangular function, base length 1.
    // template <typename T, typename Value = expr_t<T>>
    // Value tri(const T &x) {
    //     return select(abs(x) <= 0.5, 1.0 - 2.0*abs(x), x*0.f);
    // }

    // /// Sinc function basic. Sin(x)/x
    // Float sinc(Float x) const {
    //     return select(x <= math::Epsilon<Float>, sin(x)/x, 1.f);
    // }
    //
    // /// Triangular function, base length 1.
    // Float tri(Float x) const {
    //     return select(abs(x) < 0.5, 1.0 - 2.0*abs(x), x*0.f);
    // }

    // DirectionSample3f sample_wigner(Interaction3f it,
    //                                       const Point2f &sample,
    //                                       Mask active) const override {
    DirectionSample3f sample_wigner(const DirectionSample3f &ds, Wavelength wavelength,
                                          Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // We have world position, direction and wavelength. Now calculate
        // weight and store in pdf

        // This probably doesn't need to be a ds. A spectrum is probably more
        // approipriate

        DirectionSample3f ws = ds;

        // These values are global.
        // Vector3f nu_hat = ws.d;
        Float nu_0 = rcp(wavelength[0]*1e-9);
        // Vector3f k = nu_0*nu_hat;
        Float k_0 = 2*math::Pi<Float>*nu_0;

        // Make local
        // Vector3f nu_hat_local = m_to_object.transform_affine(ws.d);

        Vector3f nu_hat_local = ws.d;

        // Vector3f nu_hat_local = m_frame.to_local(ws.d);
        Float k_x = k_0*nu_hat_local.x();
        Float k_y = k_0*nu_hat_local.y();
        // Point3f p_local = m_to_object.transform_affine(ws.p);
        // Vector3f p_local = m_frame.to_local(ws.p);

        // My to local isn't working
        // ATM local is coincident with global, lets test that.
        Point3f p_local = ws.p;

        // Find widths?? Maybe x2
        // Really, I'd like a baked in answer, or a vertex sampling routine
        // When we get to patches it'll be even harder
        Float wid_x = norm(m_frame.s);
        Float wid_y = norm(m_frame.t);

        // The condition necessary is that the components of the wave remain
        // coherent throughout the whole extent of their travel.
        // Perhaps this means that on receive we need to take the wdf of
        // different sets of incoming rays and take the cross wdf.
        // Perhaps for each ray or packet, instead of summing, take the xwdf
        // then sum.
        // Showing this with single pulse temporal will be hard. Can show with
        // envelope detection, or continuous.
        // Two different experiments, beam steering and multipath.
        // Temporal phase effects in spatial ray-based radar rendering.

        // The wdf has been shown to be a useful tool in modelling spatial and
        // angular effects of scene elements in optical rendering, and easily
        // extends to beam effects which are prominent at radar wavelengths.
        // In this paper we show these effects for radar simulation. We extend
        // the framework, and show inclusion of temporal coherence/phase effects
        // which are also prominent in radar, in the form of phased-array beam
        // steering and multipath.

        // For beam steering, basically take the wdf of steering angle as a
        // function of x, and multiply wigners.

        // At each interaction (surface, emitter, receiver) apply a phase wigner
        // layer which is a phase shift proportional to travel length. This is
        // readily available and mostly unobtrusive as rays usually have path
        // lengths in this part of their calculations. ie part of the bsdf Now
        // encodes path length/phase. The counter is to just do this on
        // reception. Part of the reason we see anything when doing multipath
        // experiment with graham is that while the pulse is on, there is
        // continuous wave, so direct and indirect paths can have time together.
        // What that means is that the rays arriving in the same bin get
        // coherently summed, ie cross wignered/phased.
        // As part of paper, show 1d example in matlab/python as well as proper
        // implementation.

        // std::cout << ws.d.x() << std::endl;
        // std::cout << p_local.x() << std::endl;
        // std::cout << wid_x << std::endl;

        // Float wid_x =
        //     m_to_local.transform_affine(ScalarPoint3f(+1.f, -1.f, 0.f)) -
        //     m_to_local.transform_affine(ScalarPoint3f(-1.f, -1.f, 0.f));
        // Float wid_y =
        //     m_to_local.transform_affine(ScalarPoint3f(-1.f, +1.f, 0.f)) -
        //     m_to_local.transform_affine(ScalarPoint3f(-1.f, -1.f, 0.f));


        // Float g_static = 4.f*math::Pi<Float>*
        //     rcp(m_inv_surface_area)*(nu_0)*(nu_0);
        Float g_static = 4.f*math::Pi<Float>*
            rcp(ws.pdf)*(nu_0)*(nu_0);

        // Float g_angle = select(abs(p_local.x()/wid_x) < 0.5,
        //     1.0 - 2.0*abs(p_local.x()/wid_x), 0.f) *
        //     select(abs(p_local.y()/wid_y) < 0.5,
        //     1.0 - 2.0*abs(p_local.y()/wid_y), 0.f);

            // std::cout << g_angle << std::endl;
            // std::cout << abs(p_local.x()) << std::endl;
            // std::cout << abs(ws.p.x()) << std::endl;
        //
        Float g_angle = math::tri(p_local.x()/wid_x) *
            math::tri(p_local.y()/wid_y);
        // Float g_angle = trix*triy;

        // Float w_val = g_static;
        // Float w_val = g_angle;
        Float w_val = g_static * g_angle;

        ws.pdf = rcp(w_val);


        // jbo.K = @(z,nuxz, y,nuxy) ...
    	// 	jbo.D .* ...
    	// 	tri((z-jbo.l.z0)/jbo.l.zbw,2).*...
    	// 	sinc(2*nuxz.*tri((z-jbo.l.z0 )/jbo.l.zbw,2).*jbo.l.zbw) .*...
    	// 	tri((y-jbo.l.y0)/jbo.l.ybw,2).*...
    	// 	sinc(2*nuxy.*tri((y-jbo.l.y0 )/jbo.l.ybw,2).*jbo.l.ybw); % [W/W]

        // ds.pdf *= 4*math::Pi<Float>*nu_0^2;

        // Lets just take first wavelength element for now.

        // ds.pdf = rcp(wdf);

        // // Static gain of rectangular transmitter
        //
        // Float nu_0 = 1/si.wavelengths*1e9;
        //
        //
        //
        // Float g_static = 4*math::pi*m_surface_area*nu_0^2;
        //
        // // For single rectangular element we can use sample space to cut some
        // // corners. For more complex geometry, probably not.
        //
        // // convert ds.d to
        // // Float w_val = tri(m_to_local.transform_affine(ds.p))
        //     // * sinc(2*nu*tri(sample2.x - 0.5))
        //
        // // nuxz is the z component of the projection of the wavevector onto
        // // plane zy plane
        //
        // // this should be the dot product of u axis and wavevector
        //
        // Float w_val = tri(sample2.x - 0.5) * sinc(2*ds.d**tri(sample2.x - 0.5))

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
