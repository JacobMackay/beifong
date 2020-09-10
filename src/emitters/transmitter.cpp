/* For copyright see the LICENSE.txt file by mitsuba 2. Changes made here by
Jacob Mackay should be free */

#include "../../include/mitsuba/core/properties.h"
#include "../../include/mitsuba/core/warp.h"
#include "../../include/mitsuba/core/spectrum.h"
#include "../../include/mitsuba/render/emitter.h"
#include "../../include/mitsuba/render/medium.h"
#include "../../include/mitsuba/render/shape.h"
#include "../../include/mitsuba/render/texture.h"

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
class Transmitter final : public Emitter<Float, Spectrum> {
 public:
     MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
     MTS_IMPORT_TYPES(Scene, Shape, Texture)

     // Constructor
     Transmitter(const Properties &props) : Base(props) {
         if (props.has_property("to_world")) {
             Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");
              }

         // a 4d bitmap is useful for wigner, textures & x,y,range,doppler
         // if i can provide a bitmap texture in 4d, thats great!
         // textures are usually in uv coordinates
         // this texture algorithm could be a contribution in itself.

         // m_ele_fov = deg_to_rad(props.float_("fov_elevation", 90.0f));
         // m_azi_fov = deg_to_rad(props.float_("fov_azimuth", 90.0f));

         // Default to 94 GHz antenna centre frequency
         // m_centre_frequency = props.float_("antenna_centre_freq", 94E9);

         // TODO(Jacob): Automatically calculate SWR bandwidth based on...
         // VSWR bandwidth, resonance bandwidth, Q factor, many things
         // For now just guess for 94GHz + 6GHz classic mmw radar.
         // m_antenna_bandwidth = props.float_("antenna_bandwidth", 10E9);

         // Can be a float, a spectrum, or a texture. A texture would be like
         // the non-uniform illumination or a phase weight.
         // TODO(Jacob): replace D65 with either uniform or a tapered bandwidth
         // m_power = props.texture<Texture>("power", Texture::D65(1.f));

         // Don't know if this is correct, but we'll see, W/sr/m^2
         // This is a function in wavelength and u,v which is multiplied by
         // an incoming ray
         // TODO(Jacob): Confirm if correct
         // m_radiance = m_power/(4*math::Pi<Float>*m_shape->surface_area());
         m_radiance = props.texture<Texture>("radiance", Texture::D65(1.f));

         // should there be class antenna? That has tx and rx? or is that more
         // difficult?
         // emitter should be able to get a child object which is signal
         // m_tsignal: u(t), as list of 2 elelents
         // m_fsignal: u(f), as list of 2 elements
         // m_tfsignal: u(t,f), as list of 3 elements
         // m_power
         // m_frequency
         // m_units: 'voltage', 'power'
         // for now, just accept m_power and (m_frequency)-->(m_wavelength)

         // either specify a float for power (uniform), a given spectra,
         // or centre antenna resonance frequency and bandwidth
         // m_radiance = props.texture<Texture>("power",
         // Texture::resonance_bandwidth(1.f, f_centre, f_bw));

         // radiance(u,v,lambda)
         // What would be better here would be to specify the antenna bandwidth

         m_flags = +EmitterFlags::Surface;
         if (m_radiance->is_spatially_varying()) {
             m_flags |= +EmitterFlags::SpatiallyVarying;
         }
     }

     Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
         MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

         // Return the radiance along a ray si at locaton si.p and direction
         // si.wi (incident direction in local shading frame)
         // This asks for 'some' evaluation of a surface interaction. Return
         // m_radiance at si. Unless m_radiance has an angular component,
         // this will still just be uniform.
         return select(
             Frame3f::cos_theta(si.wi) > 0.f,
             unpolarized<Spectrum>(m_radiance->eval(si, active)),
             0.f);
         }

     std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                              const Point2f &sample2, const Point2f &sample3,
                              Mask active) const override {
         MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

         SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
         si.t = math::Infinity<Float>;

         Float pdf = 1.f;
         // Float spatial_weight = 1.f;

         // Proper WDF sampling samples a 6d function, u,v, θ,φ, λ,t
         // For now sample these separately.

         // 1. Two strategies to sample spatial component based on 'm_radiance'
         if (!m_radiance->is_spatially_varying()) {
             // Sample a position on the shape surface.
             PositionSample3f ps =
                m_shape->sample_position(time, sample2, active);

             // Radiance not spatially varying, use area-based sampling of shape
             si = SurfaceInteraction3f(ps, zero<Wavelength>());
             pdf = ps.pdf;
             // spatial_weight = select(neq(pdf, 0.f), 1.0 / ps.pdf, 0.f);

             // for tomorrow:
             // extend shape class to have sample_fourier and sample_wigner
             // do it for rectangle: we have closed form solutions
             // https://math.stackexchange.com/questions/440202/definition-of-the-fourier-transform-of-function-on-the-sphere
             // https://mathoverflow.net/questions/149692/fourier-transform-of-the-unit-sphere

             // (shouldn't need position sample)
             // m_shape->sample_fourier_direction(time,
             // sample3, wavelength, active) (this is based on ft)
             // don't need joint in traditional case cause they're separable
             // m_shape->sample_qjoint(time, sample2, sample3,
             //     wavelength, active)

             // I can simply make a shape called phased array
             // it has a base shape and an array of vectors pointing to element
             // centres
             // OR it has spacing (centre to centre) and number of elements in
             // x and y
             // also takes a phase offset. either just one, or one for each
             // element
             // can also take steering angle (theta, phi)

             // there would be some optimisation for rendering disjointed
             // pixels, but I'm not doing that.
             // just rendering each pixel separately.
             // but we are sampled over [0,1], maybe can still map that to
             // shapes

             // paper can show that it doesn't REALLY work if taken just in
             // normal space and needs wdf

             // a potential solution to diffraction around edges is at scene
             // compilation time, add a medium or material between things and
             // let the rays interact with that
             // or add object.ether box, like bounding region. if ray.intersect
             // doesn't hit, the reevaluate with the ether and count that as a
             // hit. what about parallel rays??

         } else {
             std::tie(si.uv, pdf) =
                m_radiance->sample_position(sample2, active);
             active &= neq(pdf, 0.f);

             si = m_shape->eval_parameterization(Point2f(si.uv), active);
             active &= si.is_valid();

             pdf /= norm(cross(si.dp_du, si.dp_dv));
             // spatial_weight *= norm(cross(si.dp_du, si.dp_dv));
         }

         // So I'd need to create a pos/direc 4d pdf
         // Wait, I guess pdf is just the probability @ sample2,sample3

         // For now ignore shapes like spheres and whip antennas which have
         // convex emission profiles.
         // 2. Sample directional component
         Vector3f local = warp::square_to_cosine_hemisphere(sample3);

         Wavelength wavelength;
         Spectrum spec_weight;
         Spectrum angular_weight = math::Pi<Float>;

         if constexpr (is_spectral_v<Spectrum>) {
             std::tie(wavelength, spec_weight) = m_radiance->sample_spectrum(
                 si, math::sample_shifted<Wavelength>(wavelength_sample),
                 active);
                 angular_weight *= m_shape->fourier_weight(local, wavelength);
         } else {
             wavelength = zero<Wavelength>();
             spec_weight = m_radiance->eval(si, active);
         }

         return std::make_pair(
             Ray3f(si.p, si.to_world(local), time, wavelength),
             unpolarized<Spectrum>(spec_weight) *
             unpolarized<Spectrum>(angular_weight) / pdf);
         // return std::make_pair(
         //     Ray3f(si.p, si.to_world(local), time, wavelength),
         //     spatial_weight * unpolarized<Spectrum>(spec_weight) *
         //     unpolarized<Spectrum>(angular_weight));
     }

     // Direction sample is special/different. can interpret sample as position
     // or direction on emitter. For not I will interpret as position, as it is
     // also a position.
    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it,
                     const Point2f &sample,
                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape,
            "Can't sample from an area emitter without an associated Shape.");
        DirectionSample3f ds;
        Spectrum spec;
        Spectrum angular_weight;

        // One of two very different strategies is used depending on
        // 'm_radiance'
        if (!m_radiance->is_spatially_varying()) {
            // Texture is uniform, try to importance sample the shape wrt.
            // solid angle at 'it'
            ds = m_shape->sample_direction(it, sample, active);
            active &= dot(ds.d, ds.n) < 0.f && neq(ds.pdf, 0.f);

            SurfaceInteraction3f si(ds, it.wavelengths);
            spec = m_radiance->eval(si, active) / ds.pdf;

            angular_weight = m_shape->fourier_weight(si.wi, it.wavelengths);

        } else {
            // Importance sample the texture, then map onto the shape
            auto [uv, pdf] = m_radiance->sample_position(sample, active);
            active &= neq(pdf, 0.f);

            SurfaceInteraction3f si =
                m_shape->eval_parameterization(uv, active);
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
            // the pdf is the pdf at a point in space, not direction
            ds.pdf = select(active, pdf / norm(cross(si.dp_du, si.dp_dv)) *
                                        dist_squared / -dp, 0.f);

            spec = m_radiance->eval(si, active) / ds.pdf;

            angular_weight = m_shape->fourier_weight(si.wi, it.wavelengths);
        }

        // TODO(Jacob): probably angular_weight[0]
        ds.pdf /= angular_weight[0];

        ds.object = this;
        return { ds, unpolarized<Spectrum>(spec) *
            unpolarized<Spectrum>(angular_weight) & active };
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);
        Float dp = dot(ds.d, ds.n);
        active &= dp < 0.f;

        Float value;
        Spectrum angular_weight;
        if (!m_radiance->is_spatially_varying()) {
            value = m_shape->pdf_direction(it, ds, active);
            SurfaceInteraction3f si(ds, it.wavelengths);
            angular_weight = m_shape->fourier_weight(si.wi, it.wavelengths);
        } else {
            // This surface intersection would be nice to avoid..
            SurfaceInteraction3f si =
                m_shape->eval_parameterization(ds.uv, active);
            active &= si.is_valid();

            value = m_radiance->pdf_position(ds.uv, active) * sqr(ds.dist) /
                    (norm(cross(si.dp_du, si.dp_dv)) * -dp);

            angular_weight = m_shape->fourier_weight(si.wi, it.wavelengths);
        }

        return select(active, value/angular_weight[0], 0.f);
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
    // ref<Texture> m_radiance;
    // spectrum m_power;
    ref<Texture> m_radiance;
};

MTS_IMPLEMENT_CLASS_VARIANT(Transmitter, Emitter)
MTS_EXPORT_PLUGIN(Transmitter, "Transmitter")
NAMESPACE_END(mitsuba)
