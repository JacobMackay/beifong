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

Point transmitter: Source comes from a point as in the spotlight case, but has
an attached shape which provides the beamprofile/fourier transform.
This differs from the wigner transform version which a joint position/direction
distribution, whereas this fourier transform collapses the position into
direction.

 */

template <typename Float, typename Spectrum>
class PointTransmitter final : public Emitter<Float, Spectrum> {
 public:
     MTS_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
     MTS_IMPORT_TYPES(Scene, Shape, Texture)

     // Constructor
     PointTransmitter(const Properties &props) : Base(props) {
         if (props.has_property("to_world")) {
             Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");
              }

         m_flags = +EmitterFlags::DeltaPosition;
         m_intensity = props.texture<Texture>("intensity", Texture::D65(1.f));
         m_texture = props.texture<Texture>("texture", Texture::D65(1.f));

         if (m_intensity->is_spatially_varying())
             Throw("The parameter 'intensity' cannot be spatially"
                    "varying (e.g. bitmap type)!");

         if (props.has_property("texture")) {
             if (!m_texture->is_spatially_varying())
                 Throw("The parameter 'texture' must be spatially"
                        "varying (e.g. bitmap type)!");
             m_flags |= +EmitterFlags::SpatiallyVarying;
         }

         // Atm this is symetric beam and kept just to please routine
         // Either render cutoff or null2null beamwidth
         m_cutoff_angle = props.float_("cutoff_angle", 180.0f);
         // This is -3dB beamwidth
         m_beam_width = props.float_("beam_width",
                m_cutoff_angle * 3.0f / 4.0f);
         m_cutoff_angle = deg_to_rad(m_cutoff_angle);
         m_beam_width = deg_to_rad(m_beam_width);
         m_inv_transition_width = 1.0f / (m_cutoff_angle - m_beam_width);
         m_cos_cutoff_angle = cos(m_cutoff_angle);
         m_cos_beam_width = cos(m_beam_width);
         Assert(m_cutoff_angle >= m_beam_width);
         m_uv_factor = tan(m_cutoff_angle);

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
         // m_radiance = props.texture<Texture>("radiance", Texture::D65(1.f));

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
         if (m_intensity->is_spatially_varying()) {
             m_flags |= +EmitterFlags::SpatiallyVarying;
         }
     }

     UnpolarizedSpectrum falloff_curve(const Vector3f &d,
        Wavelength wavelengths, Mask active) const {
         SurfaceInteraction3f si = zero<SurfaceInteraction3f>();
         si.wavelengths = wavelengths;
         UnpolarizedSpectrum result = m_intensity->eval(si, active);

         Vector3f local_dir = normalize(d);
         Float cos_theta = local_dir.z();

         if (m_texture->is_spatially_varying()) {
             si.uv = Point2f(.5f + .5f * local_dir.x() /
                                (local_dir.z() * m_uv_factor),
                             .5f + .5f * local_dir.y() /
                                (local_dir.z() * m_uv_factor));
             result *= m_texture->eval(si, active);
         }

         UnpolarizedSpectrum beam_res = select(cos_theta >= m_cos_beam_width,
             result, result * m_shape->fourier_weight(local_dir, wavelengths));

         return select(cos_theta <= m_cos_cutoff_angle,
             UnpolarizedSpectrum(0.0f), beam_res);
     }

     std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                           const Point2f &spatial_sample,
                                           const Point2f & /*dir_sample*/,
                                           Mask active) const override {
         MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

         // 1. Sample directional component
         // This is still uniformish
         const Transform4f &trafo = m_shape->to_world();
         Vector3f local_dir = warp::square_to_cosine_hemisphere(spatial_sample);
         Vector3f pdf_dir =
            warp::square_to_cosine_hemisphere_pdf(local_dir);

         // 2. Sample spectrum
         auto [wavelengths, spec_weight] = m_intensity->sample_spectrum(
             zero<SurfaceInteraction3f>(),
             math::sample_shifted<Wavelength>(wavelength_sample), active);

         UnpolarizedSpectrum falloff_spec =
            falloff_curve(local_dir, wavelengths, active);

         return {
             Ray3f(trafo.translation(), trafo * local_dir, time, wavelengths),
                 unpolarized<Spectrum>(falloff_spec / pdf_dir) };
     }

     std::pair<DirectionSample3f, Spectrum>
        sample_direction(const Interaction3f &it,
                         const Point2f &/*sample*/,
                         Mask active) const override {
         MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);

         const Transform4f &trafo = m_shape->to_world();

         DirectionSample3f ds;
         ds.p        = trafo.translation();
         ds.n        = 0.f;
         ds.uv       = 0.f;
         ds.pdf      = 1.0f;
         ds.time     = it.time;
         ds.delta    = true;
         ds.object   = this;
         ds.d        = ds.p - it.p;
         ds.dist     = norm(ds.d);
         Float inv_dist = rcp(ds.dist);
         ds.d        *= inv_dist;
         Vector3f local_d = trafo.inverse() * -ds.d;
         UnpolarizedSpectrum falloff_spec =
            falloff_curve(local_d, it.wavelengths, active);

         return {
             ds, unpolarized<Spectrum>(falloff_spec * (inv_dist * inv_dist)) };
     }

     Float pdf_direction(const Interaction3f &,
                         const DirectionSample3f &, Mask) const override {
         return 0.f;
     }

     Spectrum eval(const SurfaceInteraction3f &, Mask) const override {
         return 0.f;
     }

     ScalarBoundingBox3f bbox() const override {
         // return m_world_transform->translation_bounds();
         return m_shape->bbox();
     }

     void traverse(TraversalCallback *callback) override {
         callback->put_object("intensity", m_intensity.get());
         callback->put_object("texture", m_texture.get());
     }

     std::string to_string() const override {
         std::ostringstream oss;
         oss << "PointTransmitter[" << std::endl
             << "  world_transform = "
                << string::indent(m_shape->to_world()) << "," << std::endl
             << "  intensity = " << m_intensity << "," << std::endl
             << "  cutoff_angle = " << m_cutoff_angle << "," << std::endl
             << "  beam_width = " << m_beam_width << "," << std::endl
             << "  texture = " << (m_texture ? string::indent(m_texture) : "")
             << "  medium = " << (m_medium ? string::indent(m_medium) : "")
             << "]";
         return oss.str();
     }

    // std::string to_string() const override {
    //     std::ostringstream oss;
    //     oss << "AreaLight[" << std::endl
    //         << "  radiance = " << string::indent(m_radiance)
    //      << "," << std::endl
    //         << "  surface_area = ";
    //     if (m_shape) oss << m_shape->surface_area();
    //     else         oss << "  <no shape attached!>";
    //     oss << "," << std::endl;
    //     if (m_medium) oss << string::indent(m_medium);
    //     else         oss << "  <no medium attached!>";
    //     oss << std::endl << "]";
    //     return oss.str();
    // }

    MTS_DECLARE_CLASS()

 private:
    // ref<Texture> m_radiance;
    // spectrum m_power;
    // ref<Texture> m_radiance;

    ref<Texture> m_intensity;
    ref<Texture> m_texture;
    ScalarFloat m_beam_width, m_cutoff_angle, m_uv_factor;
    ScalarFloat m_cos_beam_width, m_cos_cutoff_angle, m_inv_transition_width;
};

MTS_IMPLEMENT_CLASS_VARIANT(PointTransmitter, Emitter)
MTS_EXPORT_PLUGIN(PointTransmitter, "Point Transmitter")
NAMESPACE_END(mitsuba)
