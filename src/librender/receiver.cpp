#include <mitsuba/render/receiver.h>
#include <mitsuba/render/endpoint.h>

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/render/adc.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

// =============================================================================
// Sensor interface
// =============================================================================

MTS_VARIANT Receiver<Float, Spectrum>::Receiver(const Properties &props)
    : Base(props) {
    m_adc_sampling_start = props.float_("adc_sampling_start", 0.f);
    m_adc_sampling_time  = props.float_("adc_sampling_end", 0.f) - m_adc_sampling_start;

    if (m_adc_sampling_time < 0)
        Throw("ADC sampling time must be less than or equal to the adc "
              "sampling end time!");

    for (auto &[name, obj] : props.objects(false)) {
        auto *adc = dynamic_cast<ADC *>(obj.get());
        auto *sampler = dynamic_cast<Sampler *>(obj.get());

        if (adc) {
            if (m_adc)
                Throw("Only one adc can be specified per sensor.");
            m_adc = adc;
            props.mark_queried(name);
        } else if (sampler) {
            if (m_sampler)
                Throw("Only one sampler can be specified per sensor.");
            m_sampler = sampler;
            props.mark_queried(name);
        }
    }

    auto pmgr = PluginManager::instance();
    if (!m_adc) {
        // Instantiate an high dynamic range film if none was specified
        // Properties props_film("hdrfilm");
        Properties props_adc("hdradc");
        m_adc = static_cast<ADC *>(pmgr->create_object<ADC>(props_adc));
    }

    if (!m_sampler) {
        // Instantiate an independent filter with 4 samples/pixel if none was
        // specified
        Properties props_sampler("independent");
        props_sampler.set_int("sample_count", 4);
        m_sampler = static_cast<Sampler *>
            (pmgr->create_object<Sampler>(props_sampler));
    }

    m_resolution = ScalarVector2f(m_adc->window_size());
}

MTS_VARIANT Receiver<Float, Spectrum>::~Receiver() {}

// MTS_VARIANT typename Receiver<Float, Spectrum>::Wavelength
// Receiver<Float, Spectrum>::doppler(SurfaceInteraction3f si, Mask active) const {
//     // std::cout << m_velocity.get() << std::endl;
//     // std::cout << m_world_transform.get() << std::endl;
//     // return select(active, 2*dot(si.wi, m_velocity.get()*Point3f(si.to_local(si.p))) / math::CVac<float>
//     //  * si.wavelengths, 0.f);
//     Wavelength result = 0.f;
//     // std::cout << any_or<true>(neq(m_shape, nullptr)) << std::endl;
//     if(any_or<true>(neq(m_shape, nullptr))){
//         // std::cout << m_shape << std::endl;
//         result = select(active, 2*m_shape->doppler(si, active), 0.f);
//     } else {
//         result = select(active, 2*dot(si.wi, m_velocity.get()->eval(0.f)*Point3f(si.to_local(si.p))) / math::CVac<float>
//          * si.wavelengths, 0.f);
//     }
//     return result;
// }

MTS_VARIANT std::pair<
    typename Receiver<Float, Spectrum>::RayDifferential3f, Spectrum>
Receiver<Float, Spectrum>::sample_ray_differential(Float time, Float sample1,
                                                    const Point2f &sample2,
                                                    const Point2f &sample3,
                                                    Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

    // Note, this is for testing how the output changes wrt small differences
    // in position. Changing the position on the sensor is still useful, but
    // care must be taken in temporal cases. The resolution it NOW refers to
    // is time,range/frequency,doppler bins.

    auto [temp_ray, result_spec] =
        sample_ray(time, sample1, sample2, sample3, active);

    RayDifferential result_ray(temp_ray);

    Vector2f dx(1.f / m_resolution.x(), 0.f);
    Vector2f dy(0.f, 1.f / m_resolution.y());

    // Sample a result_ray for X+1
    std::tie(temp_ray, std::ignore) =
        sample_ray(time, sample1, sample2 + dx, sample3, active);

    result_ray.o_x = temp_ray.o;
    result_ray.d_x = temp_ray.d;

    // Sample a result_ray for Y+1
    std::tie(temp_ray, std::ignore) =
        sample_ray(time, sample1, sample2 + dy, sample3, active);

    result_ray.o_y = temp_ray.o;
    result_ray.d_y = temp_ray.d;
    result_ray.has_differentials = true;

    return { result_ray, result_spec };
}

// =============================================================================
// ProjectiveCamera interface
// =============================================================================
//
// MTS_VARIANT ProjectiveCamera<Float, Spectrum>::ProjectiveCamera(const Properties &props)
//     : Base(props) {
//     /* Distance to the near clipping plane */
//     m_near_clip = props.float_("near_clip", 1e-2f);
//     /* Distance to the far clipping plane */
//     m_far_clip = props.float_("far_clip", 1e4f);
//     /* Distance to the focal plane */
//     m_focus_distance = props.float_("focus_distance", m_far_clip);
//
//     if (m_near_clip <= 0.f)
//         Throw("The 'near_clip' parameter must be greater than zero!");
//     if (m_near_clip >= m_far_clip)
//         Throw("The 'near_clip' parameter must be smaller than 'far_clip'.");
//
// }
//
// MTS_VARIANT ProjectiveCamera<Float, Spectrum>::~ProjectiveCamera() { }
//
// // =============================================================================
// // Helper functions
// // =============================================================================
//
// float parse_fov(const Properties &props, float aspect) {
//     if (props.has_property("fov") && props.has_property("focal_length"))
//         Throw("Please specify either a focal length ('focal_length') or a "
//                 "field of view ('fov')!");
//
//     float fov;
//     std::string fov_axis;
//
//     if (props.has_property("fov")) {
//         fov = props.float_("fov");
//
//         fov_axis = string::to_lower(props.string("fov_axis", "x"));
//
//         if (fov_axis == "smaller")
//             fov_axis = aspect > 1 ? "y" : "x";
//         else if (fov_axis == "larger")
//             fov_axis = aspect > 1 ? "x" : "y";
//     } else {
//         std::string f = props.string("focal_length", "50mm");
//         if (string::ends_with(f, "mm"))
//             f = f.substr(0, f.length()-2);
//
//         float value;
//         try {
//             value = std::stof(f);
//         } catch (...) {
//             Throw("Could not parse the focal length (must be of the form "
//                 "<x>mm, where <x> is a positive integer)!");
//         }
//
//         fov = 2.f *
//                 rad_to_deg(std::atan(std::sqrt(float(36 * 36 + 24 * 24)) / (2.f * value)));
//         fov_axis = "diagonal";
//     }
//
//     float result;
//     if (fov_axis == "x") {
//         result = fov;
//     } else if (fov_axis == "y") {
//         result = rad_to_deg(
//             2.f * std::atan(std::tan(.5f * deg_to_rad(fov)) * aspect));
//     } else if (fov_axis == "diagonal") {
//         float diagonal = 2.f * std::tan(.5f * deg_to_rad(fov));
//         float width    = diagonal / std::sqrt(1.f + 1.f / (aspect * aspect));
//         result = rad_to_deg(2.f * std::atan(width*.5f));
//     } else {
//         Throw("The 'fov_axis' parameter must be set to one of 'smaller', "
//                 "'larger', 'diagonal', 'x', or 'y'!");
//     }
//
//     if (result <= 0.f || result >= 180.f)
//         Throw("The horizontal field of view must be in the range [0, 180]!");
//
//     return result;
// }

MTS_IMPLEMENT_CLASS_VARIANT(Receiver, Endpoint, "receiver")
// MTS_IMPLEMENT_CLASS_VARIANT(ProjectiveCamera, Receiver)

MTS_INSTANTIATE_CLASS(Receiver)
// MTS_INSTANTIATE_CLASS(ProjectiveCamera)
NAMESPACE_END(mitsuba)
