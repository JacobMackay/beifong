#include <mitsuba/render/receiver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/python/python.h>

/// Trampoline for derived types implemented in Python
MTS_VARIANT class PyReceiver : public Receiver<Float, Spectrum> {
public:
    MTS_IMPORT_TYPES(Receiver)

    PyReceiver(const Properties &props) : Receiver(props) { }

    std::pair<Ray3f, Spectrum>
    sample_ray(Float time, Float sample1, const Point2f &sample2,
           const Point2f &sample3, Mask active) const override {
        using Return = std::pair<Ray3f, Spectrum>;
        PYBIND11_OVERLOAD_PURE(Return, Receiver, sample_ray, time, sample1, sample2, sample3,
                               active);
    }

    std::pair<RayDifferential3f, Spectrum>
    sample_ray_differential(Float time, Float sample1, const Point2f &sample2,
                            const Point2f &sample3, Mask active) const override {
        using Return = std::pair<RayDifferential3f, Spectrum>;
        PYBIND11_OVERLOAD_PURE(Return, Receiver, sample_ray_differential, time, sample1, sample2, sample3,
                               active);
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &ref,
                     const Point2f &sample,
                     Mask active) const override {
        using Return = std::pair<DirectionSample3f, Spectrum>;
        PYBIND11_OVERLOAD_PURE(Return, Receiver, sample_direction, ref, sample, active);
    }

    Float pdf_direction(const Interaction3f &ref,
                        const DirectionSample3f &ds,
                        Mask active) const override {
        PYBIND11_OVERLOAD_PURE(Float, Receiver, pdf_direction, ref, ds, active);
    }

    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        PYBIND11_OVERLOAD_PURE(Spectrum, Receiver, eval, si, active);
    }

    ScalarBoundingBox3f bbox() const override {
        PYBIND11_OVERLOAD_PURE(ScalarBoundingBox3f, Receiver, bbox,);
    }

    std::string to_string() const override {
        PYBIND11_OVERLOAD_PURE(std::string, Receiver, to_string,);
    }
};

MTS_PY_EXPORT(Receiver) {
    MTS_PY_IMPORT_TYPES(Receiver, Endpoint)
    using PyReceiver = PyReceiver<Float, Spectrum>;

    py::class_<Receiver, PyReceiver, Endpoint, ref<Receiver>>(m, "Receiver", D(Receiver))
        .def(py::init<const Properties&>())
        .def("sample_ray_differential", vectorize(&Receiver::sample_ray_differential),
            "time"_a, "sample1"_a, "sample2"_a, "sample3"_a, "active"_a = true)
        .def_method(Receiver, adc_sampling_start)
        .def_method(Receiver, adc_sampling_time)
        .def_method(Receiver, needs_aperture_sample)
        .def("adc", py::overload_cast<>(&Receiver::adc, py::const_), D(Receiver, adc))
        .def("sampler", py::overload_cast<>(&Receiver::sampler, py::const_), D(Receiver, sampler));

    MTS_PY_REGISTER_OBJECT("register_receiver", Receiver)
}
