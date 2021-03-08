#include <mitsuba/render/transmitter.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/python/python.h>

/// Trampoline for derived types implemented in Python
MTS_VARIANT class PyTransmitter : public Transmitter<Float, Spectrum> {
public:
    MTS_IMPORT_TYPES(Transmitter)

    PyTransmitter(const Properties &props) : Transmitter(props) { }

    std::pair<Ray3f, Spectrum>
    sample_ray(Float time, Float sample1, const Point2f &sample2,
           const Point2f &sample3, Mask active) const override {
        using Return = std::pair<Ray3f, Spectrum>;
        PYBIND11_OVERLOAD_PURE(Return, Transmitter, sample_ray, time, sample1, sample2, sample3,
                               active);
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &ref,
                     const Point2f &sample,
                     Mask active) const override {
        using Return = std::pair<DirectionSample3f, Spectrum>;
        PYBIND11_OVERLOAD_PURE(Return, Transmitter, sample_direction, ref, sample, active);
    }

    Float pdf_direction(const Interaction3f &ref,
                        const DirectionSample3f &ds,
                        Mask active) const override {
        PYBIND11_OVERLOAD_PURE(Float, Transmitter, pdf_direction, ref, ds, active);
    }

    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        PYBIND11_OVERLOAD_PURE(Spectrum, Transmitter, eval, si, active);
    }

    Float eval_signal(Float time, Float frequency) const {
        PYBIND11_OVERLOAD_PURE(Float, Transmitter, eval_signal, time, frequency);
    }

    ScalarBoundingBox3f bbox() const override {
        PYBIND11_OVERLOAD_PURE(ScalarBoundingBox3f, Transmitter, bbox,);
    }


    std::string to_string() const override {
        PYBIND11_OVERLOAD_PURE(std::string, Transmitter, to_string,);
    }
};

MTS_PY_EXPORT(Transmitter) {
    MTS_PY_IMPORT_TYPES()
    using PyTransmitter = PyTransmitter<Float, Spectrum>;

    auto transmitter = py::class_<Transmitter, PyTransmitter, Endpoint, ref<Transmitter>>(m, "Transmitter", D(Transmitter))
        .def(py::init<const Properties&>())
        .def_method(Transmitter, is_environment)
        .def_method(Transmitter, flags);

    if constexpr (is_cuda_array_v<Float>)
        pybind11_type_alias<UInt64, TransmitterPtr>();

    if constexpr (is_array_v<Float>) {
        transmitter.def_static("sample_ray_vec",
                            vectorize([](const TransmitterPtr &ptr, Float time, Float sample1,
                                         const Point2f &sample2, const Point2f &sample3,
                                         Mask active) {
                                 return ptr->sample_ray(time, sample1, sample2, sample3, active);
                            }),
                            "ptr"_a, "time"_a, "sample1"_a, "sample2"_a, "sample3"_a,
                            "active"_a = true,
                            D(Endpoint, sample_ray));
        transmitter.def_static("sample_direction_vec",
                            vectorize([](const TransmitterPtr &ptr, const Interaction3f &it,
                                         const Point2f &sample, Mask active) {
                                return ptr->sample_direction(it, sample, active);
                            }),
                            "ptr"_a, "it"_a, "sample"_a, "active"_a = true,
                            D(Endpoint, sample_direction));
        transmitter.def_static("pdf_direction_vec",
                            vectorize([](const TransmitterPtr &ptr, const Interaction3f &it,
                                         const DirectionSample3f &ds, Mask active) {
                                return ptr->pdf_direction(it, ds, active);
                            }),
                            "ptr"_a, "it"_a, "ds"_a, "active"_a = true,
                            D(Endpoint, pdf_direction));
        transmitter.def_static("eval_vec",
                           vectorize([](const TransmitterPtr &ptr, const SurfaceInteraction3f &si,
                                        Mask active) { return ptr->eval(si, active); }),
                           "ptr"_a, "si"_a, "active"_a = true, D(Endpoint, eval));
       transmitter.def_static("eval_signal_vec",
                          vectorize([](const TransmitterPtr &ptr, Float time,
                                       Float frequency) { return ptr->eval_signal(time, frequency); }),
                          "ptr"_a, "time"_a, "frequency"_a, D(Endpoint, eval));
    }

    MTS_PY_REGISTER_OBJECT("register_transmitter", Transmitter)
}
