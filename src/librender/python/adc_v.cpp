#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/filesystem.h>
#include <mitsuba/render/adc.h>
#include <mitsuba/render/signalblock.h>
#include <mitsuba/core/rfilter.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/spiral.h>
#include <mitsuba/python/python.h>

MTS_PY_EXPORT(ADC) {
    MTS_PY_IMPORT_TYPES(ADC)
    MTS_PY_CLASS(ADC, Object)
        .def_method(ADC, prepare, "channels"_a)
        .def_method(ADC, put, "block"_a)
        .def_method(ADC, set_destination_file, "filename"_a)
        .def("develop", py::overload_cast<>(&ADC::develop))
        .def("develop", py::overload_cast<const ScalarPoint2i &, const ScalarVector2i &,
                                            const ScalarPoint2i &, Bitmap *>(
                &ADC::develop, py::const_),
            "offset"_a, "size"_a, "target_offset"_a, "target"_a)
        .def_method(ADC, destination_exists, "basename"_a)
        .def_method(ADC, bitmap, "raw"_a = false)
        .def_method(ADC, has_high_quality_edges)
        .def_method(ADC, size)
        .def_method(ADC, window_size)
        .def_method(ADC, window_offset)
        // .def_method(ADC, bandwidth)
        // .def_method(ADC, centres)
        .def_method(ADC, set_window)
        .def_method(ADC, reconstruction_filter);
}
