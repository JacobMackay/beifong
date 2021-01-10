#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/signalblock.h>
#include <mitsuba/python/python.h>

MTS_PY_EXPORT(SignalBlock) {
    MTS_PY_IMPORT_TYPES(SignalBlock, ReconstructionFilter)
    MTS_PY_CLASS(SignalBlock, Object)
        .def(py::init<const ScalarVector2i &, size_t,
                const ReconstructionFilter *, bool, bool, bool, bool>(),
            "size"_a, "channel_count"_a, "filter"_a = nullptr,
            "warn_negative"_a = true, "warn_invalid"_a = true,
            "border"_a = true, "normalize"_a = false)
        .def("put", py::overload_cast<const SignalBlock *>(&SignalBlock::put),
            D(SignalBlock, put), "block"_a)
        .def("put", vectorize(py::overload_cast<const Point2f &,
            const wavelength_t<Spectrum> &, const Spectrum &, const Float &,
            mask_t<Float>>(&SignalBlock::put)),
            "pos"_a, "wavelengths"_a, "value"_a, "alpha"_a = 1.f, "active"_a = true,
            D(SignalBlock, put, 2))
        .def("put",
            [](SignalBlock &ib, const Point2f &pos,
                const std::vector<Float> &data, Mask mask) {
                if (data.size() != ib.channel_count())
                    throw std::runtime_error("Incompatible channel count!");
                ib.put(pos, data.data(), mask);
            }, "pos"_a, "data"_a, "active"_a = true)
        .def_method(SignalBlock, clear)
        .def_method(SignalBlock, set_offset, "offset"_a)
        .def_method(SignalBlock, offset)
        .def_method(SignalBlock, size)
        .def_method(SignalBlock, width)
        .def_method(SignalBlock, height)
        .def_method(SignalBlock, warn_invalid)
        .def_method(SignalBlock, warn_negative)
        .def_method(SignalBlock, set_warn_invalid, "value"_a)
        .def_method(SignalBlock, set_warn_negative, "value"_a)
        .def_method(SignalBlock, border_size)
        .def_method(SignalBlock, channel_count)
        .def("data", py::overload_cast<>(&SignalBlock::data, py::const_), D(SignalBlock, data));
}
