#include <mitsuba/render/transmitter.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/python/python.h>

MTS_PY_EXPORT(TransmitterExtras) {
    py::enum_<TransmitterFlags>(m, "TransmitterFlags", D(TransmitterFlags))
        .def_value(TransmitterFlags, None)
        .def_value(TransmitterFlags, DeltaPosition)
        .def_value(TransmitterFlags, DeltaDirection)
        .def_value(TransmitterFlags, Infinite)
        .def_value(TransmitterFlags, Surface)
        .def_value(TransmitterFlags, SpatiallyVarying)
        .def_value(TransmitterFlags, Delta)
        .def(py::self == py::self)
        .def(py::self | py::self)
        .def(int() | py::self)
        .def(py::self & py::self)
        .def(int() & py::self)
        .def(+py::self)
        .def(~py::self)
        .def("__pos__", [](const TransmitterFlags &f) {
            return static_cast<uint32_t>(f);
        }, py::is_operator());
}
