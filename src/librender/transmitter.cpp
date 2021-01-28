#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/transmitter.h>
#include <mitsuba/render/endpoint.h>

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT Transmitter<Float, Spectrum>::Transmitter(const Properties &props) : Base(props) { }
MTS_VARIANT Transmitter<Float, Spectrum>::~Transmitter() { }

MTS_IMPLEMENT_CLASS_VARIANT(Transmitter, Endpoint, "transmitter")
MTS_INSTANTIATE_CLASS(Transmitter)
NAMESPACE_END(mitsuba)
