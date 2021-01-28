#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/endpoint.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)


/**
 * \brief This list of flags is used to classify the different types of transmitters.
 */
enum class TransmitterFlags : uint32_t {
    // =============================================================
    //                      Transmitter types
    // =============================================================

    /// No flags set (default value)
    None                 = 0x00000,

    /// The transmitter lies at a single point in space
    DeltaPosition        = 0x00001,

    /// The transmitter emits light in a single direction
    DeltaDirection       = 0x00002,

    /// The transmitter is placed at infinity (e.g. environment maps)
    Infinite             = 0x00004,

    /// The transmitter is attached to a surface (e.g. area transmitters)
    Surface              = 0x00008,

    // =============================================================
    //!                   Other lobe attributes
    // =============================================================

    /// The emission depends on the UV coordinates
    SpatiallyVarying     = 0x00010,

    // =============================================================
    //!                 Compound lobe attributes
    // =============================================================

    /// Delta function in either position or direction
    Delta        = DeltaPosition | DeltaDirection,
};

constexpr uint32_t operator |(TransmitterFlags f1, TransmitterFlags f2)  { return (uint32_t) f1 | (uint32_t) f2; }
constexpr uint32_t operator |(uint32_t f1, TransmitterFlags f2)      { return f1 | (uint32_t) f2; }
constexpr uint32_t operator &(TransmitterFlags f1, TransmitterFlags f2)  { return (uint32_t) f1 & (uint32_t) f2; }
constexpr uint32_t operator &(uint32_t f1, TransmitterFlags f2)      { return f1 & (uint32_t) f2; }
constexpr uint32_t operator ~(TransmitterFlags f1)                   { return ~(uint32_t) f1; }
constexpr uint32_t operator +(TransmitterFlags e)                    { return (uint32_t) e; }
template <typename UInt32>
constexpr auto has_flag(UInt32 flags, TransmitterFlags f)            { return neq(flags & (uint32_t) f, 0u); }



template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER Transmitter : public Endpoint<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Endpoint)

    /// Is this an environment map light transmitter?
    bool is_environment() const {
        return has_flag(m_flags, TransmitterFlags::Infinite) && !has_flag(m_flags, TransmitterFlags::Delta);
    }

    /// Flags for all components combined.
    uint32_t flags(mask_t<Float> /*active*/ = true) const { return m_flags; }


    ENOKI_CALL_SUPPORT_FRIEND()
    MTS_DECLARE_CLASS()
protected:
    Transmitter(const Properties &props);

    virtual ~Transmitter();

protected:
    /// Combined flags for all properties of this transmitter.
    uint32_t m_flags;
};

MTS_EXTERN_CLASS_RENDER(Transmitter)
NAMESPACE_END(mitsuba)

// -----------------------------------------------------------------------
//! @{ \name Enoki support for vectorized function calls
// -----------------------------------------------------------------------

// Enable usage of array pointers for our types
ENOKI_CALL_SUPPORT_TEMPLATE_BEGIN(mitsuba::Transmitter)
    ENOKI_CALL_SUPPORT_METHOD(sample_ray)
    ENOKI_CALL_SUPPORT_METHOD(doppler)
    ENOKI_CALL_SUPPORT_METHOD(eval)
    ENOKI_CALL_SUPPORT_METHOD(sample_direction)
    ENOKI_CALL_SUPPORT_METHOD(pdf_direction)
    ENOKI_CALL_SUPPORT_METHOD(is_environment)
    ENOKI_CALL_SUPPORT_GETTER(flags, m_flags)
ENOKI_CALL_SUPPORT_TEMPLATE_END(mitsuba::Transmitter)

//! @}
// -----------------------------------------------------------------------
