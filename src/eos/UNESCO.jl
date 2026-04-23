# UNESCO equation of state following Millero & Poisson (1981).
# Reference: Millero, F.J. and Poisson, A. (1981). International one-atmosphere
#   equation of state of seawater. Deep-Sea Research, 28A(6), 625–629.
#
# Validated against Table 5 of that paper (see test/eos/test_UNESCO.jl).
# Pressure p is in bar (1 bar = 1e5 Pa). T in °C, S in psu (g/kg).

module UNESCO

export seawater_density, seawater_density_insitu

# ─── Pure-water density at 1 atm, ρ_w(T) ────────────────────────────────────
# Coefficients from Millero & Poisson (1981) Eq. (1)
const A0 =  999.842594
const A1 =    6.793952e-2
const A2 =   -9.095290e-3
const A3 =    1.001685e-4
const A4 =   -1.120083e-6
const A5 =    6.536332e-9

@inline function rho_w(T::Float64)::Float64
    return A0 + T*(A1 + T*(A2 + T*(A3 + T*(A4 + T*A5))))
end

# ─── Salinity-dependent terms at 1 atm ──────────────────────────────────────
# Coefficients from Millero & Poisson (1981) Eq. (1)
const B0 =  8.24493e-1
const B1 = -4.0899e-3
const B2 =  7.6438e-5
const B3 = -8.2467e-7
const B4 =  5.3875e-9

const C0 = -5.72466e-3
const C1 =  1.0227e-4
const C2 = -1.6546e-6

const D0 =  4.8314e-4

@inline function rho_sw_1atm(T::Float64, S::Float64)::Float64
    S32 = sqrt(S)^3  # S^(3/2)
    B   = B0 + T*(B1 + T*(B2 + T*(B3 + T*B4)))
    C   = C0 + T*(C1 + T*C2)
    return rho_w(T) + S*B + S32*C + S^2*D0
end

# ─── Secant bulk modulus K(S,T,p) ────────────────────────────────────────────
# Millero & Poisson (1981) Eqs. (2)–(5); p in bar
const KW0 = 19652.21
const KW1 =   148.4206
const KW2 =    -2.327105
const KW3 =     1.360477e-2
const KW4 =    -5.155288e-5

const KS0 =   54.6746
const KS1 =   -0.603459
const KS2 =    1.09987e-2
const KS3 =   -6.1670e-5

const KSS0 =   7.944e-2
const KSS1 =   1.6483e-2
const KSS2 =  -5.3009e-4

# Pressure terms (at given S,T)
const KP0 =  3.239908
const KP1 =  1.43713e-3
const KP2 =  1.16092e-4
const KP3 = -5.77905e-7

const KPS0 =  2.2838e-3
const KPS1 = -1.0981e-5
const KPS2 = -1.6078e-6

const KPSS0 =  1.91075e-4

const KPP0 =  8.50935e-5
const KPP1 = -6.12293e-6
const KPP2 =  5.2787e-8

const KPPS0 = -9.9348e-7
const KPPS1 =  2.0816e-8
const KPPS2 =  9.1697e-10

@inline function secant_bulk_modulus(T::Float64, S::Float64, p::Float64)::Float64
    S32 = sqrt(S)^3

    # K at zero pressure
    Kw  = KW0 + T*(KW1 + T*(KW2 + T*(KW3 + T*KW4)))
    Ks  = (KS0 + T*(KS1 + T*(KS2 + T*KS3)))*S
    Kss = (KSS0 + T*(KSS1 + T*KSS2))*S32
    K0  = Kw + Ks + Kss

    # First-order pressure terms
    Ap  = KP0 + T*(KP1 + T*(KP2 + T*KP3))
    Aps = (KPS0 + T*(KPS1 + T*KPS2))*S
    Apps = KPSS0 * S32
    A   = Ap + Aps + Apps

    # Second-order pressure terms
    Bp  = KPP0 + T*(KPP1 + T*KPP2)
    Bps = (KPPS0 + T*(KPPS1 + T*KPPS2))*S
    B   = Bp + Bps

    return K0 + p*(A + p*B)
end

# ─── In-situ density ─────────────────────────────────────────────────────────
"""
    seawater_density(T, S, p)

UNESCO (Millero & Poisson 1981) in-situ seawater density.

# Arguments
- `T`  : potential temperature (°C)
- `S`  : salinity (psu / g kg⁻¹)
- `p`  : pressure (bar; 1 bar ≈ 10 m depth)

# Returns
- density ρ (kg m⁻³)
"""
@inline function seawater_density(T::Float64, S::Float64, p::Float64)::Float64
    rho0 = rho_sw_1atm(T, S)
    K    = secant_bulk_modulus(T, S, p)
    return rho0 / (1.0 - p / K)
end

"""
    seawater_density_insitu(T, S, z)

Convenience wrapper: computes pressure from depth `z` (negative downward, m),
then calls `seawater_density`.  Uses hydrostatic approximation p ≈ ρ₀·g·|z|/1e5 bar.
"""
@inline function seawater_density_insitu(T::Float64, S::Float64, z::Float64)::Float64
    p = -z * 1025.0 * 9.81 / 1e5   # bar; z is negative, so p > 0
    return seawater_density(T, S, p)
end

# ─── Vectorised forms for 3D arrays ──────────────────────────────────────────
"""
    seawater_density!(rho, T, S, p)

In-place, element-wise computation of density on arrays.
Arrays must have the same shape `(nlon, nlat, nz)`.
"""
function seawater_density!(rho::AbstractArray{Float64,3},
                           T  ::AbstractArray{Float64,3},
                           S  ::AbstractArray{Float64,3},
                           p  ::AbstractArray{Float64,3})
    @. rho = seawater_density(T, S, p)
    return rho
end

end # module UNESCO
