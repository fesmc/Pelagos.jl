# Convective adjustment closure for Pelagos.jl.
#
# Uses Oceananigans' ConvectiveAdjustmentVerticalDiffusivity with a large
# diffusivity in unstably stratified layers, mimicking the Rahmstorf (1993)
# scheme used in the original GOLDSTEIN model.

module Convection

using Oceananigans.TurbulenceClosures
using ..Parameters: K_CONV

export convective_adjustment_closure

"""
    convective_adjustment_closure()

Returns an Oceananigans `ConvectiveAdjustmentVerticalDiffusivity` with
κ_convective = K_CONV (default 100 m² s⁻¹), which effectively mixes any
unstable density column instantly relative to the ocean dynamics timestep.
"""
function convective_adjustment_closure()
    return ConvectiveAdjustmentVerticalDiffusivity(convective_κz = K_CONV)
end

end # module Convection
