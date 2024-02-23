import CloudMicrophysics.Parameters as CMP
import Thermodynamics.Parameters as TDP

struct Empty{FT} <: CMP.ParametersType{FT} end

struct MohlerAF{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    aerosol::CMP.AerosolType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    const_dt::FT
end

struct MohlerRate{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    aerosol::CMP.AerosolType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end

struct ActivityBased{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
    r_nuc::FT
end

struct P3_dep{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct ABIFM{FT} <: CMP.ParametersType{FT}
    H₂SO₄ps::CMP.H2SO4SolutionParameters{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    aerosol::CMP.AerosolType{FT}
end

struct P3_het{FT} <: CMP.ParametersType{FT}
    ips::CMP.ParametersType{FT}
    const_dt::FT
end

struct WaterAct{FT} <: CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
    ips::CMP.ParametersType{FT}
    Jensen::Bool
end

struct P3_hom{FT} <: CMP.ParametersType{FT}
    const_dt::FT
end

struct Frostenberg{FT} <: CMP.ParametersType{FT}
    sigma::CMP.ParametersType{FT}
    drawing_interval::FT
    using_mean::Bool
end

struct CondParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end

struct DepParams{FT} <: CMP.ParametersType{FT}
    aps::CMP.ParametersType{FT}
    tps::TDP.ThermodynamicsParameters{FT}
end
