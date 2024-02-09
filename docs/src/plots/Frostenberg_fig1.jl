import Plots as PL

import CloudMicrophysics as CM
import CloudMicrophysics.HetIceNucleation as IN
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP

FT = Float64
#params = CMP.Frostenberg2023(FT)

T_range = range(-40, stop = -2, length = 500)  # air temperature
INPC_range = 10.0.^(range(-5,stop=7,length=500)) #ice nucleating particle concentration

frequency = [IN.INP_concentration_frequency(INPC,T) > 0.015 ? IN.INP_concentration_frequency(INPC,T) : missing for INPC in INPC_range, T in T_range]
mu = [exp(log(-(1*T)^9 * 10^(-9))) for T in T_range]


PL.contourf(
    T_range,
    INPC_range,
    frequency,
    xlabel = "T [°C]",
    ylabel = "INPC [m⁻³]",
    colorbar_title = "frequency",
    yaxis = :log,
    color = :lighttest,
    gridlinewidth = 3,
    ylims = (1e-5,1e7),
    lw=0,
)

PL.plot!(
    T_range,
    mu,
    label = "median INPC",
    legend = :topright,
    color = :darkred,
)

PL.plot!(
    repeat([-16.], 50),
    10.0.^(range(-2,stop=6,length=50)),
    label = "T = -16°C",
    color = :darkgreen,
    linestyle = :dash,
)
PL.savefig("Frostenberg_fig1.svg")


#plotting the distribution for T=-16°C

T = -16.
INPC_range = 10.0.^(range(-1,stop=4,length=100))
frequency = [IN.INP_concentration_frequency(INPC,T) for INPC in INPC_range]

PL.plot(
    INPC_range,
    frequency,
    xlabel = "INPC [m⁻³]",
    ylabel = "frequency",
    xaxis = :log,
    color = :darkgreen,
    legend = false,
)

PL.savefig("Frostenberg_fig1_T16.svg")