import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Initial conditions
ρₗ = wps.ρw
Nₐ = FT(0)
Nₗ = FT(500 * 1e3)
Nᵢ = FT(0)
r₀ = FT(1e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(251)
qᵥ = FT(8.1e-4)
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2) # 1.2 should be ρₐ
qᵢ = FT(0)
x_sulph = FT(0.01)

# Moisture dependent initial conditions
q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
R_v = TD.Parameters.R_v(tps)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = FT(0.5 * 1.e-4 * 1e-6)             # assumed size of nucleated particles
w = FT(0.7)                                # updraft speed
const_dt = FT(1)                           # model timestep
t_max = FT(1200)
aerosol = CMP.Illite(FT)
heterogeneous = "ImmersionFreezing"
growth_modes = ["Condensation", "Deposition"]
size_distribution_list = ["Monodisperse", "Gamma"]

# Plotting
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Supersaturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_ice [g/kg]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1], xlabel = "Time [min]", ylabel = "N_liq")
ax6 = MK.Axis(
    fig[3, 2],
    xlabel = "Time [min]",
    ylabel = "N_ice / N_tot",
    yscale = log10,
)
MK.ylims!(ax6, 1e-6, 1)

for DSD in size_distribution_list
    params = parcel_params{FT}(
        const_dt = const_dt,
        r_nuc = r_nuc,
        w = w,
        aerosol = aerosol,
        heterogeneous = heterogeneous,
        growth_modes = growth_modes,
        size_distribution = DSD,
    )
    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, params)

    # Plot results
    MK.lines!(ax1, sol.t / 60, S_i.(tps, sol[3, :], sol[1, :]) .- 1, label = DSD)
    MK.lines!(ax2, sol.t / 60, sol[3, :])
    MK.lines!(ax3, sol.t / 60, sol[6, :] * 1e3)
    MK.lines!(ax4, sol.t / 60, sol[5, :] * 1e3)
    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    MK.lines!(ax5, sol.t / 60, sol_Nₗ)
    MK.lines!(ax6, sol.t / 60, sol_Nᵢ ./ (sol_Nₗ .+ sol_Nᵢ))
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)

MK.save("immersion_freezing.svg", fig)
