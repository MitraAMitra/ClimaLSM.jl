using DiffEqCallbacks
using DataFrames
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation
using CSV
using ClimaLSM
using ClimaLSM.Domains: LSMSingleColumnDomain
using ClimaLSM.Soil
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP
#climalsm_dir= pkgdir(ClimaLSM)
climalsm_dir="/Users/mitraasadollahi/Projects/CliMA/ClimaLSM.jl"
include(joinpath(climalsm_dir, "parameters", "create_parameters.jl"))
global start_year=4#minimum(daily)
global end_year=6#maximum(daily)
global start_date=1+(start_year-1)*365#maximum(daily)
global end_date=360+(end_year-1)*365#maximum(daily)
site_name="US-Var";
global datadir="/Users/mitraasadollahi/Projects/CliMA/Data/FLUXNETSites/data_processed_L1/FLX_"*site_name*"/FLX_"*site_name*"_FLUXNET2015_FULLSET_HH.csv"
savedir = joinpath(climalsm_dir, "experiments/integrated/"*site_name*"/hydrology_only/result/")
curdir = joinpath(climalsm_dir, "experiments/integrated/"*site_name*"/")
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
df = CSV.File("/Users/mitraasadollahi/Projects/CliMA/Data/fluxnet_site_location.csv") |> DataFrame

# Filter the dataframe for the specific site_name (assuming you have a variable `site_name` with the desired name)
site_row = df[df.Site .== site_name, :]

# Extract Longitude and Latitude
global long = FT(site_row[1, :Longitude])
global lat = FT(site_row[1, :Latitude])
global spinup=0 #allows 1 year of spinup
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/"*site_name*"/met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/"*site_name*"/parameters.jl"),
)
include(joinpath(climalsm_dir, "experiments/integrated/"*site_name*"/domain.jl"))

include(
    joinpath(climalsm_dir, "experiments/integrated/"*site_name*"/hydrology_only/simulation.jl"),
)
#####################
# SECTION 0.1: setting BCs and picking models
#####################
# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain.subsurface
soil_ps = Soil.RichardsParameters{FT, vanGenuchten{FT}}(
    soil_ν,
    vanGenuchten(; α = soil_vg_α, n = soil_vg_n),
    soil_K_sat,
    soil_S_s,
    θ_r,
)

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.RichardsModel{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    radiative_transfer = Canopy.BeerLambertModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
)
# Individual Component arguments
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = BeerLambertParameters{FT}(;
        Ω = Ω,
        ld = ld,
        α_PAR_leaf = α_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        λ_γ_PAR = λ_γ_PAR,
        λ_γ_NIR = λ_γ_NIR,
    )
)
# Set up conductance
conductance_args = (;
    parameters = MedlynConductanceParameters{FT}(;
        g1 = g1,
        Drel = Drel,
        g0 = g0,
    )
)
# Set up photosynthesis
photosynthesis_args = (;
    parameters = FarquharParameters{FT}(
        Canopy.C3();
        oi = oi,
        ϕ = ϕ,
        θj = θj,
        f = f,
        sc = sc,
        pc = pc,
        Vcmax25 = Vcmax25,
        Γstar25 = Γstar25,
        Kc25 = Kc25,
        Ko25 = Ko25,
        To = To,
        ΔHkc = ΔHkc,
        ΔHko = ΔHko,
        ΔHVcmax = ΔHVcmax,
        ΔHΓstar = ΔHΓstar,
        ΔHJmax = ΔHJmax,
        ΔHRd = ΔHRd,
    )
)
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)

plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

# Canopy component args
canopy_component_args = (;
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
)
# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = land_domain.surface)

# Integrated plant hydraulics and soil model
# If you wish to pass in a prescribed transpiration function (for
# debugging purposes), you would pass it as another element
# of this tuple. The default if it is not passed is
# a diagnostic (big leaf-computed) transpiration.
soil_driver = PrognosticSoil(; soil_α_PAR = soil_α_PAR, soil_α_NIR = soil_α_NIR)
land_input = (atmos = atmos, soil_driver = soil_driver, radiation = radiation)
land = SoilPlantHydrologyModel{FT}(;
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)

#Initial conditions
Y.soil.ϑ_l = SWC[Int(round(t0 / 1800))] # Get soil water content at t0
ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        plant_ν,
        plant_S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l[i] .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

set_initial_aux_state! = make_set_initial_aux_state(land)
set_initial_aux_state!(p, Y, t0);

# Set up timestepper and jacobian for soil
update_jacobian! = make_update_jacobian(land.soil)

ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = max_iterations,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

W = RichardsTridiagonalW(Y)
jac_kwargs = (; jac_prototype = W, Wfact = update_jacobian!)

# Simulation
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)

prob = ODE.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = ODE.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLSM.dss!,
    ),
    Y,
    (t0, tf),
    p,
)
@time sol = ODE.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
)

T =[
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
#E =
#    [parent(sv.saveval[k].soil_evap)[1] for k in 1:length(sol.t)] .* (1e3 * 24 * 3600)

β = [parent(sv.saveval[k].canopy.hydraulics.β)[1] for k in 1:length(sol.t)]

daily = sol.t ./ 3600 ./ 24
root_stem_flux = [
    sum(sv.saveval[k].root_extraction) .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]

stem_leaf_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[1] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
leaf_air_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[2] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
lwp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[2] * 9800 for k in 1:length(sol.t)
]
swp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[1] * 9800 for k in 1:length(sol.t)
]
ψ_soil = [
    sum(
        parent(sv.saveval[k].soil.ψ) .*
        root_distribution.(parent(cds.subsurface.z)),
    ) / sum(root_distribution.(parent(cds.subsurface.z))) * 9800 for
    k in 1:length(sol.t)
]
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in 1:length(sv.saveval)
]
swc= [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)];

g_stomata =
    [parent(sv.saveval[k].canopy.conductance.gs)[1] for k in 1:length(sol.t)]
  
measured_ET = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)

# Calculate the combined E and T values


#######################
# SECTION 3: plotting results
########################


# Plotting


if !isdir(savedir)
    mkdir(savedir)
    println("Directory $savedir created.")
else
    println("Directory $savedir already exists.")
end
maxGPP=maximum([maximum(GPP .* 1e6),maximum(model_GPP .* 1e6)])
plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    daily,
    model_GPP .* 1e6,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0, maxGPP],
    title = "GPP [mu mol/m^2/s]",
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    GPP .* 1e6,
    label = "Data",
    ylim = [0, maxGPP],
    lalpha = 0.3,
)

plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    daily,
    model_GPP .* 1e6,
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    ylim = [0, maximum(model_GPP .* 1e6)],
    margin = 10Plots.mm,
    xlabel = "Day of year",
)
Plots.plot!(plt2, seconds ./ 3600 ./ 24, GPP .* 1e6, label = "", lalpha = 0.3)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "GPP.png"))

# SW_IN

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(SW_IN),
    label = "Data",
    title = " SWd (W/m^2)",
    xlim = [minimum(daily), maximum(daily)],
)

plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    FT.(SW_IN),
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    margin = 10Plots.mm,
    xlabel = "Day of year",
)
Plots.plot(plt1, plt2, layout = (2, 1))

Plots.savefig(joinpath(savedir, "SW_IN.png"))

# VPD

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(VPD),
    label = "Data",
    title = "VPD (Pa)",
    xlim = [minimum(daily), maximum(daily)],
)

plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    FT.(VPD),
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    ylim = [0, maximum(FT.(VPD))+10],
    margin = 10Plots.mm,
    xlabel = "Day of year",
)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "VPD.png"))




# Create the first plot to compare combined_ET against measured_ET
plt1 = Plots.plot(seconds ./ 3600 ./ 24, measured_ET, label="Data ET", margins=10Plots.mm, size = (1500, 400))

# For distinction, you can use the fillrange attribute to visually separate E and T
Plots.plot!(plt1, daily, T, fillrange=0, label="Model T", alpha=0.5)
#Plots.plot!(plt1, daily, E, fillrange=0, label="Model E", alpha=0.5)

# Create the second plot (zoomed in version)
plt2 = Plots.plot(seconds ./ 3600 ./ 24, measured_ET, label="", size=(1500, 400))
Plots.plot!(plt2, daily, T, fillrange=0, label="", alpha=0.5)

# Display the plots
Plots.plot(plt1, plt2, layout = (2, 1))


Plots.savefig(joinpath(savedir, "ET.png"))


plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    daily,
    β,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Moisture stress factor",
)
#i_week = Int(round((7 * 24 * 3600) / (sol.t[2] - sol.t[1])))
plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    daily,
    β,
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    xlabel = "Day of year",
    #ylim = [minimum(β[1:i_week]), maximum(β[1:i_week])],
    margin = 10Plots.mm,
)

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "moisture_stress.png"))


plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    daily,
    g_stomata,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Stomatal conductance (mol/m^2/s)",
)
i_week = Int(round((7 * 24 * 3600) / (sol.t[2] - sol.t[1])))
plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    daily,
    g_stomata,
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    xlabel = "Day of year",
    ylim = [minimum(g_stomata[1:i_week]), maximum(g_stomata[1:i_week])],
    margin = 10Plots.mm,
)

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "stomatal_conductance.png"))

# Current resolution is 3.333 cm per layer, our nodes are at the center of the
# layers. The second layer is ~ 5cm
plt1 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "5cm",
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0.05, 0.55],
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    color = "blue",
    margin = 10Plots.mm,
)



Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC, label = "Data")
plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    P .* (-1e3 * 24 * 3600),
    label = "Data",
    ylabel = "Precipitation [mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    margin = 10Plots.mm,
    ylim = [-200, 0],
    size = (1500, 400),
)
Plots.plot(plt2, plt1, layout = grid(2, 1, heights = [0.2, 0.8]))
Plots.savefig(joinpath(savedir, "soil_water_content.png"))

plt1 = Plots.plot(size = (1500, 700))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(H),
    label = "Data",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(H_CORR),
    label = "Data Corr",
    margins = 10Plots.mm,
)


plt1 = Plots.plot(size = (1000, 400))
dt_model = sol.t[2] - sol.t[1]
dt_data = seconds[2] - seconds[1]
# Find which index in the data our simulation starts at:
idx = argmin(abs.(seconds .- sol.t[1]))
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(measured_ET[:]) * dt_data,
    label = "Data ET",
)
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(P[:]) * dt_data * (1e3 * 24 * 3600),
    label = "Data P",
)
Plots.plot!(
    daily,
    cumsum(T) * dt_model .+ cumsum(measured_ET[:])[idx] * dt_data,
    label = "Model T",
)

Plots.plot!(ylabel = "∫ Water fluxes dt", xlabel = "Days", margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "cumul_p_et.png"))


plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    daily,
    lwp,
    label = "Model, Leaf",
    title = "Water potentials",
    xlim = [minimum(daily), maximum(daily)],
)
Plots.plot!(plt1, daily, swp, label = "Model, Stem")
Plots.plot!(plt1, daily, ψ_soil, label = "Model, Mean soil")

plt2 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt2,
    daily,
    lwp,
    label = "",
    xlim = [minimum(daily), minimum(daily) + 30],
    xlabel = "Day of year",
    ylim = [-2e6, 0],
    margin = 10Plots.mm,
)

Plots.plot!(plt2, daily, swp, label = "")
Plots.plot!(plt2, daily, ψ_soil, label = "")

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "leaf_water_potential.png"))

#################################
# save a data file csv 
###################################
#
sim_df = DataFrame(
    doy = daily,
    leaf_water_potential=lwp,
    model_GPP = model_GPP .* 1e6,
    T_model = T,
    soil_moisture_stress_factor = β,
    g_stomata = g_stomata,
    swc=swc,
    soil_water_potential=swp,
    soil_mean_potential=ψ_soil
)
# Simulation
filename = "simulated_data_"*site_name*".csv"
savepath = joinpath(savedir, filename)
CSV.write(savepath, sim_df)