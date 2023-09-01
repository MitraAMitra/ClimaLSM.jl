
##########
# simulation specification
#########
#set spinup on or off : if on, it use the first year of simulation as spinup
#using FileIO, Images
#img = load("/Users/mitraasadollahi/Projects/CliMA/Data/global_clumping_index_uncompressed.tif")
#time domain
global end_date
global start_date
global spinup
global nspinup
if spinup==1
global t0 = FT((start_date +350*nspinup)* 3600 * 24)# start mid year
else
global t0 = FT(14* 3600 * 24)# start mid year
end
global N_days = end_date-start_date
global tf = t0 + FT(3600 * 24 * N_days)
global dt = FT(300)
global Nt = FT(60)
# soil domain
# For soil column
global nelements = 60#number of grids in the domain
global zmin = FT(-8)# max depth of the domain
global zmax = FT(0)# depth at the surface,
# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
global n_stem = Int64(1);
global n_leaf = Int64(1);
println("numbrer of stem and leaf compartments are  "*string(n_stem)*" "*string(n_leaf))
global h_stem = FT(1) # m, from Wang et al.
println("height of stem compartment is "*string(h_stem))
global h_leaf = FT(6.5) # m from Wang et al.
println("height of leaf compartment is "*string(h_leaf))

##########
# simulation boundary conditions
#########

##########
# simulation initial conditions
#########

##########
# model parameters constant for all sites
#########
#soil water retention parameters

#photosynthesis
############################
# parameters to be double checked
#############################
capacity = FT(22) # kg/m^2
ϕ = FT(0.6) #The quantum yied of photosystem II : this parameter changes in simulation and as a function of GPP
Ω = FT(0.7)#clumping index can also be a function of effective LAI and number of trees a coefficient that light extiction coefficient is devided by
ψ63 = FT(-3 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa

θj = FT(0.9)#Curvature parameter, a fitting constant to compute  𝐽, this has been reported 0.7 in other works
f = FT(0.015)#Constant factor appearing the dark respiration term
Drel = FT(1.6)# relative diffusivity
Γstar25 = FT(4.275e-5)#mol/mol, CO2 compensation point, it is to compensate for photorespiration, this parameter is constant and similar to other literature
Kc25 = FT(4.049e-4)#mol/mol This is a constant value: half rate constant of carboxylase in leaf at 25 dergree C
#Kc the lower it is the higher the GPP rate will be, as it shows with lower CO2 concentration leaf reachees its maximum potentia
Ko25 = FT(0.2874)#mol/mol This is a constant value:half rate constant of oxylase in leaf at 25 dergree C
#Rubisco can either enter oxygenation process or carboxylation, if oxygenation, O2 is added to Rubisco, the process is named as photo respiration and consumes energy
# a higher value of ko indicates that it takes a higher internal concentration of O2 to reach the photorespiration process
oi = FT(0.209)#Intercellular 𝑂2 concentration unit mol/mol
#sc 5e-6 causes a faster transition
#sc very small causes values of close to 1 for all points
#beta =(1 ) ./( exp.(-sc.*(pc .- pl)) .+ 1);
#beta = (1 + exp(sc*pc)) ./( exp.(-sc.*(pc .- pl)) .+ 1)
To = FT(298.15) #Standard temperature

#but in Ozark site, I couldn't see any effect on the GPP estimates
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
n_layers = UInt64(20)
diff_perc = FT(0.2)
#temperature dependance of parameters in photosynthesis model
# has activation energies for ΔHVcmax andΔHJmax per species  https://onlinelibrary.wiley.com/doi/pdf/10.1046/j.1365-3040.2002.00898.x
ΔHkc = FT(79430)
ΔHko = FT(36380)
ΔHVcmax = FT(58520)
ΔHΓstar = FT(37830)
ΔHJmax = FT(43540)
ΔHRd = FT(46390)

# didn't seem to have effect on model performance
ld = FT(0.5)#Leaf angle distribution no unit, this is Xl in the text book, but it actually does not matter much

Weibull_param = FT(4) # unitless, Holtzman's original c param value
##########
# model parameters varying for different sites
#########
############################
#parameters from dataset
############################
# Soil parameters
soil_ν = FT(0.35) # m3/m3
soil_K_sat = FT(4.1e-7) # m/s, matches Natan
soil_vg_n = FT(1.28) # unitless
soil_vg_α = FT(1.9) # inverse meters
θ_r = FT(0.06) # m3/m3
# TwoStreamModel parameters: IT COMES FROM META DATA CLIMA

α_PAR_leaf = FT(0.11)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.35)
τ_NIR_leaf = FT(0.34)

maxLAI = FT(maximum(LAI_data.LAI))# FT(round( quantile(LAI_data.LAI,0.95),digits=2)) # m2/m2, from Wang et al.
#=
 quantile(filtered_v,0.95)
2.781066285

julia> quantile(filtered_v,0.95)*1.1
3.0591729135000003

julia> maximum(filtered_v)
4.5852633

=#
rooting_depth = FT(7) # from Natan in m

# Conductance Model
g1 = FT(80) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300. %50
g0 = FT(0.0001)#minimum stomatal conductance mol/𝑚2.s default



#Photosynthesis model
#zenith angle changes between 0 to pi
#leaf water potential of well watered is -0.2 to -1 MPa, mid to moderate stress for -1 to -2 MPa, and -3 or -4 for drought resistance plants
# for wet environment plant pc can be -1, for 
global Vcmax25 #it is read from LAI dataset= FT(9e-5) # from Yujie's paper 4.5e-5 , Natan used 9e-5 range of VCmax is 1.e-5 to 1.8e-4 mol/m2.s

############################
# parameters with high impact that are tuned
############################
#beta the moisture stress factor
sc = FT(2e-6)#FT(2e-6)#FT(2e6) # Bonan's book: range of 2-5e-6 Pa^{-1}
pc = FT(-1e6) # Bonan's book: -2e6 (Pa): this is a threshols pressure that if leaf potential drops below it, its photosynthesis will be affected by droughts
SAI = FT(maxLAI/4) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
# Plant Hydraulics and general plant parameters
f_root_to_shoot = FT(1)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM
K_sat_plant = 5e-8 # m/s # seems much too small?
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity

plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
#to be checked
soil_S_s = FT(1e-3) # 1/m, guess

##

#####################################
#
#
# Soil heat transfer parameters; not needed for hydrology only test
# The heat transfer is not important to me so it is ignored
#####################################

# Soil heat transfer parameters; not needed for hydrology only test

κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_air = FT(0.025); #W/m/K
ρp = FT(2700); # kg/m^3

################### 
# parameters not found
#####################################
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
ρc_ds = FT((1 - soil_ν) * 4e6); # J/m^3/K
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.3)
soil_α_NIR = FT(0.5)
z0_m = FT(2) # roughness length
z0_b = FT(0.5) # roughness length
########

##### functions to estimate the remaining secondary parameters
#thermal conductivity of the solid 
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
κ_dry = Soil.κ_dry(ρp, soil_ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, soil_ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, soil_ν, κ_liq);


conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
