##########
# model parameters constant for all sites
#########
#photosynthesis
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
Ω = FT(0.69)#clumping index can also be a function of effective LAI and number of trees a coefficient that light extiction coefficient is devided by
#but in Ozark site, I couldn't see any effect on the GPP estimates
##########
# model parameters varying for different sites
#########
# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021


# TwoStreamModel parameters
ld = FT(0.5)#Leaf angle distribution no unit
α_PAR_leaf = FT(0.1)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.45)
τ_NIR_leaf = FT(0.25)
n_layers = UInt64(20)
diff_perc = FT(0.2)

# Conductance Model
g1 = FT(60) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300. %50
Drel = FT(1.6)# relative humidity
g0 = FT(0.01)#minimum stomatal conductance mol/𝑚2.s default


#Photosynthesis model
#zenith angle changes between 0 to pi
ϕ = FT(0.6) #The quantum yied of photosystem II 
θj = FT(0.9)#Curvature parameter, a fitting constant to compute  𝐽, this has been reported 0.7 in other works
f = FT(0.015)#Constant factor appearing the dark respiration term
#beta the moisture stress factor
sc = FT(2e-6)#FT(2e6) # Bonan's book: range of 2-5e-6 Pa^{-1}
pc = FT(-3e6) # Bonan's book: -2e6 (Pa): this is a threshols pressure that if leaf potential drops below it, its photosynthesis will be affected by droughts
#leaf water potential of well watered is -0.2 to -1 MPa, mid to moderate stress for -1 to -2 MPa, and -3 or -4 for drought resistance plants
# for wet environment plant pc can be -1, for 
Vcmax25 = FT(9e-5) # from Yujie's paper 4.5e-5 , Natan used 9e-5 range of VCmax is 1.e-5 to 1.8e-4 mol/m2.s
#temperature dependance of parameters in photosynthesis model
# has activation energies for ΔHVcmax andΔHJmax per species  https://onlinelibrary.wiley.com/doi/pdf/10.1046/j.1365-3040.2002.00898.x
ΔHkc = FT(79430)
ΔHko = FT(36380)
ΔHVcmax = FT(58520)
ΔHΓstar = FT(37830)
ΔHJmax = FT(43540)
ΔHRd = FT(46390)

# Plant Hydraulics and general plant parameters
SAI = FT(1.0) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
maxLAI = FT(4.2) # m2/m2, from Wang et al.
f_root_to_shoot = FT(3.5)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
capacity = FT(10) # kg/m^2
plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(0.5) # from Natan in m
z0_m = FT(2)
z0_b = FT(0.2)

#####################################
#
#
# Soil heat transfer parameters; not needed for hydrology only test
# The heat transfer is not important to me so it is ignored
#####################################

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_air = FT(0.025); #W/m/K
ρp = FT(2700); # kg/m^3
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
κ_dry = Soil.κ_dry(ρp, soil_ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, soil_ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, soil_ν, κ_liq);
ρc_ds = FT((1 - soil_ν) * 4e6); # J/m^3/K
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.3)
soil_α_NIR = FT(0.5)

conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)