using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates
using Statistics
function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end
function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end
function find_site_name_dirs(root_directory,site_name)
    directories_with_site_name = String[]

    for (root, dirs, files) in walkdir(root_directory)
        for dir in dirs
            if occursin(site_name, dir)
                push!(directories_with_site_name, joinpath(root_directory, dir))
            end
        end
    end

    return directories_with_site_name
end
function find_csv_file_with_name(directory)
    for file in readdir(directory)
        path = joinpath(directory, file)
        if isfile(path) && occursin("FULLSET_HH", file) && endswith(file, ".csv")
            return path
        elseif isdir(path)
            result = find_csv_file_with_name(path)
            if result != nothing
                return result
            end
        end
    end
    return nothing
end
global site_name
global start_year
global end_year
global lat 
global long
global spinup
#read fluxnet dataset and LAI from modis, if LAI is absent remove the corresponding dates from  FLUXNET sites
root_dir_all_FLUXNET2015="/Users/mitraasadollahi/Projects/CliMA/Data/FLUXNETSites/data_processed_L1/"
root_data_files=find_site_name_dirs(root_dir_all_FLUXNET2015,site_name)
csv_file_path = find_csv_file_with_name(root_data_files[1])
driver_data = readdlm(csv_file_path, ',')
column_names = driver_data[1, :]
dates= driver_data[2:end, column_names .== "TIMESTAMP_START"]; # time yyyymmddhhmm
date_strings = vec([string(date) for date in dates])
shortened_date_strings = [date[1:4] for date in date_strings]
year=parse.(FT, shortened_date_strings)
unique_year=unique(year)
first_year=unique_year[start_year]
last_year=unique_year[end_year]
filtered_indices = findall(year .>= first_year .&& year .<= last_year).+1;
dates= driver_data[filtered_indices, column_names .== "TIMESTAMP_START"]; #
filtered_dates_date = DateTime.(string.(dates), "yyyymmddHHMM")

 first_date = Dates.format(filtered_dates_date[1], "yyyy-mm-dd-HH")#"2004-01-01-06"

if spinup.==1
    global nspinup=1#number of years spinup occurs
    spinup_indices = findall(year .>= year[1] .&& year .<= year[1+nspinup]).+1;
    spinup_indices_dates=driver_data[spinup_indices, column_names .== "TIMESTAMP_START"]; #
    spinup_dates = DateTime.(string.(spinup_indices_dates), "yyyymmddHHMM")
    driver_data =vcat(driver_data[1:1, :], driver_data[spinup_indices, :],driver_data[filtered_indices, :])
else
    filtered_indices=[1;filtered_indices]
    driver_data = driver_data[filtered_indices, :];#only selected range of data is used
end
column_names = driver_data[1, :]
####
#air temperature
TA = vec(driver_data[2:end, column_names .== "TA_F"] .+ 273.15); # convert C to K
#air aridity
VPD = driver_data[2:end, column_names .== "VPD_F"] .* 100; # convert hPa to Pa
#air pressure
PA = driver_data[2:end, column_names .== "PA_F"] .* 1000; # convert kPa to Pa
#precipitation
P = driver_data[2:end, column_names .== "P_F"] ./ (1000 * 1800); # convert mm/HH to m/s
#wind speed 
WS = driver_data[2:end, column_names .== "WS_F"]./1; # already m/s
#incoming long wave radiation
LW_IN = driver_data[2:end, column_names .== "LW_IN_F"]
#incoming short wave radiation
SW_IN = driver_data[2:end, column_names .== "SW_IN_F"]
#CO2 mole fraction
CO2_F = driver_data[2:end, column_names .== "CO2_F_MDS_QC"]
CO2 = driver_data[2:end, column_names .== "CO2_F_MDS"] .* 1e-6; # convert \mumol to mol
replace_missing_with_mean!(CO2, CO2_F)
#soil moisture
SWC_F = driver_data[2:end, column_names .== "SWC_F_MDS_1_QC"] # Most likely 5cm depth
SWC = driver_data[2:end, column_names .== "SWC_F_MDS_1"] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC, SWC_F)
SWC=vec(SWC)
#soil temperature
TS_F = driver_data[2:end, column_names .== "TS_F_MDS_1_QC"]./1 # Most likely 5cm depth
TS = driver_data[2:end, column_names .== "TS_F_MDS_1"] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS, TS_F)
#gross primary production
GPP = driver_data[2:end, column_names .== "GPP_DT_VUT_REF"] .* 1e-6 # to convert from micromol to mol.
#latent heat
LE = driver_data[2:end, column_names .== "LE_CORR"]./1
#Sensible heat flux multiple sources
H_CORR = driver_data[2:end, column_names .== "H_CORR"]./1
H = driver_data[2:end, column_names .== "H_F_MDS"]./1
H_F = driver_data[2:end, column_names .== "H_F_MDS_QC"]
replace_missing_with_mean!(H, H_F)
#soil heat flux
G = driver_data[2:end, column_names .== "G_F_MDS"]./1
G_F = driver_data[2:end, column_names .== "G_F_MDS_QC"]
replace_missing_with_mean!(G, G_F)
#outgoing long and short wave radiation
LW_OUT = driver_data[2:end, column_names .== "LW_OUT"]# This has missing data
SW_OUT = driver_data[2:end, column_names .== "SW_OUT"]# This has missing data
LOCAL_DATETIME = DateTime.(string.(driver_data[2:end, 1]), "yyyymmddHHMM")
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(8)
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )

e = @. esat - VPD
q = @. 0.622 * e ./ (PA - 0.378 * e)

#Make a bunch of splines
seconds = FT.(0:1800:((length(UTC_DATETIME) - 1) * 1800));
p_spline = Spline1D(seconds, -vec(P[:])) # m/s
atmos_q = Spline1D(seconds, q[:])
atmos_T = Spline1D(seconds, vec(TA[:]))
atmos_p = Spline1D(seconds, vec(PA[:]))
atmos_co2 = Spline1D(seconds, vec(CO2[:]))
atmos_u = Spline1D(seconds, vec(WS[:]))
LW_IN_spline = Spline1D(seconds, vec(LW_IN[:]))
SW_IN_spline = Spline1D(seconds, vec(SW_IN[:]))
atmos_h = FT(32)
precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = eltype(t)(0) # this is likely not correct


# Construct the drivers
atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    atmos_h;
    c_co2 = atmos_co2,
)

#lat = FT(38.7441) # degree
#long = FT(-92.2000) # degree
datetime_obj = DateTime(first_date, "yyyy-mm-dd-HH")

# Add 8 hours
new_datetime = datetime_obj + Hour(8)#it should be UTC time, so for USA sites it is +8 hours
function zenith_angle(
    t::FT,
    orbital_data;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC
    dt =new_datetime+ Dates.Second(round(t))
    FT(
        instantaneous_zenith_angle(
            dt,
            orbital_data,
            longitude,
            latitude,
            insol_params,
        )[1],
    )
end


radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline;
    θs = zenith_angle,
    orbital_data = Insolation.OrbitalData(),
)

# LAI data
lai_dataset_path ="/Users/mitraasadollahi/Projects/CliMA/Data/LAI/";
lai_raw_data = joinpath(lai_dataset_path, site_name*"_fSun_Vcmax_DD.csv");
LAI_data = CSV.File(lai_raw_data) |> DataFrame

tmp=[LAI_data[:,:Vcmax25_leaf_v1];LAI_data[:,:Vcmax25_leaf_v2]];
nan_vec = replace(tmp, "NA" => "NaN") #global_clumping_index_uncompressed.tif
parsed_values=parse.(Float64, nan_vec)
parsed_values = [isnan(val) ? missing : val for val in parsed_values]
global  Vcmax25=quantile(skipmissing(parsed_values),0.95)*10^-6*1.1;#minimum(skipmissing(parsed_values))*10^-6# mol/m2.s
nan_vec = replace(LAI_data[:,"LAI"], "NA" => "NaN")
tmp1=string.(LAI_data[:,"TIMESTAMP"])
tmptime = map(x -> DateTime(x, "yyyymmddHH"), tmp1)
LAI_data=DataFrame([tmptime,parse.(Float64, nan_vec)],[:time_l,:LAI])
#take monthly average of LAI across all years and use it to gapfill the data
LAI_data[!, :month] = month.(LAI_data.time_l)
monthly_means = combine(groupby(LAI_data, :month)) do sdf
    mean_value = mean(filter(x -> !isnan(x), sdf.LAI))
    DataFrame(mean_LAI = mean_value)
end
#join back mean_LAI on the LAI data
LAI_data = leftjoin(LAI_data, monthly_means, on = :month)
# Replace NaN in LAI with the corresponding monthly mean
LAI_data[!, :LAI] = ifelse.(isnan.(LAI_data[!, :LAI]), LAI_data[!, :mean_LAI], LAI_data[!, :LAI])
mask = [t in filtered_dates_date for t in LAI_data.time_l]
LAI_filtered = LAI_data[mask, "LAI"]

if spinup.==1
    spinup_indices_dates = DateTime.(string.(spinup_indices_dates), "yyyymmddHHMM")
    mask2 = [t in spinup_dates for t in LAI_data.time_l]
    spinup_LAI=LAI_data[mask2, "LAI"]
    LAI_repeated =vcat(repeat(spinup_LAI, inner = 48),repeat(LAI_filtered, inner = 48))
else
LAI_repeated = repeat(LAI_filtered, inner = 48)
end

# This has the same timestamp as the driver data, so it's ok to use the time column from that file here
LAIspline = Spline1D(seconds, LAI_repeated[:])
LAIfunction = (t) -> eltype(t)(LAIspline(t))
