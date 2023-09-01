#=t0 = FT(140 * 3600 * 24)#start day
N_days = 140#total number of days to run the simulation
tf = t0 + FT(3600 * 24 * N_days)# end day
dt = FT(225)#time step in seconds
n = 16#interval of saving data is n*dt
saveat = Array(t0:(n * dt):tf)
=#
t0 = FT(140 * 3600 * 24)#start day
N_days = 360*2#total number of days to run the simulation
tf = t0 + FT(3600 * 24 * N_days)# end day
dt = FT(1000)#time step in seconds
n = 20#interval of saving data is n*dt
saveat = Array(t0:(n * dt):tf)
timestepper = CTS.ARS222()
norm_condition = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition)
max_iterations = 50
#=global t0 
global tf 
global dt
global Nt 
println(dt)
println(Nt)
println(t0)
println(tf)
println("time step is ")
println(dt)
println("t0 is ")
println(t0)
println("Nt is ")
println(Nt)
println("tf is ")
println(tf)
saveat = Array(t0:(Nt * dt[1]):tf)
timestepper = CTS.ARS222()
norm_condition = CTS.MaximumError(FT(1e-7))
conv_checker = CTS.ConvergenceChecker(; norm_condition)
max_iterations = 20=#
