global t0 
global tf 
global dt
global Nt 
println(dt)
println(Nt)
println(t0)
println(tf)
saveat = Array(t0:(Nt * dt[1]):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)
