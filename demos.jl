include("ising2d.jl")
include("graphs.jl")

# show nearest neighbours
function demo_0()
	sys = Ising2D(L=20)
	particles = rand(1:20^2, 10)
	for particle in particles
		clear_screen()
		show_nearest_neighbours(particle, sys)
		sleep(1)
	end
end

# watch it equilibriate
function demo_1()
	sys = Ising2D(
		L=40, 
		ltype=all_up,
		H=0.0
	)

	equilibriate(sys, T₀=5.0, numIterations=10000)

	p1 = plot_m_t(sys)

	plot(p1)
end

# temperature sweep (zero external field)
function demo_2()
	sys = Ising2D(
		L=40, 
		ltype=random,
		H=0.0
	)

	sweep(sys, T₀=5.0)

	p1 = plot_m_t(sys)

	plot(p1)
end
