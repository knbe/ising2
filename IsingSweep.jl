# ising 2d.

using Plots

mutable struct Ising2D
	L::Int64
	N::Int64
	S::Vector{Int64}
	H::Float64
	J::Float64
	coords::Matrix{Float64}
	steps::Int64
	E::Vector{Float64}
	M::Vector{Float64}
	T::Vector{Float64}
end

@enum LatticeType begin
	random
	all_up
	all_down
	icing
end

function Ising2D(;L::Int64=40, ltype::LatticeType=random, 
		H::Float64=0.0, J::Float64=1.0)
	N, S, coords = initialise(L, Int(ltype))
	steps = 0
	E = Float64[]
	M = Float64[]
	T = Float64[]
	return Ising2D(L, N, S, H, J, coords, steps, E, M, T)
end

@views lpoint(n,sys) = sys[n,:]

function clear_screen()
	print("\033c")
end

function show_spins(sys::Ising2D)
	n = 1
	for j = 1:sys.L
		for i = 1:sys.L
			sys.S[n] < 0 ? print("░░") : print("██")
			n += 1
		end
		#@printf "\n"
		print("\n")
	end
	#@printf "\n"
	print("\n")
end

function print_params(sys::Ising2D)
	println(
		"T=$(round(sys.T[end], digits=4)), \
		L=$(sys.L), \
		N=$(sys.N), \
		H=$(sys.H), \
		J=$(sys.J), \
		m=$(sys.M[end] / sys.N), \
		E=$(sys.E[end])"
	)
end

function print_equilibriation_params(sys::Ising2D, numIter::Int64)
	println(
		"T=$(round(sys.T[end], digits=4)), \
		L=$(sys.L), \
		N=$(sys.N), \
		H=$(sys.H), \
		J=$(sys.J), \
		iterations=$(numIter), \
		m=$(round((sys.M[end] / sys.N), digits=4)), \
		E=$(sys.E[end])"
	)
end

function initialise(L::Int64, ltype::Int64)
	N = L^2
	S = ones(Int64, N)
	coords = zeros(N, 2)

	# random
	if ltype == 0
		for i = 1:N
			if rand() < 0.5
				S[i] = 1
			else
				S[i] = -1
			end
		end
	# allup
	elseif ltype == 1
		for i = 1:N
			S[i] = 1
		end
	# alldown
	elseif ltype == 2
		for i = 1:N
			S[i] = -1
		end
	# silly
	elseif ltype == 3
	end

	n = 1
	for j = 1:L
		for i = 1:L
			coords[n,:] = [float(i), float(j)]
			n += 1
		end
	end

	return N, S, coords
end

# returns indices of the 4 nearest neighbours of the nth particle
function get_nearest_neighbours(sys::Ising2D, n::Int64)
	coords = lpoint(n, sys.coords)

	(coords[2]-1) < 1 ? n_up = n + sys.L*(sys.L-1) : n_up = n - sys.L
	(coords[2]+1) > sys.L ? n_dn = n - sys.L*(sys.L-1) : n_dn = n + sys.L
	(coords[1]-1) < 1 ? n_lf = n + (sys.L-1) : n_lf = n - 1
	(coords[1]+1) > sys.L ? n_rt = n - (sys.L-1) : n_rt = n + 1

	return n_up, n_dn, n_lf, n_rt
end

function show_nearest_neighbours(n::Int64, sys::Ising2D)

	n_up, n_dn, n_lf, n_rt = get_nearest_neighbours(sys,n)

	clear_screen()

	println(
		"L=$(sys.L), \
		N=$(sys.N), \
		H=$(sys.H)"
	)

	println("showing particle:")
	print("\033[0;31m")
	println("\tn=$n, coords=$(sys.coords[n,:])")
	print("\033[0;37m")
	println("and nearest neighbours:")
	print("\033[1;31m")
	println("\tn_up=$(n_up), coords=$(sys.coords[n_up,:])")
	println("\tn_dn=$(n_dn), coords=$(sys.coords[n_dn,:])")
	println("\tn_lf=$(n_lf), coords=$(sys.coords[n_lf,:])")
	println("\tn_rt=$(n_rt), coords=$(sys.coords[n_rt,:])")

	print("\033[0;37m")
	index = 1
	for j = 1:sys.L
		for i = 1:sys.L
			if index == n
				print("\033[0;31m")
			elseif index == n_up || index == n_dn || 
				index == n_lf || index == n_rt
				print("\033[1;31m")
			end
			sys.S[n] < 0 ? print("▓▓") : print("██")
			index += 1
			print("\033[0;37m")
		end
		print("\n")
	end
	print("\n")
end

function ssflip_ΔE(sys::Ising2D, n::Int64)
	n_up, n_dn, n_lf, n_rt = get_nearest_neighbours(sys, n)

	return 2 * sys.J * sys.S[n] * (sys.S[n_up] +
		sys.S[n_dn] + sys.S[n_lf] + sys.S[n_rt])
end

function ss_interaction(sys::Ising2D, n::Int64)
	n_up, n_dn, n_lf, n_rt = get_nearest_neighbours(sys, n)

	return sys.S[n] * (sys.S[n_up] +
		sys.S[n_dn] + sys.S[n_lf] + sys.S[n_rt])
end

function initial_statistics!(sys::Ising2D)
	E = 0.0
	M = 0.0
	for i = 1:sys.N
		E += -sys.J * ss_interaction(sys,i) - sys.H * sys.S[i]
		M += sys.S[i]
	end
	push!(sys.E, E)
	push!(sys.M, M)
end

function update_statistics!(sys::Ising2D)
end

function equilibriate(sys::Ising2D; T₀::Float64=5.0, numIterations::Int64=10000)
	T = T₀
	push!(sys.T, T₀)
	initial_statistics!(sys)

	clear_screen()
	print_params(sys)
	show_spins(sys)
	sleep(5)


	for iteration in 1:numIterations
		n = rand(1:sys.N)

		dE = 0.0
		dM = 0.0

		ΔE = ssflip_ΔE(sys, n)

		if ΔE < 0
			sys.S[n] *= -1
			dE += sys.S[n]
			dM += 2 * sys.S[n]
		else
			if rand() < exp(-ΔE / T)
				sys.S[n] *= -1
				dE += sys.S[n]
				dM += 2 * sys.S[n]
			end
		end

		push!(sys.E, sys.E[end] + dE)
		push!(sys.M, sys.M[end] + dM)
		push!(sys.T, T)
		sys.steps = iteration

		if iteration % 10 == 0
			clear_screen()
			#print_params(sys)
			print_equilibriation_params(sys, iteration)
			show_spins(sys)
			sleep(0.01)
		end
	end
end

function sweep(sys::Ising2D; T₀::Float64=5.0)
	T = T₀
	Tstep = 0.005

	numMCIterations = 1000000
	#numTsteps = Int(T/Tstep)

	initial_statistics!(sys)
	push!(sys.T, T₀)

	clear_screen()
	print_params(sys)
	show_spins(sys)
	sleep(5)

	while T > 0.0
		dE = 0.0
		dM = 0.0
		for iteration in 1:numMCIterations
			n = rand(1:sys.N)

			ΔE = ssflip_ΔE(sys, n)

			if ΔE < 0
				sys.S[n] *= -1
				dE += sys.S[n]
				dM += 2 * sys.S[n]
			else
				if rand() < exp(-ΔE / T)
					sys.S[n] *= -1
					dE += sys.S[n]
					dM += 2 * sys.S[n]
				end
			end
		end
		push!(sys.E, sys.E[end] + dE)
		push!(sys.M, sys.M[end] + dM)
		push!(sys.T, T)

		if sys.steps % 2 == 0
			clear_screen()
			print_params(sys)
			show_spins(sys)
		end

		T -= Tstep
		sys.steps += 1
	end
end

# GRAPHS
################################################################################

function initialize_plot()
	plot(
		size=(800,800), 
		titlefontsize=12, 
		guidefontsize=12,
	)
end

function plot_m_T(sys::Ising2D)
	initialize_plot()
	plot!(
		sys.T, 
		sys.M ./ sys.N,
		xlabel="T", 
		ylabel="m",
	)
end

function plot_m_t(sys::Ising2D)
	initialize_plot()
	plot!(
		[ 0:1:sys.steps; ], 
		sys.M ./ sys.N,
		xlabel="iterations", 
		ylabel="m",
	)
end

# DEMOS
################################################################################

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
