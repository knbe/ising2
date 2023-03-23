# ising 2d.

const global J::Int64 = 1

mutable struct Lattice2D
	L::Int64
	N::Int64
	S::Vector{Int64}
	coords::Matrix{Float64}
end

@enum LatticeType begin
	square_random
	triangular_random
end

function Lattice2D(;L::Int64=40, ltype::LatticeType=square_random)
	N, S, coords = initialise(L, Int(ltype))
	return Lattice2D(L, N, S, coords)
end

@views lpoint(n,lattice) = lattice[n,:]

function clear_screen()
	print("\033c")
end

function show_spins(lattice::Lattice2D)
	n = 1
	for j = 1:lattice.L
		for i = 1:lattice.L
			lattice.S[n] < 0 ? print("░░") : print("██")
			n += 1
		end
		#@printf "\n"
		print("\n")
	end
	#@printf "\n"
	print("\n")

end

function print_params(lattice::Lattice2D, T::Float64)
	println(
		"T=$(round(T, digits=4)), \
		L=$(lattice.L), \
		N=$(lattice.N)"
	)
end

function initialise(L::Int64, ltype::Int64)
	N = L^2
	S = ones(Int64, N)
	coords = zeros(N, 2)

	# square lattice, random positions
	if ltype == 0
		for i = 1:N
			if rand() < 0.5
				S[i] = 1
			else
				S[i] = -1
			end
		end

		n = 1
		for j = 1:L
			for i = 1:L
				coords[n,:] = [float(i), float(j)]
				n += 1
			end
		end
	end

	return N, S, coords

end

# returns indices of the 4 nearest neighbours of the nth particle
function get_nearest_neighbours(l::Lattice2D, n::Int64)
	coords = lpoint(n, l.coords)

	(coords[2]-1) < 1 ? n_up = n + l.L*(l.L-1) : n_up = n - l.L
	(coords[2]+1) > l.L ? n_dn = n - l.L*(l.L-1) : n_dn = n + l.L
	(coords[1]-1) < 1 ? n_lf = n + (l.L-1) : n_lf = n - 1
	(coords[1]+1) > l.L ? n_rt = n - (l.L-1) : n_rt = n + 1

	return n_up, n_dn, n_lf, n_rt
end

function show_nearest_neighbours(n::Int64, lattice::Lattice2D)
	n_up, n_dn, n_lf, n_rt = get_nearest_neighbours(l, n)

	clear_screen()
	println("showing particle:")
	print("\033[0;31m")
	println("\tn=$n, coords=$(l.coords[n,:])")
	print("\033[0;37m")
	println("and nearest neighbours:")
	print("\033[1;31m")
	println("\tn_up=$(n_up), coords=$(l.coords[n_up,:])")
	println("\tn_dn=$(n_dn), coords=$(l.coords[n_dn,:])")
	println("\tn_lf=$(n_lf), coords=$(l.coords[n_lf,:])")
	println("\tn_rt=$(n_rt), coords=$(l.coords[n_rt,:])")

	print("\033[0;37m")
	index = 1
	for j = 1:lattice.L
		for i = 1:lattice.L
			if index == n
				print("\033[0;31m")
			elseif index == n_up || index == n_dn || 
				index == n_lf || index == n_rt
				print("\033[1;31m")
			end
			l.S[n] < 0 ? print("▓▓") : print("██")
			index += 1
			print("\033[0;37m")
		end
		print("\n")
	end
	print("\n")
end

function get_ssflip_dE(lattice::Lattice2D, n::Int64)
	n_up, n_dn, n_lf, n_rt = get_nearest_neighbours(lattice, n)

	return 2 * J * lattice.S[n] * (lattice.S[n_up] +
		lattice.S[n_dn] + lattice.S[n_lf] + lattice.S[n_rt])
end

function sweep(lattice::Lattice2D; T₀::Float64=5.0)
	T = T₀
	Tstep = 0.005

	numMCIterations = 1000000
	numTsteps = Int(T/Tstep)

	clear_screen()
	print_params(lattice, T)
	show_spins(lattice)
	sleep(5)

	count = 0
	while T > 0.0
		for iteration in 1:numMCIterations
			n = rand(1:lattice.N)

			dE = get_ssflip_dE(lattice, n)

			if dE < 0
				lattice.S[n] *= -1
			else
				if rand() < exp(-dE / T)
					lattice.S[n] *= -1
				end
			end
		end

		if count % 2 == 0
			clear_screen()
			print_params(lattice, T)
			show_spins(lattice)
		end

		T -= Tstep
		count += 1
	end
end
