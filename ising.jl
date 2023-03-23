#@views point(n, lattice) = lattice[2n-1:2n]
@views lpoint(n, lattice) = lattice[n,:]

function clear_screen()
	print("\033c")
end


function show_spins(S::Vector{Int64}, L::Int64)
	n = 1
	for j = 1:L
		for i = 1:L
			#c = ""
			#S[n] < 0 ? c = "░░" : c = "██"
			S[n] < 0 ? print("░░") : print("██")
			#@printf "%s" c
			#print(c)
			n += 1
		end
		#@printf "\n"
		print("\n")
	end
	#@printf "\n"
	print("\n")

end

function print_params(T, J, N)
	L = Int64(sqrt(N))
	println("T=$(round(T, digits=4)), J=$J, L=$L, N=$N")
end

function initialise!(S::Vector{Int64}, N::Int64, lattice::Matrix{Float64})
	L = Int(floor(sqrt(N)))
	for i = 1:N
		if rand() < 0.5
			S[i] = 1
		else
			S[i] = -1
		end
	end

	#println(S)

	n = 1
	for j = 1:L
		for i = 1:L
			lattice[n,:] = [float(i), float(j)]
			n += 1
		end
	end

	show_spins(S, L)

end

function get_ssflip_dE(n::Int64, s::Vector{Int64}, lattice::Vector{Float64})

end

function get_nearest_neighbours(N::Int64, L::Int64, S, lattice::Matrix{Float64}, n::Int64)
	ij = lattice[n,:]

	(ij[2]-1) < 1 ? n_up = n + L*(L-1) : n_up = n - L
	(ij[2]+1) > L ? n_dn = n - L*(L-1) : n_dn = n + L
	(ij[1]-1) < 1 ? n_lf = n + (L-1) : n_lf = n - 1
	(ij[1]+1) > L ? n_rt = n - (L-1) : n_rt = n + 1
	

	clear_screen()
	println("showing particle:")
	print("\033[0;31m")
	println("\tn=$n, coords=$(lattice[n,:])")
	print("\033[0;37m")
	println("and nearest neighbours:")
	print("\033[1;31m")
	println("\tn_up=$(n_up), coords=$(lattice[n_up,:])")
	println("\tn_dn=$(n_dn), coords=$(lattice[n_dn,:])")
	println("\tn_lf=$(n_lf), coords=$(lattice[n_lf,:])")
	println("\tn_rt=$(n_rt), coords=$(lattice[n_rt,:])")

	print("\033[0;37m")
	num = 1
	for j = 1:L
		for i = 1:L
			#c = ""
			#S[n] < 0 ? c = "░░" : c = "██"

			if num == n
				print("\033[0;31m")
			elseif num == n_up || num == n_dn || num == n_lf || num == n_rt
				print("\033[1;31m")
			end
#			elseif num == n_up
#				print("\033[0;35m")
#			elseif num == n_down
#				print("\033[0;36m")
#			elseif num == n_left
#				print("\033[1;33m")
#			elseif num == n_right
#				print("\033[1;32m")
#			end
			S[num] < 0 ? print("▓▓") : print("██")
			#@printf "%s" c
			#print(c)
			num += 1
			print("\033[0;37m")
		end
		#@printf "\n"
		print("\n")
	end
	#@printf "\n"
	print("\n")
end

function show_nearest_neighbours(num::Int64)
	N::Int64 = 400
	L::Int64 = 20
	S::Vector{Int64} = ones(Int64, N)
	lattice = zeros(N,2)
	initialise!(S, N, lattice)

	#show_spins(S,L)
	get_nearest_neighbours(N, L, S, lattice, num)
end

function sweep()
	T = 5.0
	Tstep = 0.01
	numSteps = Int64(floor(T/Tstep))
	J = 1

	N::Int64 = 1600
	S::Vector{Int64} = ones(Int64, N)
	lattice = zeros(N,2)

	numMCIterations::Int64 = 1000000

	clear_screen()
	print_params(T, J, N)
	initialise!(S, N, lattice)
	sleep(5)

	#while T > 0.0
		#for iteration in 1:numMCIterations
		for iteration in 1:1
			n = rand(1:N)

			dE = get_ssflip_dE(n, S, lattice)
		end
	#end
end
