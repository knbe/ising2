using Plots
using UnicodePlots

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

