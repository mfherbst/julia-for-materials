### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2bc86db6-f71a-11ed-1da5-fde04331397a
begin
	# The code in this notebook depend on unregistered packages,
	# so we disable automatic package management in Pluto and use a custom environment
	
	using Pkg
	Pkg.activate(joinpath(@__DIR__, "wannier_env"))
	Pkg.instantiate()

	using AtomsIO
	using LazyArtifacts
	using DFTK
	using PlutoUI
	using Plots
	using SpecialFunctions
	using Unitful
	using UnitfulAtomic
	using Wannier
end

# ╔═╡ b0e511c4-1fbc-42b5-9a25-030a28fe022f
LocalResource("Logos.png")

# ╔═╡ 3ec2ef1b-5edc-4510-a529-7ffb52dd8c8b
md"""
# Atomistic modelling ecosystem (2/4)

As an addendum to the previous DFTK notebook, we show here how [**DFTK**](https://dftk.org) plays nicely with the [**Wannier**](https://github.com/qiaojunfeng/Wannier.jl) Julia package for wannierisation, which has recently appeared.

!!! warning "Unstable interfaces"
    This file shows a current working version of the DFTK <-> Wannier interface.
    Names of functions and the overall API is preliminary and will likely change.
    The [DFTK documentation](https://docs.dftk.org) and the [Wannier documentation](https://wannierjl.org) will feature a well-explained update example once the interface is ready.

We will go with a graphene system. Choose discretisation according to your CPU:
  - Ecut (Hartree):    $(@bind Ecut Slider(10:1.0:40; default=10, show_value=true))
  - nkpoints:    $(@bind nkpoints Slider(2:1:12; default=5, show_value=true))
  - nbands:     $(@bind nbands Slider(10:1:25; default=15, show_value=true))
"""

# ╔═╡ 19c9abdb-9e05-4b2e-8e07-39655cebca89
begin
	d = 10u"Å"     # Graphene sheet distance
	a = 2.641u"Å"  # Graphene lattice constant
	lattice = [a  -a/2    0;
	           0  √3*a/2  0;
	           0     0    d]

	# Pseudo-Dojo standard pseudo for carbon
	C = ElementPsp(:C, psp=load_psp(artifact"pd_nc_sr_pbe_standard_0.4.1_upf/C.upf"))
	atoms     = [C, C]
	positions = [[0.0, 0.0, 0.0], [1//3, 2//3, 0.0]]
	
	model = model_PBE(lattice, atoms, positions)
	basis = PlaneWaveBasis(model; Ecut, kgrid=(nkpoints, nkpoints, 1))
	
	nbandsalg = AdaptiveBands(basis.model; n_bands_converge=nbands)
	scfres = self_consistent_field(basis; nbandsalg, tol=1e-5)
end;

# ╔═╡ 8b7d5f5a-3435-4a58-a381-4bf06bc84eb6
plot_bandstructure(scfres)

# ╔═╡ e9cc66d3-0fc5-4b6f-8321-95b31ec706fb
md"""
Now we use the `run_wannier` routine to perform the Wannierization procedure
using Wannier.jl. We first use SCDM to generate a better initial guess for the MLWFs.

Since this is an entangled case, we need a weight factor, for which we use an `erfc`
with parameters
  - `σ` = $(@bind σ Scrubbable(1e-3 : 2e-3 : 2e-2; format=".3f", default=0.01)) Hartree
  - `μ` = $(@bind μ Scrubbable(0 : 1e-3 : 1e-2; format=".3f", default=0.0)) Hartree

and construct `nwann = ` $(@bind nwann Slider(3 : 1 : nbands-4; default=5, show_value=true)) wannier functions:
"""

# ╔═╡ 89498877-ad4c-4d80-9661-54bf2ffafc4c
A = let
	# Exclude non-converged bands
	exclude_bands = DFTK._default_exclude_bands(scfres)
	basis, ψ, eigenvalues = DFTK.unfold_scfres_wannier(scfres, exclude_bands)

	function scdm_f(kpt::Kpoint)
		ik = findfirst(Ref(kpt) .== basis.kpoints)
    	εk = eigenvalues[ik]
		@. 0.5 * erfc((εk - μ) / σ)
	end

	# Builds and returns the A matrix
	A = DFTK.compute_amn_scdm(basis, ψ, nwann, scdm_f)
end;

# ╔═╡ 5983120f-9b88-4795-b469-0752c7773efd
wann_model = only(
	run_wannier(scfres;
    	fileprefix="wannier/graphene",
    	n_wann=nwann,
    	A,
    	dis_froz_max=0.1
	)
)

# ╔═╡ c5b2f3ac-02f8-4363-a925-32c132a1a0e2
md"""
!!! danger "TODO Visualisation"
    Currently WannierPlots is still private ...
"""

# ╔═╡ Cell order:
# ╠═2bc86db6-f71a-11ed-1da5-fde04331397a
# ╠═b0e511c4-1fbc-42b5-9a25-030a28fe022f
# ╟─3ec2ef1b-5edc-4510-a529-7ffb52dd8c8b
# ╠═19c9abdb-9e05-4b2e-8e07-39655cebca89
# ╠═8b7d5f5a-3435-4a58-a381-4bf06bc84eb6
# ╟─e9cc66d3-0fc5-4b6f-8321-95b31ec706fb
# ╠═89498877-ad4c-4d80-9661-54bf2ffafc4c
# ╠═5983120f-9b88-4795-b469-0752c7773efd
# ╟─c5b2f3ac-02f8-4363-a925-32c132a1a0e2
