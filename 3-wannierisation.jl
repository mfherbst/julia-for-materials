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
	using PlutoPlotly
	using Plots
	using SpecialFunctions
	using Unitful
	using UnitfulAtomic
	using Wannier
	import WannierPlots
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
Now let's use Wannier functions to interpolate band structures, and compare with DFTK bands.

We first construct a `Wannier.InterpModel` specifically for Wannier interpolations. This lets us interpolate band structures using `Wannier.interpolate`, which returns a `Brillouin.KPathInterpolant` storing the kpoint coordinates on a kpath, and a matrix that contains eigenvalues.
"""

# ╔═╡ 29a426dd-34dd-4463-943d-1f01554cdc6e
begin
	interp_model = Wannier.InterpModel(wann_model)
	kinter, bands_wann = Wannier.interpolate(interp_model)
end

# ╔═╡ aef8f91e-ff76-4420-8103-03b1a59f8403
md"""
The `Brillouin.KPathInterpolant` is also understood by DFTK, which we can leverage to compute the exact DFT band energies.
"""

# ╔═╡ 704fa624-ffd8-431b-a6cb-e2f41fce0e9b
dftk_bands = DFTK.compute_bands(scfres.basis, kinter; scfres.ρ);

# ╔═╡ a5250227-de3a-4f80-b114-45633107dbe5
md"""
Since DFTK uses atomic units, but `Wannier.interpolate` is hard-coded to electron volts, we also convert the DFTK data to `eV` and finally compare the two band structures.
"""

# ╔═╡ e02c8fc7-fb01-460d-a3b8-23c2045f6666
begin
	to_eV(x)     = ustrip(auconvert(u"eV", x))
	bands_dftk   = to_eV.(hcat(dftk_bands.λ...))
	fermi_energy = to_eV(scfres.εF)
end;

# ╔═╡ d2ba1b6b-e51c-4779-8dae-c6e390f2adfa
# wrap plot function so it works in Pluto notebook
PlutoPlot(
	WannierPlots.plot_band_diff(kinter, bands_dftk, bands_wann; fermi_energy).plot
)

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
# ╠═29a426dd-34dd-4463-943d-1f01554cdc6e
# ╟─aef8f91e-ff76-4420-8103-03b1a59f8403
# ╠═704fa624-ffd8-431b-a6cb-e2f41fce0e9b
# ╟─a5250227-de3a-4f80-b114-45633107dbe5
# ╠═e02c8fc7-fb01-460d-a3b8-23c2045f6666
# ╠═d2ba1b6b-e51c-4779-8dae-c6e390f2adfa
