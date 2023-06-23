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

# ╔═╡ aa1739bd-906e-492d-b3d7-4ea394fa8cdb
begin
	# The code in this notebook depend on non-standard package registries,
	# so we disable automatic package management in Pluto and use
	# a custom environment
	
	using Pkg
	Pkg.activate(joinpath(@__DIR__, "ace_env"))
	Pkg.instantiate()

	using AtomsBase
	using AtomsView
	using PlutoUI
	
	# These require the ACE registry, see https://github.com/ACEsuit/ACEregistry
	# for details
	using JuLIP
	using ACE1pack
	using ACEfit
end

# ╔═╡ f0d14eb6-e833-4158-b687-df4abb778ba2
md"""
!!! warning "Non-standard registries"
    This notebook depends on packages in the special [ACEregistry](https://github.com/ACEsuit/ACEregistry). It therefore does not work out of the box.
    To use it follow the instructions on the [ACEregistry readme](https://github.com/ACEsuit/ACEregistry).
"""

# ╔═╡ dbdfe11c-f716-11ed-1b54-73f192299d19
LocalResource("Logos.png")

# ╔═╡ ce164656-33de-4a9a-83db-faf409387576
md"""
# Atomistic modelling ecosystem (3/4)
"""

# ╔═╡ 82bab38d-ba86-4c60-8ade-4a59431a5c53
md"""
## JuLIP
"""

# ╔═╡ f4aa222e-95af-4b1b-95fb-f28928ae02f3
md"""
[**JuLIP**](https://github.com/JuliaMolSim/JuLIP.jl/) is a library for assembling atomistic structures and quickly testing novel interatomic potentials.

For example we can construct a Silicon vacancy.


We consider a $(@bind nx Scrubbable(1:1:3; default=2)) x $(@bind ny Scrubbable(1:1:3; default=2)) x $(@bind nz Scrubbable(1:1:3; default=1)) supercell.
"""

# ╔═╡ a438202e-11df-4adc-8bf3-682fe336a801
begin
	at = bulk(:Si, cubic=true) * (nx, ny, nz)
	deleteat!(at, 1)
end

# ╔═╡ 62371036-25df-420b-9d98-0386e81c1eff
md"Now we use the fact that JuLIP structures can be seamlessly **converted to AtomsBase** structures (which Pluto knows how to visualise) to see the outcome:"

# ╔═╡ dac59da6-b551-49f8-9284-bcf2579362e3
convert(AbstractSystem, at)

# ╔═╡ b223228c-0062-4c41-9fd3-f5b09ddcd0b4
md"""
While we could now easily employ DFTK to evaluate the energy of such a structure, e.g.
```julia
using DFTK

system = attach_psp(convert(AbstractSystem, at); 
                    Si=artifact"pd_nc_sr_pbe_standard_0.4.1_upf/Si.upf"))
model  = model_PBE(system)
basis  = PlaneWaveBasis(model; Ecut=20, kgrid=(1,1,1))
self_consistent_field(basis).energies.total
```
this is a bit too expensive for such a notebook. 

Instead we build our own **custom analytic potential** and use that:
"""

# ╔═╡ c0b88cfb-7646-4b20-98a1-e37cd2ae4e3f
begin
	r0 = rnn(:Si)
	pot = let A = 4.0, r0 = r0
	   @analytic r -> 6.0 * exp(- A * (r/r0 - 1.0)) - A * (r0/r)^6
	end
	pot = pot * SplineCutoff(2.1 * r0, 3.5 * r0)

	energy(pot, at)  # evaluate
end

# ╔═╡ 9ca6b1aa-44b7-426d-9751-b51fca84334b
md"Some standard potentials (EAM, StillingerWeber, Lennard-Jones) are also supported out of the box."

# ╔═╡ ce83e65a-d810-4f96-80a9-b0e6ea3bf869
md"""
## The ACE Julia ecosystem

The atomic cluster expansion (ACE) family of **machine-learning interatomic potentials**  has a well-developed Julia-based ecosystem, which can be found in the [**ACEsuit**](https://acesuit.github.io/) github organisation. Working with and using ACE models is closely linked with JuLIP.

Here we will follow closely one of the tutorials of the [**ACE1pack**](https://github.com/ACEsuit/ACE1pack.jl) package.

The goal of this tutorial is to **fit a TiAl potential** against pre-computed DFT data.

First we get the datafile (which is downloaded automatically via Julia's artifact system):
"""

# ╔═╡ 327fe020-6ddf-4110-be04-ede883c6a242
begin
	datafile = joinpath(ACE1pack.artifact("TiAl_tutorial"), "TiAl_tutorial.xyz")
	data = read_extxyz(datafile)
end

# ╔═╡ 29085db1-f047-4de7-ae06-675a5566f91e
md"We will not consider all data, but only a sparse subset (to speed up the process):"

# ╔═╡ 39799492-369a-44d0-95a8-ae28b3cd1304
datafit = data[1:5:end]

# ╔═╡ 65535a3e-b3c6-468b-9a3c-f7ae66d79f2a
length(data), length(datafit)

# ╔═╡ 06834bb0-d826-4010-8e59-ad33abe07b06
md"""
Next we generate a model. We take
  - `order = 3`, i.e. a 4-body potential
  - `totaldegree = 6`: a very low polynomial degree for testing
  - `rcut = 5.5`, typical cutoff radius for metals
  - Supply one-body reference potentials

(See the [ACE1pack tutoral](https://acesuit.github.io/ACE1pack.jl/dev/literate_tutorials/TiAl_model/) for details)
"""

# ╔═╡ 94c85509-7f81-461d-be80-710ed6b55b81
model = acemodel(; elements = [:Ti, :Al],
                   order = 3,
                   totaldegree = 6,
                   rcut = 5.5,
                   Eref = [:Ti => -1586.0195, :Al => -105.5954])

# ╔═╡ 355cf186-e872-4fc4-8ba0-ce72fa38af67
md"""
Next we fit the models, taking specific weights between the errors in the energies, forces and virials and employing a prior to direct us towards a smooth solution (regularisation):
"""

# ╔═╡ 6aafc12a-ad77-4022-9d76-5fd18698d0cb
weights = Dict(
	"FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 1.0 ),
	"TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 )
);

# ╔═╡ e3e5f994-5d57-4615-9212-e624b203bed0
md"""
!!! warning "Fitting currently broken"
    According to the documentation, this code works, but I was unable to get it to work:
	```julia
	acefit!(model, datafit;
			solver=ACEfit.LSQR(damp=1e-2, atol=1e-6),
			prior=smoothness_prior(model; p=3))
	```
    The rest of the code still runs, but since we have done zero fitting, the result is garbage
"""

# ╔═╡ ddea5c5d-1d1d-4604-adcf-97f56c6d4695
md"The ACEsuit provides various tools for testing the fit, e.g."

# ╔═╡ 5bc44af8-3dcd-44a0-9866-697fe5581987
begin
	test_data = data[2:10:end]  # Unseen data
	ACE1pack.linear_errors(test_data, model; weights)
end;

# ╔═╡ a61671c6-3e73-40cc-a02d-07c784bebd3e
md"Once the model is trained it can be directly used in combination with JuLIP to compute e.g. the energy of an unseen structure:"

# ╔═╡ 9ba4c5b6-d965-42e9-b2c8-3de0ba0c6a61
begin
	atoms = data[3]  # Extract an unseen structure
	energy(model.potential, atoms)
end

# ╔═╡ Cell order:
# ╟─f0d14eb6-e833-4158-b687-df4abb778ba2
# ╠═aa1739bd-906e-492d-b3d7-4ea394fa8cdb
# ╠═dbdfe11c-f716-11ed-1b54-73f192299d19
# ╟─ce164656-33de-4a9a-83db-faf409387576
# ╟─82bab38d-ba86-4c60-8ade-4a59431a5c53
# ╟─f4aa222e-95af-4b1b-95fb-f28928ae02f3
# ╟─a438202e-11df-4adc-8bf3-682fe336a801
# ╟─62371036-25df-420b-9d98-0386e81c1eff
# ╠═dac59da6-b551-49f8-9284-bcf2579362e3
# ╟─b223228c-0062-4c41-9fd3-f5b09ddcd0b4
# ╠═c0b88cfb-7646-4b20-98a1-e37cd2ae4e3f
# ╟─9ca6b1aa-44b7-426d-9751-b51fca84334b
# ╟─ce83e65a-d810-4f96-80a9-b0e6ea3bf869
# ╠═327fe020-6ddf-4110-be04-ede883c6a242
# ╟─29085db1-f047-4de7-ae06-675a5566f91e
# ╠═39799492-369a-44d0-95a8-ae28b3cd1304
# ╠═65535a3e-b3c6-468b-9a3c-f7ae66d79f2a
# ╟─06834bb0-d826-4010-8e59-ad33abe07b06
# ╠═94c85509-7f81-461d-be80-710ed6b55b81
# ╟─355cf186-e872-4fc4-8ba0-ce72fa38af67
# ╠═6aafc12a-ad77-4022-9d76-5fd18698d0cb
# ╟─e3e5f994-5d57-4615-9212-e624b203bed0
# ╟─ddea5c5d-1d1d-4604-adcf-97f56c6d4695
# ╠═5bc44af8-3dcd-44a0-9866-697fe5581987
# ╟─a61671c6-3e73-40cc-a02d-07c784bebd3e
# ╠═9ba4c5b6-d965-42e9-b2c8-3de0ba0c6a61
