### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 8f14f758-f644-11ed-18f7-83704785c15e
md"""# Atomistic modelling ecosystem (4/4)

The previous notebooks have given *a few examples* for what is currently possible with Julia-based tools for atomistic modelling in materials science. Many more is out there, for which I would like to provide you with some helpful pointers.

## Disclaimer

Using Julia for atomistic modelling is still a new thing, so new projects appear on a continuous basis to explore. Not everything survives the test of time or is even intended to be a long-running project. My bias here is towards projects with considerable recent activity or which are embedded into a long-term research activity.

On top the ecosystem has grown considerably since 2018 (when I first started with Julia), so frankly I am just loosing track of whatever is out there.

Therefore if you think I have forgotten some important packages here (all apologies!). So please, if you thing something should be added here, please open a PR.
"""

# ╔═╡ 76543c65-e1cd-4e49-877c-af7394987855
md"""
## Interfacing and community building

The Julia programming language has from the start put a strong focus on the community. Reusing and sharing code in a way that it is useful for others sits at the heard of many efforts, following the figure of speach:

> **where software composes, people compose**

As a result it comes at no surprise that a number of efforts are devoted to building a community and sharing data (e.g. structures or interatomic potential definitions) between many codes. An overview.

- [JuliaMolSim](https://github.com/JuliaMolSim): The main community cluster for all Julia-related activities in atomistic simulations. Contains pretty much all major packages discussed here:
  * [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl): The most widely understood format for representing atomistic structures in Julia, IO based on [AtomsIO](https://github.com/mfherbst/AtomsIO.jl) (thus chemfiles and ASE), visualisation, etc.
  * [DFTK](https://dftk.org): Density-functional theory code
  * [Molly](https://github.com/JuliaMolSim/Molly.jl): Molecular dynamics
  * [JuLIP](https://github.com/JuliaMolSim/JuLIP.jl): ASE-like tool for manipulating structures, compatible with AtomsBase



- [MIT CESMIX](https://github.com/cesmix-mit): Overview of Julia codes written as part of [MIT's CESMIX initiative](cesmix.mit.edu), where Julia is part of the core technolgies employed for building a multiphysics simulation pipeline from DFT all the way to MD, Kinetic Monte Carlo, Flow simulations. Some efforts include:
  * [InteratomicPotentials](https://github.com/cesmix-mit/InteratomicPotentials.jl): Abstract interface for interatomic potentials to facilitate training and use of empirical/machine-learned potentials across many packages
  * [LAMMPS](https://github.com/cesmix-mit/LAMMPS.jl): Julia wrapper for the LAMMPS simulation package
  * [PotentialLearning](https://github.com/cesmix-mit/PotentialLearning.jl): Tools for training interatomic potentials



- <https://github.com/SimonEnsemble/Xtals.jl>: Tool for importing and manipulating atomic structures (similar to JuLIP and the AtomsBase ecosystem) (partial integration with AtomsBase)
"""

# ╔═╡ 6a04f2b5-c9ca-4c09-939b-aa1c2e697f2f
md"""
## Pure-Julia simulation tools

Some pure-Julia simulation tools for Materials modelling and Quantum Chemistry.

- <https://github.com/FermiQC/Fermi.jl>: Julia-based quantum chemistry code. Has a reasonalby fast coulpled-cluster implementation
- <https://github.com/NQCD/NQCDynamics.jl>: Highly flexible package for quantum dynamics (CLassical MD, Surface-hopping, ring-polymer methods, Ehrenfest)

- <https://github.com/SimonEnsemble/PorousMaterials.jl>: Classical MD for adsorption studies in metal-organic frameworks (MOFs)

- <https://github.com/f-fathurrahman/PWDFT.jl>: Pure-Julia density-functional theory package. Slightly different goals and features than [DFTK](https://dftk.org)

"""

# ╔═╡ bdca3f1f-f2cd-4c12-8ad1-0d05d1b7e4bb
md"""
## Machine learning

- The [ACEsuit](https://github.com/ACEsuit/), which has been discussed in a previous notebook.

- The [Chemellia](https://github.com/Chemellia) suite of packages features tools for featurising crystalline materials, building atomic graphs and atomic graph neural networks, e.g.
  * <https://github.com/Chemellia/AtomicGraphNets.jl>: Crystal Graph Convolutional Neural Network package
  * <https://github.com/Chemellia/ChemistryFeaturization.jl>: Featurisation toolkit, works with AtomsBase structures.


"""

# ╔═╡ e99618f2-fd8b-4c4b-b738-64b9d11c83f9
md"""
## Workflow managers & parsers

Packages for managing high-throughput workflows (e.g using Quantum Espresso), as well as parsing input / output of standard codes.

- The [MineralsCloud](https://github.com/MineralsCloud) suite of workflow managers, including:
  * <https://github.com/MineralsCloud/Express.jl>
  * <https://github.com/MineralsCloud/QuantumESPRESSO.jl>



- <https://github.com/cesmix-mit/LAMMPS.jl>: Julia wrapper for the LAMMPS simulation package



- A number of helper packages for running calculations or managing cluster jobs (here the ones by [Louis Ponet](https://github.com/louisponet)):
  * <https://github.com/louisponet/RomeoDFT.jl> 
  * <https://github.com/louisponet/RemoteHPC.jl> 
  * <https://github.com/louisponet/DFControl.jl> 
"""

# ╔═╡ a419c57f-8230-4f56-b008-e7bbc20cc81e
md"""
## Support libraries

Some libraries that provide useful primitives for atomistic and solid-state modelling. For packages dealing with parsing / manipulating structures, see the first section above.

- <https://github.com/thchr/Brillouin.jl>: Tools for determining and visualising k-space paths (similar to SeeK-path python package

- <https://github.com/numericalEFT/BrillouinZoneMeshes.jl>: Tools for working with Brillouin zone integration / meshes etc.

- <https://github.com/azadoks/PseudoPotentialIO.jl>: Parsers and processing of standard pseudopotential data formats (UPF, psp8, ...)

-  <https://github.com/JuliaMolSim/PseudoLibrary>: Pseudopotentials as Julia artifacts (which can be downloaded automatically and transparently)

- <https://github.com/singularitti/Spglib.jl>: Julia wrapper around the spglib library for crystal symmetry determination




"""

# ╔═╡ Cell order:
# ╟─8f14f758-f644-11ed-18f7-83704785c15e
# ╟─76543c65-e1cd-4e49-877c-af7394987855
# ╟─6a04f2b5-c9ca-4c09-939b-aa1c2e697f2f
# ╟─bdca3f1f-f2cd-4c12-8ad1-0d05d1b7e4bb
# ╟─e99618f2-fd8b-4c4b-b738-64b9d11c83f9
# ╟─a419c57f-8230-4f56-b008-e7bbc20cc81e
