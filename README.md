# PySCF: The Python-based Simulations of Chemistry Framework

## Overview

PySCF is an open-source electronic structure package for computational chemistry and materials science, written primarily in Python with performance-critical kernels in C. It is distributed under the Apache 2.0 license and is developed mainly by the Garnet Chan group and a broad community of contributors. Its original design goal, stated in the 2017 announcement paper, was to emphasize code simplicity, so that new method development is easy, and flexibility, so that computational workflows can be assembled as ordinary Python programs. The project's authors argue that the combined Python/C design achieves performance comparable to established C and Fortran quantum chemistry codes while remaining far easier to read, modify, and script.

The framework is best understood as a library rather than a monolithic executable. Calculations are expressed as Python objects and function calls, which means that a PySCF run can be interleaved with NumPy, SciPy, machine learning code, plotting tools, or custom analysis scripts in the same program. This is the feature that most distinguishes PySCF from traditional packages such as Gaussian, ORCA, or Q-Chem, which are driven by text input files.

## Scope of Methods

PySCF covers a wide range of wavefunction and density-functional methods for molecules and periodic solids.

**Mean-field methods.** Restricted, unrestricted, restricted open-shell, and generalized Hartree-Fock are all available through the `scf` module. Hartree-Fock solves the self-consistent field equations

$$
F(\mathbf{C})\,\mathbf{C} = \mathbf{S}\,\mathbf{C}\,\boldsymbol{\varepsilon}
$$

where $F$ is the Fock matrix, $\mathbf{S}$ the overlap matrix, $\mathbf{C}$ the molecular orbital coefficients, and $\boldsymbol{\varepsilon}$ the orbital energies.

**Density functional theory.** The `dft` module implements Kohn-Sham calculations with a large library of exchange-correlation functionals, including local, generalized gradient, meta-GGA, hybrid, range-separated, and double-hybrid families. PySCF does not implement functionals itself; it delegates evaluation to Libxc, which should be cited whenever functionals from it are used. Dispersion corrections, such as D3, are available through extension packages.

**Post-Hartree-Fock correlation.** The package provides second-order Møller-Plesset perturbation theory (MP2) and higher-order variants, coupled cluster methods including CCSD and CCSD(T) in the `cc` module, along with equation-of-motion coupled cluster for excited states and ionization energies. Configuration interaction methods, including CISD and full configuration interaction (FCI) for small active spaces, are in the `ci` module.

**Multiconfigurational methods.** The `mcscf` module supports complete active space SCF (CASSCF) and complete active space CI (CASCI), which are essential for systems with strong static correlation such as bond breaking, transition metal complexes, and open-shell singlets.

**Excited states and response properties.** Time-dependent Hartree-Fock and time-dependent DFT are implemented in the `tdscf` module, enabling calculation of electronic excitation spectra, oscillator strengths, and related quantities. Linear response properties and polarizabilities are also accessible.

**Analytic derivatives.** Gradients are available for a substantial set of methods, which enables geometry optimization and molecular dynamics. Second derivatives (Hessians) are available for a more limited set, which supports vibrational frequency analysis.

**Relativistic effects.** Scalar relativistic corrections are supported through exact two-component (X2C) style Hamiltonians, which matter for heavy elements where spin-orbit and scalar relativity significantly alter chemistry.

**Solvation.** Implicit solvent models, including polarizable continuum model (PCM) variants and conductor-like screening approaches, can be combined with most SCF and DFT workflows.

**Periodic systems.** The `pbc` module extends many of these methods to crystalline solids using periodic boundary conditions. It supports k-point sampling, Gaussian and plane-wave style density fitting, and Goedecker-Teter-Hutter (GTH) pseudopotentials, making it applicable to band structures and solid-state energetics.

**Localization and analysis.** Orbital localization schemes (such as Boys, Pipek-Mezey, and intrinsic bond orbital approaches), population analysis, density analysis, and properties modules let users interrogate wavefunctions after a calculation.

## Architecture and Implementation

PySCF's design follows a simple principle: almost everything is written in Python, and only the hot loops are delegated to compiled code. Molecular integrals are computed through the libcint library for Gaussian-type orbitals, and dense linear algebra relies on NumPy, which in turn links to optimized BLAS and LAPACK implementations such as OpenBLAS or Intel MKL. Heavy kernels for things like integral transformations and certain correlation terms are written in C within the package's `lib` directory.

The package is organized into modules that mirror the physics. Key ones include:

- `gto`: molecules, basis set handling, and one- and two-electron integrals
- `scf`: Hartree-Fock and related self-consistent field solvers
- `dft`: Kohn-Sham functionals and grid handling
- `mp`, `cc`, `ci`: correlation methods
- `mcscf`: multiconfigurational wavefunctions
- `tdscf`: excited states and response
- `grad` and `hessian`: analytic derivatives
- `df`: density fitting, which reduces the cost of two-electron integrals
- `ao2mo`: integral transformation to molecular orbital bases
- `solvent`: implicit solvation
- `pbc`: periodic systems
- `lo`: orbital localization
- `x2c`: relativistic corrections
- `lib`: low-level utilities and compiled extensions

Because each module is an ordinary Python package, a developer can subclass a solver, override a single step such as the Fock build or the convergence accelerator, and reuse everything else. This is a central reason the code is popular in method development.

Density fitting, which approximates the four-index electron repulsion integrals with a smaller auxiliary basis, is available throughout and is the default approach for many large-system calculations. It trades a controlled accuracy loss for substantial reductions in memory and computational cost.

## Typical Usage

A minimal calculation of the hydrogen molecule with Hartree-Fock looks like this:

```python
from pyscf import gto, scf

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='cc-pvdz')
mf = scf.RHF(mol)
energy = mf.kernel()
print("RHF energy (Hartree):", energy)
```

Switching to Kohn-Sham DFT requires changing only the solver class:

```python
from pyscf import gto, dft

mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='def2-tzvp')
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()
```

Post-Hartree-Fock methods follow the same pattern, consuming the converged reference:

```python
from pyscf import gto, scf, cc

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='cc-pvdz')
mf = scf.RHF(mol).run()
mycc = cc.CCSD(mf).run()
print("CCSD correlation energy:", mycc.e_corr)
print("CCSD total energy:", mycc.e_tot)
```

Because the results are plain Python objects and NumPy arrays, one can extract the density matrix, orbital coefficients, or energies directly and pass them to downstream code, for example to build embedding calculations, train machine learning models on electronic structure data, or post-process orbitals for visualization.

## Extensions and Ecosystem

The core package is intentionally kept focused, and additional capabilities live in separately maintained extensions installable through pip. These include modules for DMRG-based active space solvers, full configuration interaction quantum Monte Carlo, selected CI and shift-correlated approaches, dispersion corrections, semiempirical methods, and geometry optimization interfaces. The core repository offers an extras syntax, `pip install pyscf[all]`, that pulls in the bundled extensions, and individual ones can be installed on their own. Newer or experimental features are also distributed through a companion package called pyscf-forge, which allows the main code to remain stable while development proceeds faster.

PySCF's development is supported by a large community, and its design has influenced a number of research tools built on top of it, including projects for quantum embedding, quantum computing algorithm prototyping, and workflow automation. For GPU acceleration, a separate project (gpu4pyscf) provides GPU-enabled implementations of several core kernels, but this is not part of the base package.

## Installation

The most common route is pip:

```
pip install pyscf
```

Conda users can install from the community channel with `conda install pyscf::pyscf`. As of the most recent listings I checked, the Anaconda build is available for linux-64 only, and the conda package listing states that macOS support ended with version 1.6.3. Users on Windows typically run PySCF inside the Windows Subsystem for Linux. Source builds require a C compiler and BLAS/LAPACK libraries, and the installation manual covers options for threading, MPI, and custom linear algebra backends. Version numbers have progressed steadily through the 2.x series; recent listings show releases in the 2.8 to 2.13 range, so users should check the current release on PyPI or GitHub before relying on a specific version.

## Parallelism

PySCF supports shared-memory threading through OpenMP and BLAS threading, and distributed-memory parallelism through MPI via the mpi4py interface. Parallel calculations are useful for large integral transformations and correlated methods. Scaling depends heavily on the method, the system size, and the linear algebra backend, so benchmarking on the target hardware is advisable before launching large production jobs.

## Strengths

The most frequently cited strengths of PySCF are its transparency, flexibility, and integration with the Python scientific stack. Its source code is readable enough that students and researchers can trace the implementation of a method from equations to code. Method developers can prototype new ideas quickly without writing a new program from scratch. Workflow automation is natural, since parameter scans, high-throughput screening, and custom post-processing are ordinary Python loops. The package also covers an unusually wide span of methods, from simple mean-field theory through multireference and periodic calculations, within one consistent interface.

## Limitations and Trade-offs

PySCF is not designed to be the fastest code for every task. Highly optimized commercial or specialized packages can outperform it on specific large-scale calculations, particularly those that rely heavily on sophisticated parallel algorithms or GPU acceleration. Its feature set for some properties and advanced methods is narrower than that of the largest quantum chemistry suites, and certain analytic derivatives are available only for a subset of methods. Users working with very large molecules or periodic cells need to understand memory requirements and choose approximations such as density fitting carefully. Because it is a library, PySCF also demands more programming familiarity than input-file-driven codes; a new user must be comfortable writing scripts to get the most out of it.

## Comparison with Other Packages

Compared with Psi4, another Python-driven quantum chemistry package, PySCF places a stronger emphasis on the Python library model, with calculations treated as objects that can be manipulated directly, whereas Psi4 uses a more input-file-centric scripting style. Compared with Gaussian, ORCA, and Q-Chem, PySCF is open source, freely modifiable, and designed for programmatic use, while those commercial or proprietary codes typically offer more polished user interfaces, broader property coverage, and extensive production-level optimization. Compared with NWChem or Molpro, PySCF is generally more accessible to users who want to script and extend methods, though the older codes have deeper histories in certain niches.

## Documentation and Citation

Documentation, tutorials, and examples are maintained alongside the code, and the GitHub repository hosts the issue tracker for bug reports and feature requests. Users publishing results should cite the 2020 paper "Recent developments in the PySCF program package" (Sun et al., *J. Chem. Phys.* 153, 024109, 2020), and the original 2018 WIREs paper for the design overview where appropriate. If density functionals are used, Libxc should be cited as well, since PySCF relies on it for functional evaluation.

## Summary

PySCF is a flexible, open-source framework that brings quantum chemistry into the Python ecosystem. Its combination of readable Python code, performance-critical C kernels, and a broad method portfolio makes it well suited for method development, teaching, custom workflows, and research that blends electronic structure with data science. It is less suited to users who want a turnkey graphical application or the absolute highest throughput on massive production runs, but for scientists who want to understand and extend their calculations, it is one of the most accessible options available.
