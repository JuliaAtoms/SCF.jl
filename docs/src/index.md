# SCF.jl

A Julia library for the self-consistent field method, mainly intended
to solve problems on the form

$$\begin{equation}
E(\vec{P},\vec{c}) = \frac{\vec{c}^H
\mat{H}\vec{c}}{\vec{c}^H\vec{c}}
\end{equation}$$

where $\vec{P}$ is a set of orbitals, $\vec{c}$ is a vector of
mixing coefficients, and

$$\begin{equation}
\mat{H}_{ij} \defd \matrixel{P_i}{\Hamiltonian}{P_j}.
\end{equation}$$

Such equations arise in the _multi-configurational
(Dirac–)Hartree–Fock_ approximation, used in theoretical atomic and
molecular physics.

This library only implements the SCF procedure. The actual equations
that are to be solved have to be setup by the user (for atoms, the
[Atoms.jl](https://github.com/JuliaAtoms/Atoms.jl) library may be
used).

## References

- Fischer, C. F. (1986). Self-Consistent-Field (SCF) and
  Multiconfiguration (MC) Hartree-Fock (HF) Methods in Atomic
  Calculations: Numerical Integration Approaches. Computer Physics
  Reports, **3**(5),
  274–325. DOI: [10.1016/0167-7977(86)90001-8](http://dx.doi.org/10.1016/0167-7977(86)90001-8)
