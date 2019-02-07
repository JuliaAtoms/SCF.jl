# Quantum Systems

A _quantum system_ is here taken to be a collection of
[single-particle
orbitals](https://en.wikipedia.org/wiki/Atomic_orbital) $\vec{P}$,
arranged into multiple
[configurations](https://en.wikipedia.org/wiki/Electron_configuration),
and a set of mixing coefficients $\vec{c}$. As an example, the helium
ground state 1s² may be approximated a linear combination of [Slater
determinants](https://en.wikipedia.org/wiki/Slater_determinant):

$$\begin{equation}
\Psi(\textrm{1s²}) \approx
\sum_i c_i \Phi(\gamma_i),
\end{equation}$$

where $\gamma_i$ denotes a configuration of single-electron orbitals
and $c_i$ its associated mixing coefficient. A low-order approximation
may be achieved with the three Slater determinants formed from the 1s
and 2s orbitals:

$$\begin{equation}
\Phi(\textrm{1s²}), \quad
\Phi(\textrm{1s 2s}), \quad
\Phi(\textrm{2s²}).
\end{equation}$$

Similar ideas can be employed for molecules, etc.

```@meta
CurrentModule = SCF
```

```@docs
AbstractQuantumSystem
coefficients(::AbstractQuantumSystem)
orbitals(::AbstractQuantumSystem)
diff(::AbstractQuantumSystem)
normalize!(::AbstractQuantumSystem)
```
