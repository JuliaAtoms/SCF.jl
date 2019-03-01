var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SCF.jl-1",
    "page": "Home",
    "title": "SCF.jl",
    "category": "section",
    "text": "A Julia library for the self-consistent field method, mainly intended to solve problems on the formbeginequation\nE(vecPvecc) = fracvecc^H\nmatHveccvecc^Hvecc\nendequationwhere vecP is a set of orbitals, vecc is a vector of mixing coefficients, andbeginequation\nmatH_ij defd matrixelP_iHamiltonianP_j\nendequationSuch equations arise in the multi-configurational (Dirac–)Hartree–Fock approximation, used in theoretical atomic and molecular physics.This library only implements the SCF procedure. The actual equations that are to be solved have to be setup by the user (for atoms, the Atoms.jl library may be used)."
},

{
    "location": "#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Fischer, C. F. (1986). Self-Consistent-Field (SCF) and Multiconfiguration (MC) Hartree-Fock (HF) Methods in Atomic Calculations: Numerical Integration Approaches. Computer Physics Reports, 3(5), 274–325. DOI: 10.1016/0167-7977(86)90001-8"
},

{
    "location": "quantum_systems/#",
    "page": "Quantum Systems",
    "title": "Quantum Systems",
    "category": "page",
    "text": ""
},

{
    "location": "quantum_systems/#SCF.AbstractQuantumSystem",
    "page": "Quantum Systems",
    "title": "SCF.AbstractQuantumSystem",
    "category": "type",
    "text": "AbstractQuantumSystem\n\nAny implementation of the AbstractQuantumSystem interface, used for self-consistent field calculations, must implement the functions coefficients, orbitals, diff, normalize!.\n\n\n\n\n\n"
},

{
    "location": "quantum_systems/#SCF.coefficients-Tuple{AbstractQuantumSystem}",
    "page": "Quantum Systems",
    "title": "SCF.coefficients",
    "category": "method",
    "text": "coefficients(quantum_system)\n\nRetrieves the mixing coefficients of quantum_system. Must return a view that the scf! routine can modify. To be overloaded by the implementation of AbstractQuantumSystem.\n\n\n\n\n\n"
},

{
    "location": "quantum_systems/#SCF.orbitals-Tuple{AbstractQuantumSystem}",
    "page": "Quantum Systems",
    "title": "SCF.orbitals",
    "category": "method",
    "text": "orbitals(quantum_system)\n\nRetrieves the orbitals of quantum_system. Must return a view that the scf! routine can modify. To be overloaded by the implementation of AbstractQuantumSystem.\n\n\n\n\n\n"
},

{
    "location": "quantum_systems/#Base.diff-Tuple{AbstractQuantumSystem}",
    "page": "Quantum Systems",
    "title": "Base.diff",
    "category": "method",
    "text": "diff(quantum_system[; kwargs...])\n\nVaries the the quantum_system with respect to all orbitals. Used to derive the multi-configurational Hartree–Fock equations. To be overloaded by the implementation of AbstractQuantumSystem.\n\n\n\n\n\n"
},

{
    "location": "quantum_systems/#LinearAlgebra.normalize!-Tuple{AbstractQuantumSystem}",
    "page": "Quantum Systems",
    "title": "LinearAlgebra.normalize!",
    "category": "method",
    "text": "normalize!(quantum_system, v)\n\nNormalize the orbital v of quantum_system. To be overloaded by the implementation of AbstractQuantumSystem.\n\n\n\n\n\n"
},

{
    "location": "quantum_systems/#Quantum-Systems-1",
    "page": "Quantum Systems",
    "title": "Quantum Systems",
    "category": "section",
    "text": "A quantum system is here taken to be a collection of single-particle orbitals vecP, arranged into multiple configurations, and a set of mixing coefficients vecc. As an example, the helium ground state 1s² may be approximated a linear combination of Slater determinants:beginequation\nPsi(textrm1s²) approx\nsum_i c_i Phi(gamma_i)\nendequationwhere gamma_i denotes a configuration of single-electron orbitals and c_i its associated mixing coefficient. A low-order approximation may be achieved with the three Slater determinants formed from the 1s and 2s orbitals:beginequation\nPhi(textrm1s²) quad\nPhi(textrm1s 2s) quad\nPhi(textrm2s²)\nendequationSimilar ideas can be employed for molecules, etc.CurrentModule = SCFAbstractQuantumSystem\ncoefficients(::AbstractQuantumSystem)\norbitals(::AbstractQuantumSystem)\ndiff(::AbstractQuantumSystem)\nnormalize!(::AbstractQuantumSystem)"
},

{
    "location": "fock_operators/#",
    "page": "Fock Operators",
    "title": "Fock Operators",
    "category": "page",
    "text": ""
},

{
    "location": "fock_operators/#Fock-Operators-1",
    "page": "Fock Operators",
    "title": "Fock Operators",
    "category": "section",
    "text": "CurrentModule = SCFThe Fock operator consists of an AbstractQuantumSystem and a set of coupled integro-differential equations, the solution of which is the objective of the self-consistent field procedure. For the solution process to work, the set of equations need to implement a few methods:It must fulfil Julia\'s iteration interface, i.e. each element must be the equation for a single orbital, which is solved independently from the other equations, but with the other orbitals as inputs.\nenergy which calculates the orbital energy, i.e. matrixelchiHamiltonianchi.\nenergy_matrix! which calculates the energy matrix matH_i for orbital equation i, where vecc^HmatH_ivecc gives the orbital energy for the correspond orbital. The overall energy matrix matH=sum_imatH_i is used to solve the secular problem for the mixing coefficients.\nupdate! which recomputes all orbital-dependent integrals, shared among the equations of the equation system.Fock\nnorm_rot!\nrotate_max_lobe!\nenergy\nenergy_matrix!\nhamiltonian\nupdate!"
},

{
    "location": "self_consistent_iteration/#",
    "page": "Self-Consistent Iteration",
    "title": "Self-Consistent Iteration",
    "category": "page",
    "text": ""
},

{
    "location": "self_consistent_iteration/#SCF.scf!",
    "page": "Self-Consistent Iteration",
    "title": "SCF.scf!",
    "category": "function",
    "text": "scf!([fun!, ]fock[; ω=0, max_iter=200, tol=1e-8, verbosity=1])\n\nPerform the SCF iterations for the fock operator. One iteration is performed by scf_iteration!. The optional fun!(P̃,c̃) argument allows for extra operations to be performed every SCF cycle (such as plotting, etc).\n\nω is a relaxation parameter; the orbitals and mixing coefficients are updated as wᵢ₊₁ = (1-ω)w̃ + ωwᵢ where w̃ is the solution in the current iteration and wᵢ the previous solution. The default (ω=0) is to only use this solution. A larger value of ω makes it easier to achieve convergence in oscillatory problems, but a too large value decreases the convergence rate. It is therefore desireable to reduce the value of ω if the convergence is deemed to be monotonous, the criterion for which is whether the change over the monotonous_window last iterations is steadily decreasing. If this is the case, ω will be reduced by multiplying by ωfactor. Conversely, if the change is non-monotonous, ω will be increased, but not beyond ωmax.\n\nThe SCF procedure continues until either the amount of iterations equals max_iter or the change in the coefficients is below tol.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#SCF.scf_iteration!",
    "page": "Self-Consistent Iteration",
    "title": "SCF.scf_iteration!",
    "category": "function",
    "text": "scf_iteration!(fock, P, c[; verbosity=0, method=:arnoldi, tol=1e-10,\n               update_mixing_coefficients=true])\n\nPerform one step of the SCF iteration. This will\n\nImprove each of the orbitals in turn, using the values of the orbitals P and mixing coefficients c from the previous step as input.\nSolve the secular problem to find new values for the mixing coefficients, c, unless update_mixing_coefficients==false.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#SCF.solve_secular_problem!",
    "page": "Self-Consistent Iteration",
    "title": "SCF.solve_secular_problem!",
    "category": "function",
    "text": "solve_secular_problem!(H, c, fock)\n\nForm the energy matrix, store it in H, and then solve the secular problem Hc = Ec for the lowest eigenvalue.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#Self-Consistent-Iteration-1",
    "page": "Self-Consistent Iteration",
    "title": "Self-Consistent Iteration",
    "category": "section",
    "text": "CurrentModule = SCFscf!\nscf_iteration!\nsolve_secular_problem!"
},

{
    "location": "self_consistent_iteration/#SCF.KrylovWrapper",
    "page": "Self-Consistent Iteration",
    "title": "SCF.KrylovWrapper",
    "category": "type",
    "text": "KrylovWrapper(hamiltonian)\n\nProxy object used in the Krylov iterations, during orbital improvement. This is useful, since hamiltonian may be defined to act on objects such as vectors living in function spaces (as e.g. implemented using ContinuumArrays.jl), whereas the SCF iterations act on the coefficients directly.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#Base.size-Tuple{SCF.KrylovWrapper}",
    "page": "Self-Consistent Iteration",
    "title": "Base.size",
    "category": "method",
    "text": "size(::KrylovWrapper)\n\nReturns the dimension of the KrylovWrapper. For Hamiltonians which are not <:AbstractMatrix, this needs to be overloaded.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#LinearAlgebra.mul!",
    "page": "Self-Consistent Iteration",
    "title": "LinearAlgebra.mul!",
    "category": "function",
    "text": "mul!(y, ::KrylovWrapper, x)\n\nCompute the action of the wrapped Hamiltonian on x and store it in y. For Hamiltonians which are not <:AbstractMatrix, this needs to be overloaded.\n\n\n\n\n\n"
},

{
    "location": "self_consistent_iteration/#KrylovWrapper-1",
    "page": "Self-Consistent Iteration",
    "title": "KrylovWrapper",
    "category": "section",
    "text": "The orbital improvements are performed via diagonalization of the integro-differential equations for each orbital. The diagonalization procedure is the Arnoldi method, which requires the repeated action of the Hamiltonian on a vector. The KrylovWrapper object is used to wrap Hamiltonians which are not simple <:AbstractMatrix objects, but which require a specialized implementation of mul! to act on a <:AbstractVector of coefficients. An example is the kind of Hamiltonians implemented in the Atoms.jl library, which are based on the function space algebra from ContinuumArrays.jl.KrylovWrapper\nBase.size(::KrylovWrapper)\nmul!"
},

]}
