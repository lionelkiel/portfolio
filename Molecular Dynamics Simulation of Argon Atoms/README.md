# Extracts from [Report](Report.pdf)
## Introduction and Method
The behaviour of two interacting particles is well understood. Classically, the way to solve such problems for arbitrary potential has been known for centuries. However, when we consider three particles, things become very complicated very quickly: no closed-form solution exists for this problem, even classically. This is problematic, since at our human scale, all objects are made up of incredible amounts of interacting particles. This means we cannot hope to possibly find exact solutions for the behaviour of a gas of interacting particles, for instance, and we must resort to using statistical physics.

However, all hope is not yet lost: where calculus fails, numerical analysis succeeds. By using a computer to numerically integrate the equations of motion for many-body problems, we can approximately calculate the behaviour of such systems, theoretically to arbitrary precision. We will explore this method by using it to study the dynamics of a system of argon atoms in a box with periodic boundary conditions, with the end goal of finding the pair correlation function (a measure of how frequently certain particle-particle distances appear) and pressure of the system for the gas, liquid and solid phases of argon.

The classic numerical integration method is the Euler method. However, this has one big flaw. It does not conserve the energy. Instead we thus use the velocity-verlet algorithm, which does make sure that energy is conserved. We also want the system to be at equilibrium when we perform measurements. For this we use what we called the 'lambda algorithm', where we force the system into equilibrium by scaling the velocities.

## Conclusion
It can be concluded that the simulation works well. Firstly, it conserves the total energy very well. Secondly, the experimental kinetic energy is very close to the energy in equilibrium that is to be expected. Thirdly, there is a definite increase of the organization of the atom distances as you go from gas so solid, as is to be expected.  Lastly, the values for the pressures for each combination of temperature and density are of the same order of magnitude as the pressure values one would actually measure.

# Run notebook ([MDSArgon](MDSArgon.ipynb))
- The jupyter notebook should be executed from top to bottom.
- The first cell shows all the dependencies for the notebook and for [MDSlib](MDSlib.py), which contains all the functions we made for the simulation.
- The second cell consist of the function used to simulate with the given conditions with all the sub functions it uses
- Cell 2-5 contain the analysis of a gas state
- Cell 6-9 contain the analysis of a liquid state
- Cell 10-13 contain the analysis of a solid state