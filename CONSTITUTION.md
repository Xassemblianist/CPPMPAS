# CPPMPAS constitution

1. **Goal: bit identity.** CPPMPAS reproduces MPAS-Atmosphere v8.4.2, built and run as the pinned reference recipe,
   bit for bit: tolerance 0 in every compared field and output time, on the CPU and on the GPU. A difference of one
   bit is a bug until it is found and fixed.
2. **Proof: comparison against the reference.** A capability counts only when its output matches the reference
   runs at tolerance 0. When it does not, the first differing step and field are found by bisection and the cause
   by dumps of internal arrays in both codes.
3. **Port: from the Fortran, statement by statement.** Each C++ routine names the Fortran it reproduces and keeps
   its operation order, parenthesisation, branches and loop order as the reference executable computes them.
