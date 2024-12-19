# C5G7-2D MOX benchmark
A C++ solution of the C5G7-2D MOX benchmark.
- Discrete-Ordinates Method (SN method) is used to discrete the angle term.
- Source iteration method is used to solve the k-eigenvalue problem.
- The convergence rate is accelerated by Anderson Acceleration (AA) and diffusion Synthetic Acceleration (DSA) method.
- The result is compared with which of AZTRAN code.

Compilation method:
- The code depends on the [adept library](https://github.com/rjhogan/Adept-2).
- <code>g++ C5G72D_sn.cpp -o C5G72D_sn -O3 -I "/path/to/adept"</code>

Input format (in <code>input.txt</code>):
```
nx ny N
/*
The 1.26cm * 1.26cm square area is divided into cx * cy grids; cx must be equal to cy.
N in the SN method, 2, 4, 6, 8, 12, 16 is available.
*/
l_1
l_2
...
l_nx  // Length of mesh. l_1 + l_2 + ... + l_nx = 1.26.
D // Size of coarse mesh. Positive integer, typical value is 1, 2, 3, 4, 5.
```

l_1, l_2, ..., l_nx can be generated by <code>circle_fitting.cpp</code>.
Input a positive integer not exceeding 8 and a real number c, and you can get nx and a set of l, each of which does not exceed 2c.
