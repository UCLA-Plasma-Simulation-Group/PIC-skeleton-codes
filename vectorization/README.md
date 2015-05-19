# PIC Skeleton Codes:  Vectorization

These codes illustrate how to use vectorization with the Intel Processors. Two approaches are illustrated. One uses the Intel SSE2 vector intrinsics, which is a low level data parallel language closely related to the native assembly instructions. This gives the best performance but requires substantial effort and expertise. The other approach uses compiler directives and often requires reorganization of the data structures and loops, but is much simpler.

For the 2D electrostatic:
  * no-vec = 35 nsec/particle/timestep
  * compiler vec = 18 nsec/particle/timestep
  * SSE2 = 12 nsec/particle/timestep

For the 2-1/2D electromagnetic:
  * no-vec = 100 nsec/particle/timestep
  * compiler vec = 60 nsec
  * SSE2 = 34 nsec

With SSE2 intrinsics one typically obtains about 3x speedup compared to no vectorization. Compiler vectorization achieves about 2x speedup. (All timings are on a 2.67GHz Intel Nehalem processor.)

 

1. 2D Parallel Electrostatic Spectral code:  vpic2
2. 2-1/2D Parallel Electromagnetic Spectral code:  vbpic2



### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


