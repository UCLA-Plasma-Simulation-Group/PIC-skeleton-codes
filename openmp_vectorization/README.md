# PIC Skeleton Codes:  OpenMP/Vectorization

These codes illustrate how to use hybrid shared memory/vectorization algorithm, with a tiled scheme on each shared memory multi-core node implemented wit OpenMP and vectorization implemented with both SSE vector intrinsics and compiler vectorization. The tiling scheme is described in detail in Ref.[4]. The Intel SSE2 vector intrinsics are a low level data parallel language closely related to the native assembly instructions. The compiler vectorization uses compiler directives and often requires reorganization of the data structures and loops.

For the 2D electrostatic with 12 processing cores:
no-vec = 2.7 nsec/particle/timestep
compiler vec = 2.0 nsec/particle/timestep
SSE2 = 1.6 nsec/particle/timestep

For the 2-1/2D electromagnetic with 12 processing cores:
no-vec = 9.2 nsec/particle/timestep
compiler vec = 6.1 nsec/particle/timestep
SSE2 = 4.2 nsec/particle/timestep

With SSE2 intrinsics one typically obtains about 2x speedup compared to no vectorization. Compiler vectorization achieves about 1.5x speedup.

For the 3D code, the SSE2 vector intrinsics have not been implemented.
 

1. 2D Parallel Electrostatic Spectral code:  vmpic2
2. 3D Parallel Electrostatic Spectral code:  vmpic3
3. 2-1/2D Parallel Electromagnetic Spectral code:  vmbpic2

### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


