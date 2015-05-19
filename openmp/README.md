# PIC Skeleton Codes:  OpenMP

These codes illustrate how to use shared memory parallelism with OpenMP. The algorithm used here is the same as that used on the GPU: it divides particles into small tiles and reorders them every time step to avoid data collisions when performing the deposit. Each tile is controlled by a single thread. The algorithm is described in detail in Ref. [4]. For the 2D electrostatic code, a typical execution time for the particle part of this code is about 3 ns/particle/time-step with 12 processing cores. For the 2-1/2D electromagnetic code, a typical execution time for the particle part of this code is about 10 ns/particle/time-step with 12 processing cores. For the 2-1/2D Darwin code, a typical execution time for the particle part of this code is about 22 ns/particle/time-step with 12 processing cores. (All timings are on a 2.67GHz Intel Nehalem processor.)

###Electrostatic:
1. 1D Parallel Electrostatic Spectral code:  mpic1
2. 2D Parallel Electrostatic Spectral code:  mpic2
3. 3D Parallel Electrostatic Spectral code:  mpic3

###Electromagnetic:
4. 1-2/2D Parallel Electromagnetic Spectral code:  mbpic1
5. 2-1/2D Parallel Electromagnetic Spectral code:  mbpic2
6. 3D Parallel Electromagnetic Spectral code:  mbpic3

###Darwin:
7. 1-2/2D Parallel Darwin Spectral code:  mdpic1
8. 2-1/2D Parallel Darwin Spectral code:  mdpic2
9. 3D Parallel Darwin Spectral code:  mdpic3

### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


