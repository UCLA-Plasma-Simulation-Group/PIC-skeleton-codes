# PIC Skeleton Codes:  MPI

These codes illustrate how to use domain decomposition with message-passing (MPI). This is the dominant programming paradigm for PIC codes today. These codes use only a simple 1d domain decomposition. The algorithm is described in detail in Refs. [2-3]. For the 2D electrostatic code, a typical execution time for the particle part of this code is about 140 ps/particle/time-step with 256 MPI nodes. For the 2-1/2D electromagnetic code, a typical execution time for the particle part of this code is about 400 ps/particle/time-step with 256 MPI nodes.  For the 2-1/2D Darwin code, a typical execution time for the particle part of this code is about 460 ps/particle/time-step with 256 MPI nodes. The CPUs (2.67GHz Intel Nehalem processors) were throttled down to 1.6 GHz for these benchmarks.

1. 2D Parallel Electrostatic Spectral code:  ppic2
2. 2-1/2D Parallel Electromagnetic Spectral code:  pbpic2
2. 2-1/2D Parallel Darwin Spectral code:  pdpic2


### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


