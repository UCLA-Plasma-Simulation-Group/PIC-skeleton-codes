# Particle-In-Cell Skeleton Codes
Succinct PIC Implementation 

__All source code in this repository is confidential.__

If you're viewing this and don't have explicit permission from the UCLA Plasma Physics Simulation Group to access this repository please contact Warren Mori (mori@physics.ucla.edu).

## Quick Start

Welcome! There is a large collection of PIC codes here with a range of complexity. They are all deigned to be easy to read, understand, and adapt. 

1. Each directory is a complete, standalone PIC implemenation.

2. *Where to Start*? The directory `pic2` is the basic, classic serial Particle-In-Cell code. For those new to PIC it is the perfect place to start.. no bells or whistles.. just a simple way to see how PIC codes work. For people already experienced with PIC, `pic2` can be used in benchmarking or as unit-testing scaffolding.

3. These codes have no external dependecies - You should be able to set basic coniguration and then `make`.. For full instructions see the readme files in each directory.

##Introduction

...We are providing a hierarchy of skeleton PIC codes, from very basic serial codes for beginning students to very sophisticated codes with 3 levels of parallelism for HPC experts.  The codes are designed for high performance, but not at the expense of code obscurity.  They have multiple purposes.  The first is educational, by providing example codes for computational physics researchers of what a PIC code is and how to parallelize it with different parallel strategies and languages.  The second is to provide benchmark codes for computer science researchers to evaluate new computer hardware and software with a non-trivial computation.  Third, for researchers with already existing codes, these skeleton codes can be mined for useful ideas and procedures that can be incorporated into their own codes.  Finally, the skeleton codes can be expanded by adding diagnostics, initial conditions, and so forth, to turn them into production codes.

All the codes in the initial group are 2 dimensional, and for simplicity are designed for uniform plasmas.   They are divided into layers with increasing levels of parallelism.  There are two types of codes in each layer, simple electrostatic codes and more complex electromagnetic codes. Both codes use the same algorithms, so if one’s interest is primarily in the algorithms, the electrostatic codes are recommended.  The codes are executing the same calculation for each level of parallelism, but the execution is done differently.  Each level of parallelism incrementally increases the complexity, and enables the program to execute faster.  This hierarchy of codes could be used in a semester course on parallel PIC algorithms, so that a beginning student could become an HPC expert by the end.  Spectral field solvers based on Fast Fourier Transforms (FFTs) are used, but the main focus of these skeleton codes is on processing particles, which normally dominates the CPU time.

This work is also supported in part by the National Science Foundation, Department of Energy SciDAC program and the UCLA Institute for Digital Research and Innovation (IDRE).  It is also closely co-ordinated with the activities of the NSF funded Particle-in-Cell and Kinetic Simulation Center (PICKSC) at UCLA.  The codes are available for download on the following web site:  https://idre.ucla.edu/hpc/parallel-plasma-pic-codes.

## Software Design

Most of the codes are written in two languages, Fortran and C.  Language interoperability is also maintained, so that Fortran codes can call C libraries and vice-versa.  Fortran2003 has built-in support for interoperability with C.  C procedures can be called by simply writing a file containing interface statements (similar to a header file in C).  To allow C to call Fortran procedures, a simple wrapper library must be built in Fortran to allow interoperability with Fortran arrays (which is functionally similar to array objects in C++).  Interoperability is important for a number of reasons.  The first is to give Fortran access to features available only in C, for example, to permit a Fortran program to call a CUDA C procedure on a GPU.  The second is to allow a C program access to procedures already available in Fortran.  Our goal is to support multi-lingual codes, where researchers write in whatever language they are fluent in, and the pieces can be integrated seamlessly.   Because the features of Fortran2003 are not yet well known in the scientific community, an older style of interoperability is also currently provided.   The details of this approach are discussed in Ref. [1], but the short story is that one line wrapper functions are written in C, one function for calling from Fortran to C, and another function for calling C from Fortran.   Most of the libraries in the skeleton codes have interface functions to make them interoperable.

The computationally intensive procedures in Fortran are written in a Fortran77 style subset of Fortran90.  This is primarily because programs written in such a style tend to execute faster.  But an additional reason is that Fortran has a large legacy of codes written in this style and it is useful that students be aware of such codes.  All libraries written in such a style have a corresponding interface library to enable type checking when used with Fortran90.  The C codes adhere to the C99 standard, primarily to make use of float complex types.

Most of the codes maintain a simple two level hierarchy, with a main code and an encapsulated library.  In a some cases, additional layers are needed.  One example is a library which makes use of MPI communications.  In such a case, the MPI features are usually completely contained within the library and the MPI details (or MPI header files) are not exposed to other layers.  This kind of isolation more easily enables integrating different programming layers, such as GPU and MPI programming.

A flat data structure is used, using only basic integer and floating point types, without the use of derived types or structs.  This is primarily to expose as much as possible without hiding features in a multi-layered hierarchy.  We believe this makes it easier to understand the changes made when increasing the levels of parallelism.  But another reason is that the PIC community has not yet agreed on standard objects and it would be more difficult to reuse software from the skeleton codes if the types or structs we chose conflicted.

###Overview of What's Here
In the following discussion, timings refer to the 2.67 GHz Intel i7 CPUs (Dual Intel Xeon-X5650), unless otherwise mentioned.
###Skeleton Codes
##Basic Serial Codes
The basic serial codes do not make use of any parallelism, and are the base codes for students or researchers who are unfamiliar with PIC codes.  For the electrostatic code there are 8 procedures in the inner loop, only 3 of which are computationally intensive, and this code contains about 2,000 lines of code in each language.   A typical execution time for the particle part of this code is about 33 ns/particle/time-step. For the electromagnetic code there are 15 procedures in the inner loop, only 4 of which are computationally intensive, and this code contains approximately 4,000 lines of code in each language.  A typical execution time for the particle part of this code is about 110 ns/particle/time-step.
##MPI Codes with one level of parallelism
These codes illustrate how to use domain decomposition with message-passing (MPI).  This is the dominate programming paradigm for PIC codes today.  These codes use only a simple 1d domain decomposition.  The algorithm is described in detail in Refs. [2-3].  For the electrostatic code, a typical execution time for the particle part of this code is about 140 ps/particle/time-step with 256 MPI nodes.  For the electromagnetic code, a typical execution time for the particle part of this code is about 400 ps/particle/time-step with 256 MPI nodes.  The CPUs were throttled down to 1.6 GHz for these benchmarks.
##OpenMP Codes with one level of parallelism
These codes illustrate how to use shared memory parallelism with OpenMP.  The algorithm used here is the same as that used on the GPU: it divides particles into small 2d tiles and reorders them every time step to avoid data collisions when performing the deposit.  Each tile is controlled by a single thread.  The algorithm is described in detail in Ref. [4],  For the electrostatic code, a typical execution time for the particle part of this code is about 3 ns/particle/time-step with 12 processing cores.  For the electromagnetic code, a typical execution time for the particle part of this code is about 10 ns/particle/time-step with 12 processing cores.
##MPI/OpenMP Codes with two levels of parallelism
These codes illustrate how to use a hybrid shared/distributed memory algorithm, with a tiled scheme on each shared memory multi-core node implemented with OpenMP, and domain decomposition connecting such nodes implemented with MPI.  The algorithms are described in detail in Refs. [2-4],  For the electrostatic code, a typical execution time for the particle part of this code is about 80 ps/particle/time-step with 576 processing cores.  For the electromagnetic code, a typical execution time for the particle part of this code is about 230 ps/particle/time-step with 576 processing cores.  The CPUs were throttled down to 1.6 GHz for these benchmarks.
##GPU Codes with two levels of parallelism
These codes illustrate how to implement an optimal PIC algorithm on a GPU.  It uses a hybrid tiling scheme with SIMD vectorization, both written with NVIDIA’s CUDA programming environment.  The tiling algorithm used within a thread block on the GPU is the same as that used with OpenMP [4].  Unlike the OpenMP implementation, however, each tile is controlled by a block of threads rather than a single thread, which requires a vectorized or data parallel implementation. Both CUDA C and CUDA Fortran interoperable versions are available, where a Fortran code can call the CUDA C libraries and a C code can call the CUDA Fortran libraries.  For the electrostatic code, a typical execution time for the particle part of this code on the M2090 GPU is about 0.9 ns/particle/time-step.  For the electromagnetic code, a typical execution time for the particle part of this code is about 2.0 ns/particle/time-step.
##GPU-MPI Codes with three levels of parallelism
These codes illustrate how to implement an optimal PIC algorithm on multiple GPUs.  It uses a hybrid tiling scheme with SIMD vectorization on each GPU, and domain decomposition connecting such GPUs implemented with MPI.  The algorithms used are described in Refs. [2-4].  Both CUDA C and CUDA Fortran interoperable versions are available, where a Fortran code can call the CUDA C libraries and a C code can call the CUDA Fortran libraries.  For the electrostatic code, a typical execution time for the particle part of this code on 108 M2090 GPUs is about 13 ps/particle/time-step.  For the electromagnetic code, a typical execution time for the particle part of this code is about 30 ps/particle/time-step.  To put this into perspective, on 96 GPUs, a 2-1/2d electromagnetic simulation with 10 billion particles and a 16384x16384 grid takes about 0.5 sec/time step, including the field solver.  The particle calculation is 3600 times faster than in the original serial code.
##Future Work
      The next focus of this project is to study the application of vectorizing technologies for PIC, with the immediate focus on the Intel Phi co-processor.  After that, there are 3 different topics we plan to study.  One is to create advanced versions of some of these skeleton codes to include features such as dynamic load balancing.  Another is to implement 3D versions of some of the skeleton codes.  Finally, we plan to explore alternative programming paradigms, such as OpenCL, OpenACC, or co-array Fortran, perhaps with other partners.
###References
[1]  Viktor K. Decyk, "A Method for Passing Data Between C and Opaque Fortran90 Pointers," ACM Fortran Forum, vol. 27, no. 2, p. 2 (2008).

[2] P. C. Liewer and V. K. Decyk, “A General Concurrent Algorithm for Plasma Particle-in-Cell Codes,” J. Computational Phys. 85, 302 (1989).

[3] V. K. Decyk, "Skeleton PIC Codes for Parallel Computers," Computer Physics Communications, 87, 87 (1995).

[4] Viktor K. Decyk and Tajendra V. Singh, “Particle-in-Cell algorithms for emerging computer architectures,” Computer Physics Communications 185, 708 (2014).



