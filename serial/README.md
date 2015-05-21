# PIC Skeleton Codes:  Serial

The basic serial codes do not make use of any parallelism, and are the base codes for students or researchers who are unfamiliar with PIC codes. For the 2D electrostatic code there are 8 procedures in the inner loop, only 3 of which are computationally intensive, and this code contains about 2,000 lines of code in each language. A typical execution time for the particle part of this code is about 33 ns/particle/time-step. For the 2-1/2D electromagnetic code there are 14 procedures in the inner loop, only 4 of which are computationally intensive, and this code contains approximately 4,000 lines of code in each language. A typical execution time for the particle part of this code is about 110 ns/particle/time-step.  For the 2-1/2D Darwin code there are 21 procedures in the inner loop, only 6 of which are computationally intensive, and this code contains approximately 4,000 lines of code in each language. A typical execution time for the particle part of this code is about 256 ns/particle/time-step.  (All timings are on a 2.67GHz Intel Nehalem processor.)

###Electrostatic:
1.  1D Serial Electrostatic Spectral code:  pic1
2.  2D Serial Electrostatic Spectral code:  pic2
3.  3D Serial Electrostatic Spectral code:  pic3

###Electromagnetic:
4.  1-2/2D Serial Electromagnetic Spectral code:  bpic1
5.  2-1/2D Serial Electromagnetic Spectral code:  bpic2
6.  3D Serial Electromagnetic Spectral code:  bpic3

###Darwin:
7.  1-2/2D Serial Darwin Spectral code:  dpic1
8.  2-1/2D Serial Darwin Spectral code:  dpic2
9.  3D Serial Darwin Spectral code:  dpic3



### Want to contact the developer?

Send mail to Viktor Decyk – decyk@physics.ucla.edu 


