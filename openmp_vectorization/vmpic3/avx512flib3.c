/* AVX512 Fortran90 utility Library */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <complex.h>
#include "avx512lib3.h"

/* Prototypes for Fortran90 functions called by C */

void getf3cptr_(unsigned long *carrayref, float *carray, int *nx,
                int *ny, int *nz);

void getf4cptr_(unsigned long *carrayref, float *carray, int *ndim, 
                int *nx, int *ny, int *nz);

void getc3cptr_(unsigned long *carrayref, float complex *carray,
                int *nx, int *ny, int *nz);

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void avx512_f3allocatex_(unsigned long *sp_f, int *nx, int *ny, int *nz,
                         int *nd, int *irc) {
/* allocate aligned 3d real memory, assign data   */
/* pointer to Fortran pointer object sp_f         */
/* nd = (1,2) = (default,double precision)        */
/* This procedure needs an interface in Fortran90 */
/* interface                                              */
/*    subroutine avx512_f3allocatex(sp_f,nx,ny,nz,nd,irc) */
/*    implicit none                                       */
/*    integer :: nx, ny, nz, nd, irc                      */
/*    real, dimension(:,:,:), pointer :: sp_f             */
/*    end subroutine                                      */
/* end interface                                          */
/* local data */
   int nsize;
   float *fptr;
   nsize = (*nx)*(*ny)*(*nz);
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_fallocate(&fptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getf3cptr_(sp_f,fptr,nx,ny,nz);
   return;
}

/*--------------------------------------------------------------------*/
void avx512_f4allocatex_(unsigned long *sp_f, int *ndim, int *nx,
                         int *ny, int *nz, int *nd, int *irc) {
/* allocate aligned vector 3d real memory, assign data   */
/* pointer to Fortran pointer object sp_f                */
/* nd = (1,2) = (default,double precision)               */
/* This procedure needs an interface in Fortran90        */
/* interface                                                   */
/*    subroutine avx512_f4allocatex(sp_f,ndim,nx,ny,nz,nd,irc) */
/*    implicit none                                            */
/*    integer :: ndim, nx, ny, nz, nd, irc                     */
/*    real, dimension(:,:,:,:), pointer :: sp_f                */
/*    end subroutine                                           */
/* end interface                                               */
/* local data */
   int nsize;
   float *fptr;
   nsize = (*ndim)*(*nx)*(*ny)*(*nz);
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_fallocate(&fptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getf4cptr_(sp_f,fptr,ndim,nx,ny,nz);
   return;
}

/*--------------------------------------------------------------------*/
void avx512_c3allocatex_(unsigned long *sp_c, int *nx, int *ny, 
                         int *nz, int *nd, int *irc) {
/* allocate aligned 3d complex memory, assign data */
/* pointer to Fortran pointer object sp_c          */
/* nd = (1,2) = (default,double precision)         */
/* This procedure needs an interface in Fortran90  */
/* interface                                              */
/*    subroutine avx512_c3allocatex(sp_c,nx,ny,nz,nd,irc) */
/*    implicit none                                       */
/*    integer :: nx, ny, nz, nd, irc                      */
/*    complex, dimension(:,:,:), pointer :: sp_c          */
/*    end subroutine                                      */
/* end interface                                          */
/* local data */
   int nsize;
   float complex *cptr;
   nsize = (*nx)*(*ny)*(*nz);
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_callocate(&cptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getc3cptr_(sp_c,cptr,nx,ny,nz);
   return;
}

