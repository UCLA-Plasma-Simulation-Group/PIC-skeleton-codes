/* AVX512 Fortran90 utility Library */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <complex.h>
#include "avx512lib3.h"

/* Prototypes for Fortran90 functions called by C */

void getf2cptr_(unsigned long *carrayref, float *carray, int *nx,
                int *ny);

void getf3cptr_(unsigned long *carrayref, float *carray, int *nx,
                int *ny, int *nz);

void getf4cptr_(unsigned long *carrayref, float *carray, int *ndim, 
                int *nx, int *ny, int *nz);

void getc2cptr_(unsigned long *carrayref, float complex *carray,
                int *nx, int *ny);

void getc3cptr_(unsigned long *carrayref, float complex *carray,
                int *nx, int *ny, int *nz);

void geti1cptr_(unsigned long *carrayref, int *carray, int *nx);

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void avx512_f2allocatex_(unsigned long *sp_f, int *nx, int *ny, int *nd,
                         int *irc) {
/* allocate aligned 2d real memory, assign data   */
/* pointer to Fortran pointer object sp_f         */
/* nd = (1,2) = (default,double precision)        */
/* This procedure needs an interface in Fortran90 */
/* interface                                           */
/*    subroutine avx512_f2allocatex(sp_f,nx,ny,nd,irc) */
/*    implicit none                                    */
/*    integer :: nx, ny, nd, irc                       */
/*    real, dimension(:,:), pointer :: sp_f            */
/*    end subroutine                                   */
/* end interface                                       */
/* local data */
   int nsize;
   float *fptr;
   nsize = (*nx)*(*ny);
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_fallocate(&fptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getf2cptr_(sp_f,fptr,nx,ny);
   return;
}

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
void avx512_c2allocatex_(unsigned long *sp_c, int *nx, int *ny, int *nd,
                         int *irc) {
/* allocate aligned 2d complex memory, assign data */
/* pointer to Fortran pointer object sp_c          */
/* nd = (1,2) = (default,double precision)         */
/* This procedure needs an interface in Fortran90  */
/* interface                                           */
/*    subroutine avx512_c2allocatex(sp_c,nx,ny,nd,irc) */
/*    implicit none                                    */
/*    integer :: nx, ny, nd, irc                       */
/*    complex, dimension(:,:), pointer :: sp_c         */
/*    end subroutine                                   */
/* end interface                                       */
/* local data */
   int nsize;
   float complex *cptr;
   nsize = (*nx)*(*ny);
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_callocate(&cptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getc2cptr_(sp_c,cptr,nx,ny);
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

/*--------------------------------------------------------------------*/
void avx512_i1allocatex_(unsigned long *sp_i, int *nx,  int *nd,
                         int *irc) {
/* allocate aligned 1d integer memory, assign data */
/* pointer to Fortran pointer object sp_i          */ 
/* This procedure needs an interface in Fortran90  */
/* nd = (1,2) = (default,double precision)         */
/* interface                                        */
/*    subroutine avx512_i1allocatex(sp_i,nx,nd,irc) */
/*    implicit none                                 */
/*    integer :: nx, nd, irc                        */
/*    real, dimension(:), pointer :: sp_i           */
/*    end subroutine                                */
/* end interface                                    */
/* local data */
   int nsize;
   int *iptr;
   nsize = *nx;
   if ((*nd >= 0) && (*nd <= 2))
      nsize = (*nd)*nsize;
   avx512_iallocate(&iptr,nsize,irc);
/* set reference to C data in integer Fortran pointer object */
   geti1cptr_(sp_i,iptr,nx);
   return;
}
