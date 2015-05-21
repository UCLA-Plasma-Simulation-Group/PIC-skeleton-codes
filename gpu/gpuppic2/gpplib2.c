/*--------------------------------------------------------------------*/
/* Basic parallel PIC library for GPU-MPI communications
   gpplib2.c contains basic communications procedures for 1d partitions:
   cgppcaguard2l accumulates guard cells and copies to scalar field
   cgppccguard2l replicates guard cells and copies to 2 component vector
!                field
   cwappfft2rcs performs real to complex asynchronous fft for scalar
                array
   cwappfft2rcsn performs real to complex asynchronous fft for vector
                 array
   gpuppfft2rrcu performs real to complex asynchronous fft for scalar
                 array, based on NVIDIA FFT
   gpuppfft2rrcun performs real to complex asynchronous fft for vector
                  array, based on NVIDIA FFT
   cgpporder2l sorts partiles by tiles
   cgpptpose performs a transpose of a complex scalar array, distributed
             in y, to a complex scalar array, distributed in x.
             data from GPU is sent asynchronous, overlapping with MPI
   cgpptposen performs a transpose of an n component complex vector array,
              distributed in y, to an n component complex vector array,
              distributed in x.
              data from GPU is sent asynchronous, overlapping with MPI
   written by viktor k. decyk, ucla
   copyright 2013, regents of the university of california
   update: april 28, 2014                                              */

#include <complex.h>
#include <sys/time.h>
#include "gpulib2.h"
#include "gpuppush2.h"
#include "gpupfft2.h"
#include "pplib2.h"
#include "gpplib2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

/*--------------------------------------------------------------------*/
void cgppcaguard2l(float complex g_q[], float g_qe[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int nxe, int nypmx, int nxvh, int kypd) {
/* this subroutine copies scalar field and accumulates guard cells */
/* in y from remote GPU into scalar field                          */
   cgpuppcaguard2xl(g_q,g_scs,g_qe,nyp,nx,nxe,nypmx,nxvh,kypd);
   gpu_fcopyout(scs,g_scs,2*nxvh);
   cpppnaguard2l(scs,scr,kstrt,nvp,2*nxvh);
   gpu_fcopyin(scr,g_scs,2*nxvh);
   cgpuppcaguard2yl(g_q,g_scs,nx,nxvh,kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cgppccguard2l(float complex g_fxy[], float g_fxye[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int ndim, int nxe, int nypmx, int nxvh,
                   int kypd) {
/* this subroutine copies 2 component vector field and adds additional */
/* guard cells in y from remote GPU into extended vector field         */
   int nnxe;
   nnxe = ndim*nxe;
   cgpuppccguard2xl(g_fxy,g_scs,g_fxye,nyp,nx,nxe,nypmx,nxvh,kypd);
   gpu_fcopyout(scs,g_scs,2*nxvh*ndim);
   cpppncguard2l(scs,scr,kstrt,nvp,nnxe);
   gpu_fcopyin(scr,g_scs,2*nxvh*ndim);
   cgpuppccguard2yl(g_fxye,g_scs,nyp,nx,nxe,nxvh,nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cwappfft2rcs(float complex g_f[], float complex g_g[],
                  float complex g_bsm[], float complex g_brm[],
                  float complex bsm[], float complex brm[], int isign,
                  int g_mixup[], float complex g_sct[], float ttp[],
                  int indx, int indy, int kstrt, int nvp, int kxpd,
                  int kyp, int nxhd, int nyd, int kypd, int nxhyd,
                  int nxyhd) {
/* wrapper function for gpu-mpi parallel real to complex fft,
   without packed data
   if isign = -1, g_f = input, g_g = output
   if isign = 1, g_g = input, g_f = output
   nxhd must be = nx/2 + 1
local data */
/* kasync = (0,1) = (no,yes) use asynchronous communications */
#define KASYNC             1
   int nx, nxh1, ny;
   float ani, time;
   double dtime;
   struct timeval itime;
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ani = 1.0/(((float) nx)*((float) ny));
/* inverse fourier transform */
   if (isign < 0) {
/* first transpose in x */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsx(g_f,g_bsm,isign,g_mixup,g_sct,indx,indy,kstrt,
                      nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd);
/* transpose on local GPU */
      cgpuppltpose(g_f,g_g,nxhd,ny,kxpd,kyp,kstrt,nxhd,nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptpose(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,nvp);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(bsm,g_bsm,kxpd*kyp*(nvp-1));
            cppptpose(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp);
            gpu_ccopyin(brm,g_brm,kxpd*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in y */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsy(g_g,g_brm,isign,g_mixup,g_sct,indx,indy,kstrt,
                      nvp,kxpd,kyp,nyd,nxhyd,nxyhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
/* forward fourier transform */
   else if (isign > 0) {
/* first transpose in y */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsy(g_g,g_brm,isign,g_mixup,g_sct,indx,indy,kstrt,
                      nvp,kxpd,kyp,nyd,nxhyd,nxyhd);
/* transpose on local GPU */
      cgpuppltpose(g_g,g_f,ny,nxhd,kyp,kxpd,kstrt,nyd,nxhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptpose(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,nvp);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(brm,g_brm,kxpd*kyp*(nvp-1));
            cppptpose(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp);
            gpu_ccopyin(bsm,g_bsm,kxpd*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in x */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsx(g_f,g_bsm,isign,g_mixup,g_sct,indx,indy,kstrt,
                      nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
   return;
#undef KASYNC
}

/*--------------------------------------------------------------------*/
void cwappfft2rcsn(float complex g_fn[], float complex g_gn[],
                   float complex g_bsm[], float complex g_brm[],
                   float complex bsm[], float complex brm[], int isign,
                   int g_mixup[], float complex g_sct[], float ttp[], 
                   int indx, int indy, int kstrt, int nvp, int ndim,
                   int kxpd, int kyp, int nxhd, int nyd, int kypd,
                   int nxhyd, int nxyhd) {
/* wrapper function for multiple gpu-mpi parallel real to complex ffts,
   without packed data
   if isign = -1, g_fn = input, g_gn = output
   if isign = 1, g_gn = input, g_fn = output
   ndim = vector dimension
   nxhd must be = nx/2 + 1
local data */
/* kasync = (0,1) = (no,yes) use asynchronous communications */
#define KASYNC             1
   int nx, nxh1, ny;
   float ani, time;
   double dtime;
   struct timeval itime;
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ani = 1.0/(((float) nx)*((float) ny));
/* inverse fourier transform */
   if (isign < 0) {
/* first transpose in x */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsxn(g_fn,g_bsm,isign,g_mixup,g_sct,indx,indy,ndim,
                       kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd);
/* transpose on local GPU */
      cgpuppltposen(g_fn,g_gn,nxhd,ny,kxpd,kyp,kstrt,ndim,nxhd,nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptposen(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,nvp,
                       ndim);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(bsm,g_bsm,kxpd*ndim*kyp*(nvp-1));
            cppptposen(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp,ndim);
            gpu_ccopyin(brm,g_brm,kxpd*ndim*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in y */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsyn(g_gn,g_brm,isign,g_mixup,g_sct,indx,indy,ndim,
                       kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
/* forward fourier transform */
   else if (isign > 0) {
/* first transpose in y */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsyn(g_gn,g_brm,isign,g_mixup,g_sct,indx,indy,ndim,
                       kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd);
/* transpose on local GPU */
      cgpuppltposen(g_gn,g_fn,ny,nxhd,kyp,kxpd,kstrt,ndim,nyd,nxhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptposen(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,nvp,
                       ndim);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(brm,g_brm,kxpd*ndim*kyp*(nvp-1));
            cppptposen(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp,ndim);
            gpu_ccopyin(bsm,g_bsm,kxpd*ndim*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in x */
      dtimer(&dtime,&itime,-1);
      cgpuwppfft2rcsxn(g_fn,g_bsm,isign,g_mixup,g_sct,indx,indy,ndim,
                       kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
   return;
#undef KASYNC
}

/*--------------------------------------------------------------------*/
void gpuppfft2rrcu(float complex g_f[], float complex g_g[],
                   float complex g_bsm[], float complex g_brm[],
                   float complex bsm[], float complex brm[], int isign,
                   float ttp[], int indx, int indy, int kstrt, int nvp,
                   int kxpd, int kyp, int nxhd, int nyd, int kypd) {
/* wrapper function for gpu-mpi parallel real to complex fft,
   based on NVIDIA FFT, without packed data
   if isign = -1, g_f = input, g_g = output
   if isign = 1, g_g = input, g_f = output
   nxhd must be = nx/2 + 1
local data */
/* kasync = (0,1) = (no,yes) use asynchronous communications */
#define KASYNC             1
   int nx, nxh1, ny;
   float ani, time;
   double dtime;
   struct timeval itime;
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ani = 1.0/(((float) nx)*((float) ny));
/* inverse fourier transform */
   if (isign < 0) {
/* first transpose in x */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcux(g_f,g_bsm,isign,indx,indy,kstrt,nvp,kxpd,kyp,
                    nxhd,kypd);
/* transpose on local GPU with scaling */
      cgpuppsltpose(g_f,g_g,ani,nxhd,ny,kxpd,kyp,kstrt,nxhd,nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptpose(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,nvp);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(bsm,g_bsm,kxpd*kyp*(nvp-1));
            cppptpose(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp);
            gpu_ccopyin(brm,g_brm,kxpd*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in y */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuy(g_g,g_brm,isign,indx,indy,kstrt,nvp,kxpd,kyp,nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
/* forward fourier transform */
   else if (isign > 0) {
/* first transpose in y */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuy(g_g,g_brm,isign,indx,indy,kstrt,nvp,kxpd,kyp,nyd);
/* transpose on local GPU */
      cgpuppltpose(g_g,g_f,ny,nxhd,kyp,kxpd,kstrt,nyd,nxhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptpose(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,nvp);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(brm,g_brm,kxpd*kyp*(nvp-1));
            cppptpose(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp);
            gpu_ccopyin(bsm,g_bsm,kxpd*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in x */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcux(g_f,g_bsm,isign,indx,indy,kstrt,nvp,kxpd,kyp,
                    nxhd,kypd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
   return;
#undef KASYNC
}

/*--------------------------------------------------------------------*/
void gpuppfft2rrcun(float complex g_fn[], float complex g_gn[],
                    float complex g_bsm[], float complex g_brm[],
                    float complex bsm[], float complex brm[], int isign,
                    float ttp[], int indx, int indy, int kstrt, int nvp,
                    int ndim, int kxpd, int kyp, int nxhd, int nyd,
                    int kypd) {
/* wrapper function for multiple gpu-mpi parallel real to complex ffts,
   based on NVIDIA FFT, without packed data
   if isign = -1, g_fn = input, g_gn = output
   if isign = 1, g_gn = input, g_fn = output
   ndim = vector dimension
   nxhd must be = nx/2 + 1
local data */
/* kasync = (0,1) = (no,yes) use asynchronous communications */
#define KASYNC             1
   int nx, nxh1, ny;
   float ani, time;
   double dtime;
   struct timeval itime;
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ani = 1.0/(((float) nx)*((float) ny));
/* inverse fourier transform */
   if (isign < 0) {
/* first transpose in x */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuxn(g_fn,g_bsm,isign,indx,indy,ndim,kstrt,nvp,kxpd,
                     kyp,nxhd,kyp);
/* transpose on local GPU with scaling */
      cgpuppsltposen(g_fn,g_gn,ani,nxhd,ny,kxpd,kyp,kstrt,ndim,nxhd,
                     nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptposen(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,nvp,
                       ndim);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(bsm,g_bsm,kxpd*ndim*kyp*(nvp-1));
            cppptposen(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp,ndim);
            gpu_ccopyin(brm,g_brm,kxpd*ndim*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in y */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuyn(g_gn,g_brm,isign,indx,indy,ndim,kstrt,nvp,kxpd,
                     kyp,nyd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
/* forward fourier transform */
   else if (isign > 0) {
/* first transpose in y */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuyn(g_gn,g_brm,isign,indx,indy,ndim,kstrt,nvp,kxpd,
                     kyp,nyd);
/* transpose on local GPU */
      cgpuppltposen(g_gn,g_fn,ny,nxhd,kyp,kxpd,kstrt,ndim,nyd,nxhd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
/* transpose between GPUs */
      if (nvp > 1) {
         dtimer(&dtime,&itime,-1);
/* use asynchronous communication */
         if (KASYNC==1) {
            cgpptposen(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,nvp,
                       ndim);
         }
/* use synchronous communication */
         else {
            gpu_ccopyout(brm,g_brm,kxpd*ndim*kyp*(nvp-1));
            cppptposen(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp,ndim);
            gpu_ccopyin(bsm,g_bsm,kxpd*ndim*kyp*(nvp-1));
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         ttp[1] += time;
      }
/* then transpose in x */
      dtimer(&dtime,&itime,-1);
      gpupfft2rrcuxn(g_fn,g_bsm,isign,indx,indy,ndim,kstrt,nvp,kxpd,
                     kyp,nxhd,kypd);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      ttp[0] += time;
   }
   return;
#undef KASYNC
}

/*--------------------------------------------------------------------*/
void cgpporder2l(float g_ppart[], float g_ppbuff[], float g_sbufl[],
                 float g_sbufr[], int g_kpic[], int g_ncl[],
                 int g_ihole[], int g_ncll[], int g_nclr[],
                 float sbufl[], float sbufr[], float rbufl[],
                 float rbufr[], int ncll[], int nclr[], int mcll[],
                 int mclr[], float ttp[], int noff, int nyp, int kstrt,
                 int nvp, int idimp, int nppmx, int nx, int ny, int mx,
                 int my, int mx1, int myp1, int npbmx, int ntmax,
                 int nbmax, int *g_irc) {
/* this subroutine performs an mpi-gpu particle sort by x,y grid in tiles
   of mx, my
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
local data */
   float time;
   double dtime;
   struct timeval itime;
/* first part of particle reorder on x and y cell with mx, my tiles */
   dtimer(&dtime,&itime,-1);
   cgpupppord2la(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,g_ihole,
                 g_ncll,g_nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,
                 myp1,npbmx,ntmax,nbmax,g_irc);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[0] += time;
/* move particles on GPU into appropriate spatial regions */
   dtimer(&dtime,&itime,-1);
   gpu_icopyout(ncll,g_ncll,3*mx1);
   gpu_icopyout(nclr,g_nclr,3*mx1);
   gpu_fcopyout(sbufl,g_sbufl,idimp*ncll[3*mx1-1]);
   gpu_fcopyout(sbufr,g_sbufr,idimp*nclr[3*mx1-1]);
   cpppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,
             idimp,nbmax,mx1);
   gpu_icopyin(mcll,g_ncll,3*mx1);
   gpu_icopyin(mclr,g_nclr,3*mx1);
   gpu_fcopyin(rbufl,g_sbufl,idimp*mcll[3*mx1-1]);
   gpu_fcopyin(rbufr,g_sbufr,idimp*mclr[3*mx1-1]);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[1] += time;
/* second part of particle reorder on x and y cell with mx, my tiles */
   dtimer(&dtime,&itime,-1);
   cgpupppord2lb(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,g_ihole,
                 g_ncll,g_nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,
                 g_irc);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[0] += time;
   return;
}

/*--------------------------------------------------------------------*/
void cgpptpose(float complex g_bsm[], float complex g_btm[],
               float complex sm[], float complex tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp) {
/* this subroutine sends and receives data between GPUS on different MPI
   nodes to perform a transpose of a matrix distributed in y, to another
   matrix distributed in x.
   one message is sent and received at a time.
   data from GPU is sent asynchronous, overlapping with MPI
   g_bsm/g_btm are complex buffers on GPU to be sent/received
   sm/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   kxyp = kxp*kyp;
/* special case for one processor  */
/* better to use a kernel function */
   if (nvp==1) {
      gpu_ccopyout(sm,g_bsm,kxyp);
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      gpu_ccopyin(tm,g_btm,kxyp);
      return;
   }
   nn = -1;
   ns = kxyp;
   stp = 1; st = 2;
/* send first group to host from GPU */
   gpu_cascopyout(sm,g_bsm,0,ns,stp);
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* send next group to host from GPU */
         if (nn < nvp) {
            st = stp + 1;
            if (st > 2)
               st -= 2;
            gpu_cascopyout(&sm[kxyp*(nn+1)],g_bsm,ns*(nn+1),ns,st);
         }
/* wait for previous MPI sends and receives to complete */
         if (nn > 0) {
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
            cacsndrec(sm,0,0,0,3);
/* copy received group from host to GPU */
            gpu_cascopyin(&tm[kxyp*(nn-1)],g_btm,ns*(nn-1),ns,3);
         }
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
/* ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&mrid); */
         cacsndrec(&tm[kxyp*nn],id,kxyp,n+1,1);
/* wait for previous group to arrive from GPU */
         gpu_waitstream(stp);
         stp = st;
/* send data */
/* ierr = MPI_Isend(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp,&msid); */
         cacsndrec(&sm[kxyp*nn],id,ld,n+1,2);
      }
   }
/* wait for last sends and receives to complete */
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
   cacsndrec(sm,0,0,0,3);
   gpu_cascopyin(&tm[kxyp*nn],g_btm,ns*nn,ns,3);
/* wait for last group item to arrive */
   gpu_waitstream(3);
   return;
}

/*--------------------------------------------------------------------*/
void cgpptposen(float complex g_bsm[], float complex g_btm[],
                float complex sm[], float complex tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim) {
/* this subroutine sends and receives data between GPUS on different MPI
   nodes to perform a transpose of an n component  matrix distributed in
   y, to another an n component matrix distributed in x.
   one message is sent and received at a time.
   data from GPU is sent asynchronous, overlapping with MPI.
   g_bsm/g_btm are complex buffers on GPU to be sent/received
   sm/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int  j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = ndim*(kyp < kyps ? kyp : kyps);
   kxyp = kxp*ndim*kyp;
/* special case for one processor  */
/* better to use a kernel function */
   if (nvp==1) {
      gpu_ccopyout(sm,g_bsm,kxyp);
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      gpu_ccopyin(tm,g_btm,kxyp);
      return;
   }
   nn = -1;
   ns = kxyp;
   stp = 1; st = 2;
/* send first group to host from GPU */
   gpu_cascopyout(sm,g_bsm,0,ns,stp);
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* send next group to host from GPU */
         if (nn < nvp) {
            if (st > 2)
               st -= 2;
            gpu_cascopyout(&sm[kxyp*(nn+1)],g_bsm,ns*(nn+1),
                              ns,st);
         }
/* wait for previous MPI sends and receives to complete */
         if (nn > 0) {
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
            cacsndrec(sm,0,0,0,3);
/* copy received group from host to GPU */
            gpu_cascopyin(&tm[kxyp*(nn-1)],g_btm,ns*(nn-1),
                             ns,3);
         }
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
/* ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&mrid); */
         cacsndrec(&tm[kxyp*nn],id,kxyp,n+1,1);
/* wait for previous group to arrive from GPU */
         gpu_waitstream(stp);
         stp = st;
/* send data */
/* ierr = MPI_Isend(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp,&msid); */
         cacsndrec(&sm[kxyp*nn],id,ld,n+1,2);
      }
   }
/* wait for last sends and receives to complete */
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
   cacsndrec(sm,0,0,0,3);
   gpu_cascopyin(&tm[kxyp*nn],g_btm,ns*nn,ns,3);
/* wait for last group item to arrive */
   gpu_waitstream(3);
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cgppcaguard2l_(unsigned long *gp_q, unsigned long *gp_qe,
                    unsigned long *gp_scs, float *scs, float *scr,
                    int *nx, int *nyp, int *kstrt, int *nvp, int *nxe,
                    int *nypmx, int *nxvh, int *kypd) {
   float complex *g_q;
   float *g_qe, *g_scs;
   g_q = (float complex *)*gp_q;
   g_qe = (float *)*gp_qe;
   g_scs = (float *)*gp_scs;
   cgppcaguard2l(g_q,g_qe,g_scs,scs,scr,*nx,*nyp,*kstrt,*nvp,*nxe,
                 *nypmx,*nxvh,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cgppccguard2l_(unsigned long *gp_fxy, unsigned long *gp_fxye,
                    unsigned long *gp_scs, float *scs, float *scr,
                    int *nx, int *nyp, int *kstrt, int *nvp, int *ndim,
                    int *nxe, int *nypmx, int *nxvh, int *kypd) {
   float complex *g_fxy;
   float *g_fxye, *g_scs;
   g_fxy = (float complex *)*gp_fxy;
   g_fxye = (float *)*gp_fxye;
   g_scs = (float *)*gp_scs;
   cgppccguard2l(g_fxy,g_fxye,g_scs,scs,scr,*nx,*nyp,*kstrt,*nvp,*ndim,
                 *nxe,*nypmx,*nxvh,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cwappfft2rcs_(unsigned long *gp_f, unsigned long *gp_g,
                   unsigned long *gp_bsm, unsigned long *gp_brm,
                   float complex *bsm, float complex *brm, int *isign,
                   unsigned long *gp_mixup, unsigned long *gp_sct,
                   float *ttp, int *indx, int *indy, int *kstrt,
                   int *nvp, int *kxpd, int *kyp, int *nxhd, int *nyd,
                   int *kypd, int *nxhyd, int *nxyhd) {
   float complex *g_f, *g_g, *g_bsm, *g_brm, *g_sct;
   int *g_mixup;
   g_f = (float complex *)*gp_f;
   g_g = (float complex *)*gp_g;
   g_bsm = (float complex *)*gp_bsm;
   g_brm = (float complex *)*gp_brm;
   g_mixup = (int *)*gp_mixup;
   g_sct = (float complex *)*gp_sct;
   cwappfft2rcs(g_f,g_g,g_bsm,g_brm,bsm,brm,*isign,g_mixup,g_sct,ttp,
                *indx,*indy,*kstrt,*nvp,*kxpd,*kyp,*nxhd,*nyd,*kypd,
                *nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwappfft2rcsn_(unsigned long *gp_fn, unsigned long *gp_gn,
                    unsigned long *gp_bsm, unsigned long *gp_brm,
                    float complex *bsm, float complex *brm, int *isign,
                    unsigned long *gp_mixup, unsigned long *gp_sct,
                    float *ttp, int *indx, int *indy, int *kstrt,
                    int *nvp, int *ndim, int *kxpd, int *kyp, int *nxhd,
                    int *nyd, int *kypd, int *nxhyd, int *nxyhd) {
   float complex *g_fn, *g_gn, *g_bsm, *g_brm, *g_sct;
   int *g_mixup;
   g_fn = (float complex *)*gp_fn;
   g_gn = (float complex *)*gp_gn;
   g_bsm = (float complex *)*gp_bsm;
   g_brm = (float complex *)*gp_brm;
   g_mixup = (int *)*gp_mixup;
   g_sct = (float complex *)*gp_sct;
   cwappfft2rcsn(g_fn,g_gn,g_bsm,g_brm,bsm,brm,*isign,g_mixup,g_sct,ttp,
                 *indx,*indy,*kstrt,*nvp,*ndim,*kxpd,*kyp,*nxhd,*nyd,
                 *kypd,*nxhyd,*nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void gpuppfft2rrcu_(unsigned long *gp_f, unsigned long *gp_g,
                    unsigned long *gp_bsm, unsigned long *gp_brm,
                    float complex *bsm, float complex *brm, int *isign,
                    float *ttp, int *indx, int *indy, int *kstrt,
                    int *nvp, int *kxpd, int *kyp, int *nxhd, int *nyd,
                    int *kypd) {
   float complex *g_f, *g_g, *g_bsm, *g_brm;
   g_f = (float complex *)*gp_f;
   g_g = (float complex *)*gp_g;
   g_bsm = (float complex *)*gp_bsm;
   g_brm = (float complex *)*gp_brm;
   gpuppfft2rrcu(g_f,g_g,g_bsm,g_brm,bsm,brm,*isign,ttp,*indx,*indy,
                 *kstrt,*nvp,*kxpd,*kyp,*nxhd,*nyd,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
void gpuppfft2rrcun_(unsigned long *gp_fn, unsigned long *gp_gn,
                     unsigned long *gp_bsm, unsigned long *gp_brm,
                     float complex *bsm, float complex *brm, int *isign,
                     float *ttp, int *indx, int *indy, int *kstrt,
                     int *nvp, int *ndim, int *kxpd, int *kyp, int *nxhd,
                     int *nyd, int *kypd) {
   float complex *g_fn, *g_gn, *g_bsm, *g_brm;
   g_fn = (float complex *)*gp_fn;
   g_gn = (float complex *)*gp_gn;
   g_bsm = (float complex *)*gp_bsm;
   g_brm = (float complex *)*gp_brm;
   gpuppfft2rrcun(g_fn,g_gn,g_bsm,g_brm,bsm,brm,*isign,ttp,*indx,*indy,
                  *kstrt,*nvp,*ndim,*kxpd,*kyp,*nxhd,*nyd,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cgpporder2l_(unsigned long *gp_ppart, unsigned long *gp_ppbuff,
                  unsigned long *gp_sbufl, unsigned long *gp_sbufr,
                  unsigned long *gp_kpic, unsigned long *gp_ncl,
                  unsigned long *gp_ihole, unsigned long *gp_ncll,
                  unsigned long *gp_nclr, float *sbufl, float *sbufr,
                  float *rbufl, float *rbufr, int *ncll, int *nclr,
                  int *mcll, int *mclr, float *ttp, int *noff, int *nyp,
                  int *kstrt, int *nvp, int *idimp, int *nppmx, int *nx,
                  int *ny, int *mx, int *my, int *mx1, int *myp1,
                  int *npbmx, int *ntmax, int *nbmax,
                  unsigned long *gp_irc) {
   float *g_ppart, *g_ppbuff, *g_sbufl, *g_sbufr;
   int *g_kpic, *g_ncl, *g_ihole, *g_ncll, *g_nclr, *g_irc;
   g_ppart = (float *)*gp_ppart;
   g_ppbuff = (float *)*gp_ppbuff;
   g_sbufl = (float *)*gp_sbufl;
   g_sbufr = (float *)*gp_sbufr;
   g_kpic = (int *)*gp_kpic;
   g_ncl = (int *)*gp_ncl;
   g_ihole = (int *)*gp_ihole;
   g_ncll = (int *)*gp_ncll;
   g_nclr = (int *)*gp_nclr;
   g_irc = (int *)*gp_irc;
   cgpporder2l(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,g_ihole,
               g_ncll,g_nclr,sbufl,sbufr,rbufl,rbufr,ncll,nclr,mcll,
               mclr,ttp,*noff,*nyp,*kstrt,*nvp,*idimp,*nppmx,*nx,*ny,
               *mx,*my,*mx1,*myp1,*npbmx,*ntmax,*nbmax,g_irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpptpose_(unsigned long *gp_bsm, unsigned long *gp_btm,
                float complex *sm, float complex *tm, int *nx, int *ny,
                int *kxp, int *kyp, int *kstrt, int *nvp) {
   float complex *g_bsm, *g_btm;
   g_bsm = (float complex *)*gp_bsm;
   g_btm = (float complex *)*gp_btm;
   cgpptpose(g_bsm,g_btm,sm,tm,*nx,*ny,*kxp,*kyp,*kstrt,*nvp);
   return;
}

/*--------------------------------------------------------------------*/
void cgpptposen_(unsigned long *gp_bsm, unsigned long *gp_btm,
                 float complex *sm, float complex *tm, int *nx, int *ny,
                 int *kxp, int *kyp, int *kstrt, int *nvp, int *ndim) {
   float complex *g_bsm, *g_btm;
   g_bsm = (float complex *)*gp_bsm;
   g_btm = (float complex *)*gp_btm;
   cgpptposen(g_bsm,g_btm,sm,tm,*nx,*ny,*kxp,*kyp,*kstrt,*nvp,*ndim);
   return;
}
