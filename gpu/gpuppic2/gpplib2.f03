!-----------------------------------------------------------------------
! Basic parallel PIC library for GPU-MPI communications
! gpplib2.f03 contains basic communications procedures for 1d partitions
! GPPCAGUARD2L accumulates guard cells and copies to scalar field
! GPPCCGUARD2L replicates guard cells and copies to 2 component vector
!              field
! WAPPFFT2RCS performs real to complex asynchronous fft for scalar
!             array
! WAPPFFT2RCSN performs real to complex asynchronous fft for vector
!              array
! GPUPPFFT2RRCU performs real to complex asynchronous fft for scalar
!               array, based on NVIDIA FFT
! GPUPPFFT2RRCUN performs real to complex asynchronous fft for vector
!                array, based on NVIDIA FFT
! GPPORDER2L sorts partiles by tiles
! GPPTPOSE performs a transpose of a complex scalar array, distributed
!          in y, to a complex scalar array, distributed in x.
!          data from GPU is sent asynchronous, overlapping with MPI
! GPPTPOSEN performs a transpose of an n component complex vector array,
!           distributed in y, to an n component complex vector array,
!           distributed in x.
!           data from GPU is sent asynchronous, overlapping with MPI
! written by viktor k. decyk, ucla
! copyright 2013, regents of the university of california
! update: june 12, 2014
      module gpplib2_f03
      use iso_c_binding
      use gpulib2_c
      use gpuppush2_c
      use gpupfft2_c
      use pplib2       ! use with pplib2.f90
!     use pplib2_h     ! use with pplib2.f
      use dtimer_c
      implicit none
      private
      public :: GPPCAGUARD2L, GPPCCGUARD2L
      public :: WAPPFFT2RCS, WAPPFFT2RCSN, GPUPPFFT2RRCU, GPUPPFFT2RRCUN
      public :: GPPORDER2L
!
      contains
!
!-----------------------------------------------------------------------
      subroutine GPPCAGUARD2L(g_q,g_qe,g_scs,scs,scr,nx,nyp,kstrt,nvp,  &
     &nxe,nypmx,nxvh,kypd)
! this subroutine copies scalar field and accumulates guard cells in y
! from remote GPU into scalar field
      implicit none
      integer, intent(in) :: nx, nyp, kstrt, nvp, nxe, nypmx, nxvh, kypd
      type (c_ptr) :: g_q, g_qe, g_scs
!     complex dimension g_q(nxvh,kypd)
!     real dimension g_qe(nxe,nypmx)
!     real dimension g_scs(nxe)
      real, dimension(nxe), target :: scs, scr
      call cgpuppcaguard2xl(g_q,g_scs,g_qe,nyp,nx,nxe,nypmx,nxvh,kypd)
      call gpu_fcopyout(c_loc(scs(1)),g_scs,2*nxvh)
      call PPPNAGUARD2L(scs,scr,kstrt,nvp,2*nxvh)
      call gpu_fcopyin(c_loc(scr(1)),g_scs,2*nxvh)
      call cgpuppcaguard2yl(g_q,g_scs,nx,nxvh,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPPCCGUARD2L(g_fxy,g_fxye,g_scs,scs,scr,nx,nyp,kstrt,  &
     &nvp,ndim,nxe,nypmx,nxvh,kypd)
! this subroutine copies 2 component vector field and adds additional
! guard cells in y from remote GPU into extended vector field
      implicit none
      integer, intent(in) :: nx, nyp, kstrt, nvp, ndim, nxe, nypmx
      integer, intent(in) :: nxvh, kypd
      type (c_ptr) :: g_fxy, g_fxye, g_scs
!     complex dimension g_fxy(nxvh,ndim,kypd)
!     real dimension g_fxye(ndim,nxe,nypmx)
!     real dimension g_scs(nxe*ndim)
      real, dimension(nxe*ndim), target :: scs, scr
! local data
      integer nnxe
      nnxe = ndim*nxe
      call cgpuppccguard2xl(g_fxy,g_scs,g_fxye,nyp,nx,nxe,nypmx,nxvh,   &
     &kypd)
      call gpu_fcopyout(c_loc(scs(1)),g_scs,2*nxvh*ndim)
      call PPPNCGUARD2L(scs,scr,kstrt,nvp,nnxe)
      call gpu_fcopyin(c_loc(scr(1)),g_scs,2*nxvh*ndim)
      call cgpuppccguard2yl(g_fxye,g_scs,nyp,nx,nxe,nxvh,nypmx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine WAPPFFT2RCS(g_f,g_g,g_bsm,g_brm,bsm,brm,isign,g_mixup, &
     &g_sct,ttp,indx,indy,kstrt,nvp,kxpd,kyp,nxhd,nyd,kypd,nxhyd,nxyhd)
! wrapper function for gpu-mpi parallel real to complex fft,
! without packed data
! if isign = -1, g_f = input, g_g = output
! if isign = 1, g_g = input, g_f = output
! (g_bsm,g_brm)/(bsm,brm) are temporary scratch arrays on GPU/host
! nxhd must be = nx/2 + 1
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxpd, kyp
      integer, intent(in) :: nxhd, nyd, kypd, nxhyd, nxyhd
      type (c_ptr) :: g_f, g_g, g_bsm, g_brm, g_mixup, g_sct
!     complex dimension g_f(nxhd,kypd), g_g(nyd,kxpd)
!     complex dimension g_bsm(kxpd*kyp,nvp), g_brm(kxpd*kyp,nvp)
!     integer dimension g_mixup(nxhyd)
!     complex dimension g_sct(nxyhd)
      real, dimension(2) :: ttp
      complex, dimension(kxpd*kyp,nvp), target :: bsm, brm
! local data
! kasync = (0,1) = (no,yes) use asynchronous communications
      integer, parameter :: kasync = 1
      integer :: nx, nxh1, ny, nxyp
      real :: ani, time
      double precision :: dtime
      type (timeval) :: itime
      nx = 2**indx
      nxh1 = nx/2 + 1
      ny = 2**indy
      nxyp = kxpd*kyp*(nvp-1)
      ani = 1.0/(real(nx)*real(ny))
! inverse fourier transform
      if (isign < 0) then
         call dtimer(dtime,itime,-1)
! first transpose in x
         call cgpuwppfft2rcsx(g_f,g_bsm,isign,g_mixup,g_sct,indx,indy,  &
     &kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd)
! transpose on local GPU
         call cgpuppltpose(g_f,g_g,nxhd,ny,kxpd,kyp,kstrt,nxhd,nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSE(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,  &
     &nvp)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(bsm(1,1)),g_bsm,nxyp)
               call PPPTPOSE(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp)
               call gpu_ccopyin(c_loc(brm(1,1)),g_brm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in y
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsy(g_g,g_brm,isign,g_mixup,g_sct,indx,indy,  &
     &kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! forward fourier transform
      else if (isign > 0) then
         call dtimer(dtime,itime,-1)
! first transpose in y
         call cgpuwppfft2rcsy(g_g,g_brm,isign,g_mixup,g_sct,indx,indy,  &
     &kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd)
! transpose on local GPU
         call cgpuppltpose(g_g,g_f,ny,nxhd,kyp,kxpd,kstrt,nyd,nxhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
        if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSE(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,  &
     &nvp)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(brm(1,1)),g_brm,nxyp)
               call PPPTPOSE(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp)
               call gpu_ccopyin(c_loc(bsm(1,1)),g_bsm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in x
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsx(g_f,g_bsm,isign,g_mixup,g_sct,indx,indy,  &
     &kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine WAPPFFT2RCSN(g_fn,g_gn,g_bsm,g_brm,bsm,brm,isign,      &
     &g_mixup,g_sct,ttp,indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,nyd,kypd,&
     &nxhyd,nxyhd)
! wrapper function for multiple gpu-mpi parallel real to complex ffts,
! without packed data
! if isign = -1, g_fn = input, g_gn = output
! if isign = 1, g_gn = input, g_fn = output
! (g_bsm,g_brm)/(bsm,brm) are temporary scratch arrays on GPU/host
! ndim = vector dimension
! nxhd must be = nx/2 + 1
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, ndim, kxpd
      integer, intent(in) :: kyp, nxhd, nyd, kypd, nxhyd, nxyhd
      type (c_ptr) :: g_fn, g_gn, g_bsm, g_brm, g_mixup, g_sct
!     complex dimension g_fn(nxhd,ndim,kypd), g_gn(nyd,ndim,kxpd)
!     complex dimension g_bsm(kxpd*ndim*kyp,nvp)
!     complex dimension g_brm(kxpd*ndim*kyp,nvp)
!     integer dimension g_mixup(nxhyd)
!     complex dimension g_sct(nxyhd)
      real, dimension(2) :: ttp
      complex, dimension(kxpd*ndim*kyp,nvp), target :: bsm, brm
! local data
! kasync = (0,1) = (no,yes) use asynchronous communications
      integer, parameter :: kasync = 1
      integer :: nx, nxh1, ny, nxyp
      real :: ani, time
      double precision :: dtime
      type (timeval) :: itime
      nx = 2**indx
      nxh1 = nx/2 + 1
      ny = 2**indy
      nxyp = kxpd*ndim*kyp*(nvp-1)
      ani = 1.0/(real(nx)*real(ny))
! inverse fourier transform
      if (isign < 0) then
! first transpose in x
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsxn(g_fn,g_bsm,isign,g_mixup,g_sct,indx,indy,&
     &ndim,kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd)
! transpose on local GPU
         call cgpuppltposen(g_fn,g_gn,nxhd,ny,kxpd,kyp,kstrt,ndim,nxhd, &
     &nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSEN(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt, &
     &nvp,ndim)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(bsm(1,1)),g_bsm,nxyp)
               call PPPTPOSEN(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp,ndim)
               call gpu_ccopyin(c_loc(brm(1,1)),g_brm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in y
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsyn(g_gn,g_brm,isign,g_mixup,g_sct,indx,indy,&
     &ndim,kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! forward fourier transform
      else if (isign > 0) then
! first transpose in y
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsyn(g_gn,g_brm,isign,g_mixup,g_sct,indx,indy,&
     &ndim,kstrt,nvp,kxpd,kyp,nyd,nxhyd,nxyhd)
! transpose on local GPU
         call cgpuppltposen(g_gn,g_fn,ny,nxhd,kyp,kxpd,kstrt,ndim,nyd,  &
     &nxhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSEN(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt, &
     &nvp,ndim)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(brm(1,1)),g_brm,nxyp)
               call PPPTPOSEN(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp,ndim)
               call gpu_ccopyin(c_loc(bsm(1,1)),g_bsm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in x
         call dtimer(dtime,itime,-1)
         call cgpuwppfft2rcsxn(g_fn,g_bsm,isign,g_mixup,g_sct,indx,indy,&
     &ndim,kstrt,nvp,kxpd,kyp,nxhd,kypd,nxhyd,nxyhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPUPPFFT2RRCU(g_f,g_g,g_bsm,g_brm,bsm,brm,isign,ttp,   &
     &indx,indy,kstrt,nvp,kxpd,kyp,nxhd,nyd,kypd)
! wrapper function for gpu-mpi parallel real to complex fft,
! based on NVIDIA FFT, without packed data
! if isign = -1, g_f = input, g_g = output
! if isign = 1, g_g = input, g_f = output
! (g_bsm,g_brm)/(bsm,brm) are temporary scratch arrays on GPU/host
! nxhd must be = nx/2 + 1
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kxpd, kyp
      integer, intent(in) :: nxhd, nyd, kypd
      type (c_ptr) :: g_f, g_g, g_bsm, g_brm
!     complex dimension g_f(nxhd,kypd), g_g(nyd,kxpd)
!     complex dimension g_bsm(kxpd*kyp,nvp), g_brm(kxpd*kyp,nvp)
      real, dimension(2) :: ttp
      complex, dimension(kxpd*kyp,nvp), target :: bsm, brm
! local data
! kasync = (0,1) = (no,yes) use asynchronous communications
      integer, parameter :: kasync = 1
      integer :: nx, nxh1, ny, nxyp
      real :: ani, time
      double precision :: dtime
      type (timeval) :: itime
      nx = 2**indx
      nxh1 = nx/2 + 1
      ny = 2**indy
      nxyp = kxpd*kyp*(nvp-1)
      ani = 1.0/(real(nx)*real(ny))
! inverse fourier transform
      if (isign < 0) then
! first transpose in x
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcux(g_f,g_bsm,isign,indx,indy,kstrt,nvp,kxpd,kyp&
     &,nxhd,kypd)
! transpose on local GPU with scaling
         call cgpuppsltpose(g_f,g_g,ani,nxhd,ny,kxpd,kyp,kstrt,nxhd,nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSE(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt,  &
     &nvp)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(bsm(1,1)),g_bsm,nxyp)
               call PPPTPOSE(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp)
               call gpu_ccopyin(c_loc(brm(1,1)),g_brm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in y
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuy(g_g,g_brm,isign,indx,indy,kstrt,nvp,kxpd,kyp&
     &,nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! forward fourier transform
      else if (isign > 0) then
! first transpose in y
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuy(g_g,g_brm,isign,indx,indy,kstrt,nvp,kxpd,kyp&
     &,nyd)
! transpose on local GPU
         call cgpuppltpose(g_g,g_f,ny,nxhd,kyp,kxpd,kstrt,nyd,nxhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
        if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSE(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt,  &
     &nvp)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(brm(1,1)),g_brm,nxyp)
               call PPPTPOSE(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp)
               call gpu_ccopyin(c_loc(bsm(1,1)),g_bsm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in x
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcux(g_f,g_bsm,isign,indx,indy,kstrt,nvp,kxpd,kyp&
     &,nxhd,kypd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPUPPFFT2RRCUN(g_fn,g_gn,g_bsm,g_brm,bsm,brm,isign,    &
     &ttp,indx,indy,kstrt,nvp,ndim,kxpd,kyp,nxhd,nyd,kypd)
! wrapper function for multiple gpu-mpi parallel real to complex ffts,
! based on NVIDIA FFT, without packed data
! if isign = -1, g_fn = input, g_gn = output
! if isign = 1, g_gn = input, g_fn = output
! (g_bsm,g_brm)/(bsm,brm) are temporary scratch arrays on GPU/host
! ndim = vector dimension
! nxhd must be = nx/2 + 1
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, ndim, kxpd
      integer, intent(in) :: kyp, nxhd, nyd, kypd
      type (c_ptr) :: g_fn, g_gn, g_bsm, g_brm
!     complex dimension g_fn(nxhd,ndim,kypd), g_gn(nyd,ndim,kxpd)
!     complex dimension g_bsm(kxpd*ndim*kyp,nvp)
!     complex dimension g_brm(kxpd*ndim*kyp,nvp)
      real, dimension(2) :: ttp
      complex, dimension(kxpd*ndim*kyp,nvp), target :: bsm, brm
! local data
! kasync = (0,1) = (no,yes) use asynchronous communications
      integer, parameter :: kasync = 1
      integer :: nx, nxh1, ny, nxyp
      real :: ani, time
      double precision :: dtime
      type (timeval) :: itime
      nx = 2**indx
      nxh1 = nx/2 + 1
      ny = 2**indy
      nxyp = kxpd*ndim*kyp*(nvp-1)
      ani = 1.0/(real(nx)*real(ny))
! inverse fourier transform
      if (isign < 0) then
! first transpose in x
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuxn(g_fn,g_bsm,isign,indx,indy,ndim,kstrt,nvp, &
     &kxpd,kyp,nxhd,kyp)
! transpose on local GPU with scaling
         call cgpuppsltposen(g_fn,g_gn,ani,nxhd,ny,kxpd,kyp,kstrt,ndim, &
     &nxhd,nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSEN(g_bsm,g_brm,bsm,brm,nx,ny,kxpd,kyp,kstrt, &
     &nvp,ndim)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(bsm(1,1)),g_bsm,nxyp)
               call PPPTPOSEN(bsm,brm,nxh1,ny,kxpd,kyp,kstrt,nvp,ndim)
               call gpu_ccopyin(c_loc(brm(1,1)),g_brm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in y
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuyn(g_gn,g_brm,isign,indx,indy,ndim,kstrt,nvp, &
     &kxpd,kyp,nyd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! forward fourier transform
      else if (isign > 0) then
! first transpose in y
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuyn(g_gn,g_brm,isign,indx,indy,ndim,kstrt,nvp, &
     &kxpd,kyp,nyd)
! transpose on local GPU
         call cgpuppltposen(g_gn,g_fn,ny,nxhd,kyp,kxpd,kstrt,ndim,nyd,  &
     &nxhd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
! transpose between GPUs
         if (nvp > 1) then
            call dtimer(dtime,itime,-1)
! use asynchronous communication
            if (kasync==1) then
               call GPPTPOSEN(g_brm,g_bsm,brm,bsm,ny,nx,kyp,kxpd,kstrt, &
     &nvp,ndim)
! use synchronous communication
            else
               call gpu_ccopyout(c_loc(brm(1,1)),g_brm,nxyp)
               call PPPTPOSEN(brm,bsm,ny,nxh1,kyp,kxpd,kstrt,nvp,ndim)
               call gpu_ccopyin(c_loc(bsm(1,1)),g_bsm,nxyp)
            endif
            call dtimer(dtime,itime,1)
            time = real(dtime)
            ttp(2) = ttp(2) + time
         endif
! then transpose in x
         call dtimer(dtime,itime,-1)
         call gpupfft2rrcuxn(g_fn,g_bsm,isign,indx,indy,ndim,kstrt,nvp, &
     &kxpd,kyp,nxhd,kypd)
         call dtimer(dtime,itime,1)
         time = real(dtime)
         ttp(1) = ttp(1) + time
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPPORDER2L(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,    &
     &g_ncl,g_ihole,g_ncll,g_nclr,sbufl,sbufr,rbufl,rbufr,ncll,nclr,mcll&
     &,mclr,ttp,noff,nyp,kstrt,nvp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,    &
     &npbmx,ntmax,nbmax,g_irc)
! this subroutine performs an mpi-gpu particle sort by x,y grid in tiles
! of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
      implicit none
      integer, intent(in) :: noff, nyp, kstrt, nvp, idimp, nppmx, nx, ny
      integer, intent(in) :: mx, my, mx1, myp1, npbmx, ntmax, nbmax
      type (c_ptr) :: g_ppart, g_ppbuff, g_sbufl, g_sbufr
      type (c_ptr) :: g_kpic, g_ncl, g_ihole, g_ncll, g_nclr, g_irc
!     real dimension g_ppart(nppmx,idimp,mx1*myp1)
!     real dimension g_ppbuff(npbmx,idimp,mx1*myp1)
!     real dimension g_sbufl(idimp,nbmax), g_sbufr(idimp,nbmax)
!     integer dimension g_kpic(mx1*myp1), g_ncl(8,mx1*myp1)
!     integer dimension g_ihole(2,ntmax+1,mx1*myp1)
!     integer dimension g_ncll(3,mx1), g_nclr(3,mx1)
      real, dimension(idimp,nbmax), target :: sbufl, sbufr, rbufl, rbufr
      real, dimension(2) :: ttp
      integer, dimension(3,mx1), target :: ncll, nclr, mcll, mclr
! local data
      real :: time
      double precision :: dtime
      type (timeval) :: itime
! first part of particle reorder on x and y cell with mx, my tiles
      call dtimer(dtime,itime,-1)
      call cgpupppord2la(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl, &
     &g_ihole,g_ncll,g_nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,  &
     &npbmx,ntmax,nbmax,g_irc)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      ttp(1) = ttp(1) + time
! move particles on GPU into appropriate spatial regions:
! updates rbufr, rbufl, mcll, mclr
      call dtimer(dtime,itime,-1)
      call gpu_icopyout(c_loc(ncll(1,1)),g_ncll,3*mx1)
      call gpu_icopyout(c_loc(nclr(1,1)),g_nclr,3*mx1)
      call gpu_fcopyout(c_loc(sbufl(1,1)),g_sbufl,idimp*ncll(3,mx1))
      call gpu_fcopyout(c_loc(sbufr(1,1)),g_sbufr,idimp*nclr(3,mx1))
      call PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,  &
     &nvp,idimp,nbmax,mx1)
      call gpu_icopyin(c_loc(mcll(1,1)),g_ncll,3*mx1)
      call gpu_icopyin(c_loc(mclr(1,1)),g_nclr,3*mx1)
      call gpu_fcopyin(c_loc(rbufl(1,1)),g_sbufl,idimp*mcll(3,mx1))
      call gpu_fcopyin(c_loc(rbufr(1,1)),g_sbufr,idimp*mclr(3,mx1))
      call dtimer(dtime,itime,1)
      time = real(dtime)
      ttp(2) = ttp(2) + time
! second part of particle reorder on x and y cell with mx, my tiles:
! updates g_ppart, g_kpic, g_irc
      call dtimer(dtime,itime,-1)
      call cgpupppord2lb(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl, &
     &g_ihole,g_ncll,g_nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,g_irc&
     &)
      call dtimer(dtime,itime,1)
      time = real(dtime)
      ttp(1) = ttp(1) + time
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPPTPOSE(g_bsm,g_btm,sm,tm,nx,ny,kxp,kyp,kstrt,nvp)
! this subroutine sends and receives data between GPUS on different MPI
! nodes to perform a transpose of a matrix distributed in y, to another
! matrix distributed in x.
! one message is sent and received at a time.
! data from GPU is sent asynchronous, overlapping with MPI
! g_bsm/g_btm are complex buffers on GPU to be sent/received
! sm/tm = complex buffers on host to be sent/received
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp
      type (c_ptr) :: g_bsm, g_btm
!     complex dimension g_bsm(kxp*kyp,nvp), g_btm(kxp*kyp,nvp)
      complex, dimension(kxp*kyp,nvp), target :: sm, tm
! local data
      integer :: j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp
      ks = kstrt - 1
      kyps = min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*kyp
! special case for one processor
! better to use a kernel function
      if (nvp==1) then
         call gpu_ccopyout(c_loc(sm(1,1)),g_bsm,kxyp)
         do j = 1, kxyp
            tm(j,1) = sm(j,1)
         enddo
         call gpu_ccopyin(c_loc(tm(1,1)),g_btm,kxyp)
         return
      endif
      nn = 0
      ns = kxyp
      stp = 1
! send first group to host from GPU
      call gpu_cascopyout(c_loc(sm(1,1)),g_bsm,0,ns,stp)
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
         if (id /= ks) then
! adjust counter to omit data sent to oneself
            nn = nn + 1
! send next group to host from GPU
            if ((nn+1) < nvp) then
               st = stp + 1
               if (st > 2) st = st - 2 
               call gpu_cascopyout(c_loc(sm(1,nn+1)),g_bsm,ns*nn,ns,st)
            endif
! wait for previous MPI sends and receives to complete
            if (nn > 1) then
! call MPI_WAIT(msid,istatus,ierr); call MPI_WAIT(mrid,istatus,ierr)
               call ACSNDREC(sm,0,0,0,3)
! copy received group from host to GPU
               call gpu_cascopyin(c_loc(tm(1,nn-1)),g_btm,ns*(nn-2),ns,3&
     &)
            endif
! calculate length of data to send
            joff = kxp*id
            ld = kyps*min(kxp,max(0,nx-joff))
! post receive
! call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,mrid,ierr)
            call ACSNDREC(tm(1,nn),id,kxyp,n,1)
! wait for previous group to arrive from GPU
            call gpu_waitstream(stp)
            stp = st
! send data
! call MPI_ISEND(sm(1,nn),ld,mcplx,id,n,lgrp,msid,ierr)
            call ACSNDREC(sm(1,nn),id,ld,n,2)
         endif
      enddo
! wait for last sends and receives to complete
! call MPI_WAIT(msid,istatus,ierr); call MPI_WAIT(mrid,istatus,ierr)
      call ACSNDREC(sm,0,0,0,3)
      call gpu_cascopyin(c_loc(tm(1,nn)),g_btm,ns*(nn-1),ns,3)
! wait for last group item to arrive
      call gpu_waitstream(3)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine GPPTPOSEN(g_bsm,g_btm,sm,tm,nx,ny,kxp,kyp,kstrt,nvp,   &
     &ndim)
! this subroutine sends and receives data between GPUS on different MPI
! nodes to perform a transpose of an n component  matrix distributed in
! y, to another an n component matrix distributed in x.
! one message is sent and received at a time.
! data from GPU is sent asynchronous, overlapping with MPI.
! g_bsm/g_btm are complex buffers on GPU to be sent/received
! sm/tm = complex buffers on host to be sent/received
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      type (c_ptr) :: g_bsm, g_btm
!     complex dimension g_bsm(kxp*ndim*kyp,nvp), g_btm(kxp*ndim*kyp,nvp)
      complex, dimension(kxp*ndim*kyp,nvp), target :: sm, tm
! local data
      integer :: j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp
      ks = kstrt - 1
      kyps = ndim*min(kyp,max(0,ny-kyp*ks))
      kxyp = kxp*ndim*kyp
! special case for one processor
! better to use a kernel function
      if (nvp==1) then
         call gpu_ccopyout(c_loc(sm(1,1)),g_bsm,kxyp)
         do j = 1, kxyp
            tm(j,1) = sm(j,1)
          enddo
         call gpu_ccopyin(c_loc(tm(1,1)),g_btm,kxyp)
         return
      endif
      nn = 0
      ns = kxyp
      stp = 1
! send first group to host from GPU
      call gpu_cascopyout(c_loc(sm(1,1)),g_bsm,0,ns,stp)
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
         if (id /= ks) then
! adjust counter to omit data sent to oneself
            nn = nn + 1
! send next group to host from GPU
            if ((nn+1) < nvp) then
               st = stp + 1
               if (st > 2) st = st - 2
               call gpu_cascopyout(c_loc(sm(1,nn+1)),g_bsm,ns*nn,ns,st)
            endif
! wait for previous MPI sends and receives to complete
            if (nn > 1) then
! call MPI_WAIT(msid,istatus,ierr); call MPI_WAIT(mrid,istatus,ierr)
               call ACSNDREC(sm,0,0,0,3)
! copy received group from host to GPU
               call gpu_cascopyin(c_loc(tm(1,nn-1)),g_btm,ns*(nn-2),ns,3&
     &)
            endif
! calculate length of data to send
            joff = kxp*id
            ld = kyps*min(kxp,max(0,nx-joff))
! post receive
! call MPI_IRECV(tm(1,nn),kxyp,mcplx,id,n,lgrp,mrid,ierr)
            call ACSNDREC(tm(1,nn),id,kxyp,n,1)
! wait for previous group to arrive from GPU
            call gpu_waitstream(stp)
            stp = st
! send data
! call MPI_ISEND(sm(1,nn),ld,mcplx,id,n,lgrp,msid,ierr)
            call ACSNDREC(sm(1,nn),id,ld,n,2)
         endif
      enddo
! wait for last sends and receives to complete
! call MPI_WAIT(msid,istatus,ierr); call MPI_WAIT(mrid,istatus,ierr)
      call ACSNDREC(sm,0,0,0,3)
      call gpu_cascopyin(c_loc(tm(1,nn)),g_btm,ns*(nn-1),ns,3)
! wait for last group item to arrive
      call gpu_waitstream(3)
      end subroutine
!
      end module
