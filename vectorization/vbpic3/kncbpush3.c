/* KNC C Library for Skeleton 3D Electromagnetic Vector PIC Code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include "kncbpush3.h"

/*--------------------------------------------------------------------*/
void ckncxiscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method, exclusive scan. */
/* requires KNC, isdata needs to be 64 byte aligned */
/* local data */
   int j, ns, isum, ist;
   __m512i v_m1, v_m2, v_it0, v_it, v_is, v_ioff;
   __attribute__((aligned(64))) int ll[16];
   ns = 16*(nths/16);
   v_m1 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
   v_m2 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
   isum = 0;
   v_ioff = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
   for (j = 0; j < ns; j+=16) {
/* load data */
      v_it0 = _mm512_load_epi32(&isdata[j]);
/* first pass */
      v_is = _mm512_shuffle_epi32(v_it0,177);
      v_it = _mm512_mask_add_epi32(v_it0,_mm512_int2mask(43690),v_it0,
             v_is);
/* second pass */
      v_is = _mm512_shuffle_epi32(v_it,80);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(52428),v_it,
             v_is);
/* third pass */
      v_is = _mm512_permutevar_epi32(v_m1,v_it);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(61680),v_it,
             v_is);
/* fourth pass */
      v_is = _mm512_permutevar_epi32(v_m2,v_it);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(65280),v_it,
             v_is);
/* add offset */
      v_it = _mm512_add_epi32(v_it,v_ioff);
/* next offset */
      v_ioff = _mm512_shuffle_epi32(v_it,255);
      v_ioff = _mm512_permute4f128_epi32(v_ioff,255);
/* subtract for exclusive scan */
      v_it = _mm512_sub_epi32(v_it,v_it0);
/* write data */
      _mm512_store_epi32(&isdata[j],v_it);
   }
   if (ns > 0) {
      _mm512_store_epi32(ll,v_ioff);
      isum = ll[0];
   }
/* loop over remaining elements */
   for (j = ns; j < nths; j++) {
      ist = isdata[j];
      isdata[j] = isum;
      isum += ist;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncgbpush3lt(float part[], float fxyz[], float bxyz[], float qbm,
                   float dt, float dtc, float *ek, int idimp, int nop,
                   int npe, int nx, int ny, int nz, int nxv, int nyv,
                   int nzv, int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with magnetic field.  Using the Boris Mover.
   scalar version using guard cells
   190 flops/particle, 1 divide, 54 loads, 6 stores
   input: all, output: part, ek
   velocity equations used are:
   vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
   omz = (q/m)*bz(x(t),y(t),z(t)).
   position equations used are:
   x(t+dt)=x(t) + vx(t+dt/2)*dt
   y(t+dt)=y(t) + vy(t+dt/2)*dt
   z(t+dt)=z(t) + vz(t+dt/2)*dt
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   part[3][n] = velocity vx of particle n
   part[4][n] = velocity vy of particle n
   part[5][n] = velocity vz of particle n
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   dtc = time interval between successive co-ordinate calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
        (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
        .25*(vz(t+dt/2) + vz(t-dt/2))**2)
   idimp = size of phase space = 6
   nop = number of particles
   npe = first dimension of particle array
   nx/ny/nz = system length in x/y/z direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
   int j, nps, nn, mm, ll, nm, nxyv;
   float qtmh, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9;
   float x, y, z, vx, vy, vz;
   double sum1, d0;
   __m512i v_nxv4, v_nxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_perm;
   __m512 v_qtmh, v_dt, v_dtc, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   qtmh = 0.5f*qbm*dt;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   v_nxv4 = _mm512_set1_epi32(4*nxv);
   v_nxyv4 = _mm512_set1_epi32(4*nxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_sum1 = _mm512_set1_pd(0.0);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = x - (float) nn; */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_sub_ps(v_x,v_dxp);
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nm = 4*(nn + nxv*mm + nxyv*ll); */
      v_it = _mm512_mullo_epi32(v_nxyv4,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv4,v_mm));
      v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_one,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*    nn = nm; */
      _mm512_store_epi32(kk,v_nm);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*    dx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_dx = _mm512_mul_ps(v_amx,a);
      v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*    dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_dy = _mm512_mul_ps(v_amx,b);
      v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*    dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_dz = _mm512_mul_ps(v_amx,c);
      v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*    mm = nn + 4*nxv;                                     */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxv;                                    */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxv;                                      */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*    dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);          */
      v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
      v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
      v_dx = _mm512_mul_ps(v_amz,v_dx);
/*    dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);      */
      v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
      v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
      v_dy = _mm512_mul_ps(v_amz,v_dy);
/*    dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);      */
      v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
      v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
      v_dz = _mm512_mul_ps(v_amz,v_dz);
/*    nn += 4*nxyv;                                           */
      v_nn = _mm512_add_epi32(v_nm,v_nxyv4);
      _mm512_store_epi32(kk,v_nn);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
/* where nn = nn + 4*nxyv;                                   */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*nxyv;                                     */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*    vx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_vx = _mm512_mul_ps(v_amx,a);
      v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*    vy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_vy = _mm512_mul_ps(v_amx,b);
      v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*    vz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_vz = _mm512_mul_ps(v_amx,c);
      v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*    mm = nn + 4*nxv;                                        */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*    dx = dx + dzp*(vx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);     */
      v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
      v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
      v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*    dy = dy + dzp*(vy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]); */
      v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
      v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
      v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*    dz = dz + dzp*(vz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]); */
      v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
      v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
      v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*    nn = nm; */
      _mm512_store_epi32(kk,v_nm);
/* load bxyz[nn:nn+3] and bxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[nn:nn+3] field components */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[nn+4:nn+7] field components */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*    ox = amx*bxyz[nn] + amy*bxyz[nn+4];                     */
      v_ox = _mm512_mul_ps(v_amx,a);
      v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*    oy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];                 */
      v_oy = _mm512_mul_ps(v_amx,b);
      v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*    oz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];                 */
      v_oz = _mm512_mul_ps(v_amx,c);
      v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*    mm = nn + 4*nxv;                                     */
/* load bxyz[mm:mm+3] and bxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxv;                                    */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxv;                                      */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*    ox = amz*(ox + dyp*bxyz[mm] + dx1*bxyz[mm+4]);          */
      v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
      v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
      v_ox = _mm512_mul_ps(v_amz,v_ox);
/*    oy = amz*(oy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);      */
      v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
      v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
      v_oy = _mm512_mul_ps(v_amz,v_oy);
/*    oz = amz*(oz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);      */
      v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
      v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
      v_oz = _mm512_mul_ps(v_amz,v_oz);
/*    nn += 4*nxyv;                                           */
      v_nn = _mm512_add_epi32(v_nm,v_nxyv4);
      _mm512_store_epi32(kk,v_nn);
/* load bxyz[nn:nn+3] and bxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[nn:nn+3] field components */
/* where nn = nn + 4*nxyv;                                   */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*nxyv;                                     */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*    vx = amx*bxyz[nn] + amy*bxyz[nn+4];                     */
      v_vx = _mm512_mul_ps(v_amx,a);
      v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*    vy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];                 */
      v_vy = _mm512_mul_ps(v_amx,b);
      v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*    vz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];                 */
      v_vz = _mm512_mul_ps(v_amx,c);
      v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*    mm = nn + 4*nxv;                                        */
/* load bxyz[mm:mm+3] and bxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*    ox = ox + dzp*(vx + dyp*bxyz[mm] + dx1*bxyz[mm+4]);     */
      v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
      v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
      v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*    oy = oy + dzp*(vy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]); */
      v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
      v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
      v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*    oz = oz + dzp*(vz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]); */
      v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
      v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
      v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*    dx *= qtmh; */
/*    dy *= qtmh; */
/*    dz *= qtmh; */
      v_dx = _mm512_mul_ps(v_dx,v_qtmh);
      v_dy = _mm512_mul_ps(v_dy,v_qtmh);
      v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*    acx = part[j+3*npe] + dx; */
/*    acy = part[j+4*npe] + dy; */
/*    acz = part[j+5*npe] + dz; */
      a = _mm512_add_ps(v_dx,_mm512_load_ps(&part[j+3*npe]));
      b = _mm512_add_ps(v_dy,_mm512_load_ps(&part[j+4*npe]));
      c = _mm512_add_ps(v_dz,_mm512_load_ps(&part[j+5*npe]));
/* time-centered kinetic energy */
/*    sum1 += (acx*acx + acy*acy + acz*acz); */
      v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
      v_at = _mm512_fmadd_ps(c,c,v_at);
/* convert to double precision before accumulating */
      v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
      v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* calculate cyclotron frequency */
/*    omxt = qtmh*ox; */
/*    omyt = qtmh*oy; */
/*    omzt = qtmh*oz; */
      e = _mm512_mul_ps(v_qtmh,v_ox);
      f = _mm512_mul_ps(v_qtmh,v_oy);
      g = _mm512_mul_ps(v_qtmh,v_oz);
/* calculate rotation matrix */
/*    vx = omxt*omxt; */
      v_vx = _mm512_mul_ps(e,e);
/*    vy = omyt*omyt; */
      v_vy = _mm512_mul_ps(f,f);
/*    vz = omzt*omzt; */
      v_vz = _mm512_mul_ps(g,g);
/*    omt = omxt*omxt + omyt*omyt + omzt*omzt; */
      v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*    anorm = 2.0f/(1.0f + omt); */
      d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*    omt = 0.5f*(1.0f - omt); */
      h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*    vx = (omt + vx)*acx; */
      v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*    vy = (omt + vy)*acy; */
      v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*    vz = (omt + vz)*acz; */
      v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*    omt = omxt*omyt; */
      h = _mm512_mul_ps(e,f);
/*    vx = vx + (omzt + omt)*acy; */
      v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*    vy = vy + (omt - omzt)*acx; */
      v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*    omt = omxt*omzt;  */
      h = _mm512_mul_ps(e,g);
/*    vx = vx + (omt - omyt)*acz; */
      v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*    vz = vz + (omt + omyt)*acx; */
      v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*    omt = omyt*omzt; */
      h = _mm512_mul_ps(f,g);
/*    vy = vy + (omt + omxt)*acz; */
      v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*    vz = vz + (omt - omxt)*acy; */
      v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new velocity */
/*    vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*    vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*    vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
      v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
      v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
      v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* new position */
/*    dx = x + vx*dtc; */
/*    dy = y + vy*dtc; */
/*    dz = z + vz*dtc; */
      v_dx = _mm512_fmadd_ps(v_vx,v_dtc,v_x);
      v_dy = _mm512_fmadd_ps(v_vy,v_dtc,v_y);
      v_dz = _mm512_fmadd_ps(v_vz,v_dtc,v_z);
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         v_dx = _mm512_mask_add_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         v_dx = _mm512_mask_sub_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dy < edgely) dy += edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         v_dy = _mm512_mask_add_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dy >= edgery) dy -= edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         v_dy = _mm512_mask_sub_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if ((dz < edgelz) || (dz >= edgerz)) { */
/*          dz = z;                             */
/*          vz = -vz;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
               _MM_CMPINT_GE));
         v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
         v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* set new position */
/*    part[j] = dx;       */
/*    part[j+npe] = dy;   */
/*    part[j+2*npe] = dz; */
      _mm512_store_ps(&part[j],v_dx);
      _mm512_store_ps(&part[j+npe],v_dy);
      _mm512_store_ps(&part[j+2*npe],v_dz);
/* set new velocity */
/*    part[j+3*npe] = vx; */
/*    part[j+4*npe] = vy; */
/*    part[j+5*npe] = vz; */
      _mm512_store_ps(&part[j+3*npe],v_vx);
      _mm512_store_ps(&part[j+4*npe],v_vy);
      _mm512_store_ps(&part[j+5*npe],v_vz);
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nm = 4*(nn + nxv*mm + nxyv*ll);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amz = 1.0f - dzp;
      amy = dxp*amy;
/* find electric field */
      nn = nm;
      dx = amx*fxyz[nn] + amy*fxyz[nn+4];
      dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
      nn += 4*nxyv;
      acx = amx*fxyz[nn] + amy*fxyz[nn+4];
      acy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      acz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = dx + dzp*(acx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = dy + dzp*(acy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = dz + dzp*(acz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
/* find magnetic field */
      nn = nm;
      ox = amx*bxyz[nn] + amy*bxyz[nn+4];
      oy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];
      oz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];
      mm = nn + 4*nxv;
      ox = amz*(ox + dyp*bxyz[mm] + dx1*bxyz[mm+4]);
      oy = amz*(oy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);
      oz = amz*(oz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);
      nn += 4*nxyv;
      acx = amx*bxyz[nn] + amy*bxyz[nn+4];
      acy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];
      acz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];
      mm = nn + 4*nxv;
      ox = ox + dzp*(acx + dyp*bxyz[mm] + dx1*bxyz[mm+4]);
      oy = oy + dzp*(acy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);
      oz = oz + dzp*(acz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[j+3*npe] + dx;
      acy = part[j+4*npe] + dy;
      acz = part[j+5*npe] + dz;
/* time-centered kinetic energy */
      sum1 += (acx*acx + acy*acy + acz*acz);
/* calculate cyclotron frequency */
      omxt = qtmh*ox;
      omyt = qtmh*oy;
      omzt = qtmh*oz;
/* calculate rotation matrix */
      omt = omxt*omxt + omyt*omyt + omzt*omzt;
      anorm = 2.0f/(1.0f + omt);
      omt = 0.5f*(1.0f - omt);
      rot4 = omxt*omyt;
      rot7 = omxt*omzt;
      rot8 = omyt*omzt;
      rot1 = omt + omxt*omxt;
      rot5 = omt + omyt*omyt;
      rot9 = omt + omzt*omzt;
      rot2 = omzt + rot4;
      rot4 -= omzt;
      rot3 = -omyt + rot7;
      rot7 += omyt;
      rot6 = omxt + rot8;
      rot8 -= omxt;
/* new velocity */
      vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm;
      vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm;
      vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* new position */
      dx = x + vx*dtc;
      dy = y + vy*dtc;
      dz = z + vz*dtc;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            vy = -vy;
         }
         if ((dz < edgelz) || (dz >= edgerz)) {
            dz = z;
            vz = -vz;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            vy = -vy;
         }
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
      part[j+2*npe] = dz;
/* set new velocity */
      part[j+3*npe] = vx;
      part[j+4*npe] = vy;
      part[j+5*npe] = vz;
   }
/* normalize kinetic energy */
/* *ek += 0.5f*sum1; */
   d0 = _mm512_reduce_add_pd(v_sum1);
   *ek += 0.5f*(sum1 + d0);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrbpush3lt(float part[], float fxyz[], float bxyz[], float qbm,
                    float dt, float dtc, float ci, float *ek, int idimp,
                    int nop, int npe, int nx, int ny, int nz, int nxv,
                    int nyv, int nzv, int ipbc) {
/* for 3d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, for relativistic particles with magnetic field
   Using the Boris Mover.
   vector version using guard cells
   202 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
   input: all, output: part, ek
   momentum equations used are:
   px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fx(x(t),y(t),z(t))*dt)
   py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fy(x(t),y(t),z(t))*dt)
   pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
      rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
      rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
      .5*(q/m)*fz(x(t),y(t),z(t))*dt)
   where q/m is charge/mass, and the rotation matrix is given by:
      rot[0] = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[1] = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[2] = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[3] = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
      rot[4] = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
      rot[5] = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[6] = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[7] = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
      rot[8] = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
   and om**2 = omx**2 + omy**2 + omz**2
   the rotation matrix is determined by:
   omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
   omy = (q/m)*by(x(t),y(t),z(t))*gami,
   omz = (q/m)*bz(x(t),y(t),z(t))*gami,
   where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
   position equations used are:
   x(t+dt) = x(t) + px(t+dt/2)*dtg
   y(t+dt) = y(t) + py(t+dt/2)*dtg
   z(t+dt) = z(t) + pz(t+dt/2)*dtg
   where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
   pz(t+dt/2)*pz(t+dt/2))*ci*ci)
   fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
   bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
   are approximated by interpolation from the nearest grid points:
   fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
                  + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
             + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
                  + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   part[3][n] = momentum px of particle n
   part[4][n] = momentum py of particle n
   part[5][n] = momentum pz of particle n
   fxyz[l][k][j][0] = x component of force/charge at grid (j,k,l)
   fxyz[l][k][j][1] = y component of force/charge at grid (j,k,l)
   fxyz[l][k][j][2] = z component of force/charge at grid (j,k,l)
   that is, convolution of electric field over particle shape
   bxyz[l][k][j][0] = x component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][1] = y component of magnetic field at grid (j,k,l)
   bxyz[l][k][j][2] = z component of magnetic field at grid (j,k,l)
   that is, the convolution of magnetic field over particle shape
   qbm = particle charge/mass ratio
   dt = time interval between successive force calculations
   dtc = time interval between successive co-ordinate calculations
   ci = reciprocal of velocity of light
   kinetic energy/mass at time t is also calculated, using
   ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
       (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
        (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
   idimp = size of phase space = 6
   nop = number of particles
   npe = first dimension of particle array
   nx/ny/nz = system length in x/y/z direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   nzv = fourth dimension of field array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
   fxyz needs to have 4 components, although one is not used
local data                                                            */
   int j, nps, nn, mm, ll, nm, nxyv;
   float qtmh, ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, ox, oy, oz, dx1;
   float acx, acy, acz, p2, gami, qtmg, omxt, omyt, omzt, omt, anorm;
   float rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, dtg;
   float x, y, z, vx, vy, vz;
   double sum1, d0;
   __m512i v_nxv4, v_nxyv4;
   __m512i v_nn, v_mm, v_ll, v_nm, v_it, v_perm;
   __m512 v_qtmh, v_ci2, v_dt, v_dtc, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s;
   __m512 v_two, v_half, v_ox, v_oy, v_oz;
   __m512d v_sum1, v_d;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   qtmh = 0.5f*qbm*dt;
   ci2 = ci*ci;
   sum1 = 0.0;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   v_nxv4 = _mm512_set1_epi32(4*nxv);
   v_nxyv4 = _mm512_set1_epi32(4*nxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qtmh = _mm512_set1_ps(qtmh);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_dtc = _mm512_set1_ps(dtc);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_two = _mm512_set1_ps(2.0f);
   v_half = _mm512_set1_ps(0.5f);
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
   v_sum1 = _mm512_set1_pd(0.0);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = x - (float) nn; */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_sub_ps(v_x,v_dxp);
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nm = 4*(nn + nxv*mm + nxyv*ll); */
      v_it = _mm512_mullo_epi32(v_nxyv4,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv4,v_mm));
      v_nm = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_one,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/* find electric field */
/*    nn = nm; */
      _mm512_store_epi32(kk,v_nm);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of electric field */
/*    dx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_dx = _mm512_mul_ps(v_amx,a);
      v_dx = _mm512_fmadd_ps(v_amy,p,v_dx);
/*    dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_dy = _mm512_mul_ps(v_amx,b);
      v_dy = _mm512_fmadd_ps(v_amy,q,v_dy);
/*    dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_dz = _mm512_mul_ps(v_amx,c);
      v_dz = _mm512_fmadd_ps(v_amy,r,v_dz);
/*    mm = nn + 4*nxv;                                     */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxv;                                    */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxv;                                      */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of electric field */
/*    dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);          */
      v_dx = _mm512_fmadd_ps(v_dyp,a,v_dx);
      v_dx = _mm512_fmadd_ps(v_dx1,p,v_dx);
      v_dx = _mm512_mul_ps(v_amz,v_dx);
/*    dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);      */
      v_dy = _mm512_fmadd_ps(v_dyp,b,v_dy);
      v_dy = _mm512_fmadd_ps(v_dx1,q,v_dy);
      v_dy = _mm512_mul_ps(v_amz,v_dy);
/*    dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);      */
      v_dz = _mm512_fmadd_ps(v_dyp,c,v_dz);
      v_dz = _mm512_fmadd_ps(v_dx1,r,v_dz);
      v_dz = _mm512_mul_ps(v_amz,v_dz);
/*    nn += 4*nxyv;                                           */
      v_nn = _mm512_add_epi32(v_nm,v_nxyv4);
      _mm512_store_epi32(kk,v_nn);
/* load fxyz[nn:nn+3] and fxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[nn:nn+3] field components */
/* where nn = nn + 4*nxyv;                                   */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*nxyv;                                     */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of electric field */
/*    vx = amx*fxyz[nn] + amy*fxyz[nn+4];                     */
      v_vx = _mm512_mul_ps(v_amx,a);
      v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*    vy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];                 */
      v_vy = _mm512_mul_ps(v_amx,b);
      v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*    vz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];                 */
      v_vz = _mm512_mul_ps(v_amx,c);
      v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*    mm = nn + 4*nxv;                                        */
/* load fxyz[mm:mm+3] and fxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &fxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &fxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &fxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &fxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for fxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for fxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of electric field */
/*    dx = dx + dzp*(vx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);     */
      v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
      v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
      v_dx = _mm512_fmadd_ps(v_dzp,v_vx,v_dx);
/*    dy = dy + dzp*(vy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]); */
      v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
      v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
      v_dy = _mm512_fmadd_ps(v_dzp,v_vy,v_dy);
/*    dz = dz + dzp*(vz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]); */
      v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
      v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
      v_dz = _mm512_fmadd_ps(v_dzp,v_vz,v_dz);
/* find magnetic field */
/*    nn = nm; */
      _mm512_store_epi32(kk,v_nm);
/* load bxyz[nn:nn+3] and bxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[nn:nn+3] field components */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[nn+4:nn+7] field components */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find first part of magnetic field */
/*    ox = amx*bxyz[nn] + amy*bxyz[nn+4];                     */
      v_ox = _mm512_mul_ps(v_amx,a);
      v_ox = _mm512_fmadd_ps(v_amy,p,v_ox);
/*    oy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];                 */
      v_oy = _mm512_mul_ps(v_amx,b);
      v_oy = _mm512_fmadd_ps(v_amy,q,v_oy);
/*    oz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];                 */
      v_oz = _mm512_mul_ps(v_amx,c);
      v_oz = _mm512_fmadd_ps(v_amy,r,v_oz);
/*    mm = nn + 4*nxv;                                     */
/* load bxyz[mm:mm+3] and bxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxv;                                    */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxv;                                      */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find second part of magnetic field */
/*    ox = amz*(ox + dyp*bxyz[mm] + dx1*bxyz[mm+4]);          */
      v_ox = _mm512_fmadd_ps(v_dyp,a,v_ox);
      v_ox = _mm512_fmadd_ps(v_dx1,p,v_ox);
      v_ox = _mm512_mul_ps(v_amz,v_ox);
/*    oy = amz*(oy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);      */
      v_oy = _mm512_fmadd_ps(v_dyp,b,v_oy);
      v_oy = _mm512_fmadd_ps(v_dx1,q,v_oy);
      v_oy = _mm512_mul_ps(v_amz,v_oy);
/*    oz = amz*(oz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);      */
      v_oz = _mm512_fmadd_ps(v_dyp,c,v_oz);
      v_oz = _mm512_fmadd_ps(v_dx1,r,v_oz);
      v_oz = _mm512_mul_ps(v_amz,v_oz);
/*    nn += 4*nxyv;                                           */
      v_nn = _mm512_add_epi32(v_nm,v_nxyv4);
      _mm512_store_epi32(kk,v_nn);
/* load bxyz[nn:nn+3] and bxyz[nn+4:nn+7] field components */
/* first block of 4 particles */
      mm = kk[0];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14];
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15];
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[nn:nn+3] field components */
/* where nn = nn + 4*nxyv;                                   */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[nn+4:nn+7] field components */
/* where nn = nn + 4*nxyv;                                     */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find third part of magnetic field */
/*    vx = amx*bxyz[nn] + amy*bxyz[nn+4];                     */
      v_vx = _mm512_mul_ps(v_amx,a);
      v_vx = _mm512_fmadd_ps(v_amy,p,v_vx);
/*    vy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];                 */
      v_vy = _mm512_mul_ps(v_amx,b);
      v_vy = _mm512_fmadd_ps(v_amy,q,v_vy);
/*    vz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];                 */
      v_vz = _mm512_mul_ps(v_amx,c);
      v_vz = _mm512_fmadd_ps(v_amy,r,v_vz);
/*    mm = nn + 4*nxv;                                        */
/* load bxyz[mm:mm+3] and bxyz[mm+4:mm+7] field components */
/* first block of 4 particles */
      mm = kk[0] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[1] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[2] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[3] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      p = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* second block of 4 particles */
      mm = kk[4] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[5] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[6] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[7] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      b = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* third block of 4 particles */
      mm = kk[8] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[9] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[10] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[11] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      c = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      r = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* fourth block of 4 particles */
      mm = kk[12] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(255),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[13] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(255),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(255),
          &bxyz[mm+16]);
      mm = kk[14] + 4*nxv;
      e = _mm512_mask_loadunpacklo_ps(e,_mm512_int2mask(65280),
          &bxyz[mm]);
      e = _mm512_mask_loadunpackhi_ps(e,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      mm = kk[15] + 4*nxv;
      f = _mm512_mask_loadunpacklo_ps(f,_mm512_int2mask(65280),
          &bxyz[mm]);
      f = _mm512_mask_loadunpackhi_ps(f,_mm512_int2mask(65280),
          &bxyz[mm+16]);
      d = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),f,177);
      s = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(3855),e,177);
/* perform 16x3 transpose for bxyz[mm:mm+3] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(61680),b,177);
      f = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(3855),a,177);
      g = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(61680),d,177);
      b = _mm512_mask_permute4f128_ps(d,_mm512_int2mask(3855),c,177);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      c = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      b = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),b,78);
/* perform 16x3 transpose for bxyz[mm+4:mm+7] field components */
/* where mm = nn + 4*nxyv + 4*nxv;                           */
      p = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)p);
      q = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)q);
      r = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)r);
      s = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)s);
      e = _mm512_mask_permute4f128_ps(p,_mm512_int2mask(61680),q,177);
      f = _mm512_mask_permute4f128_ps(q,_mm512_int2mask(3855),p,177);
      g = _mm512_mask_permute4f128_ps(r,_mm512_int2mask(61680),s,177);
      q = _mm512_mask_permute4f128_ps(s,_mm512_int2mask(3855),r,177);
      p = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(65280),g,78);
      r = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(255),e,78);
      q = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(65280),q,78);
/* find fourth part of magnetic field */
/*    ox = ox + dzp*(vx + dyp*bxyz[mm] + dx1*bxyz[mm+4]);     */
      v_vx = _mm512_fmadd_ps(v_dyp,a,v_vx);
      v_vx = _mm512_fmadd_ps(v_dx1,p,v_vx);
      v_ox = _mm512_fmadd_ps(v_dzp,v_vx,v_ox);
/*    oy = oy + dzp*(vy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]); */
      v_vy = _mm512_fmadd_ps(v_dyp,b,v_vy);
      v_vy = _mm512_fmadd_ps(v_dx1,q,v_vy);
      v_oy = _mm512_fmadd_ps(v_dzp,v_vy,v_oy);
/*    oz = oz + dzp*(vz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]); */
      v_vz = _mm512_fmadd_ps(v_dyp,c,v_vz);
      v_vz = _mm512_fmadd_ps(v_dx1,r,v_vz);
      v_oz = _mm512_fmadd_ps(v_dzp,v_vz,v_oz);
/* calculate half impulse */
/*    dx *= qtmh; */
/*    dy *= qtmh; */
/*    dz *= qtmh; */
      v_dx = _mm512_mul_ps(v_dx,v_qtmh);
      v_dy = _mm512_mul_ps(v_dy,v_qtmh);
      v_dz = _mm512_mul_ps(v_dz,v_qtmh);
/* half acceleration */
/*    acx = part[j+3*npe] + dx; */
/*    acy = part[j+4*npe] + dy; */
/*    acz = part[j+5*npe] + dz; */
      a = _mm512_add_ps(v_dx,_mm512_load_ps(&part[j+3*npe]));
      b = _mm512_add_ps(v_dy,_mm512_load_ps(&part[j+4*npe]));
      c = _mm512_add_ps(v_dz,_mm512_load_ps(&part[j+5*npe]));
/* find inverse gamma */
/*    p2 = acx*acx + acy*acy + acz*acz; */
      v_at = _mm512_fmadd_ps(b,b,_mm512_mul_ps(a,a));
      v_at = _mm512_fmadd_ps(c,c,v_at);
/*    gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
      v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
      v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*    v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* time-centered kinetic energy */
/*    sum1 += gami*p2/(1.0f + gami); */
      v_at = _mm512_mul_ps(v_gami,v_at);
      v_at = _mm512_div_ps(v_at,_mm512_add_ps(v_one,v_gami));
/* convert to double precision before accumulating */
      v_sum1 = _mm512_add_pd(v_sum1,_mm512_cvtpslo_pd(v_at));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_at,78));
      v_sum1 = _mm512_add_pd(v_sum1,v_d);
/* renormalize magnetic field */
/*    qtmg = qtmh*gami; */
      v_at = _mm512_mul_ps(v_qtmh,v_gami);
/* calculate cyclotron frequency */
/*    omxt = qtmg*ox; */
/*    omyt = qtmg*oy; */
/*    omzt = qtmg*oz; */
      e = _mm512_mul_ps(v_at,v_ox);
      f = _mm512_mul_ps(v_at,v_oy);
      g = _mm512_mul_ps(v_at,v_oz);
/* calculate rotation matrix */
/*    vx = omxt*omxt; */
      v_vx = _mm512_mul_ps(e,e);
/*    vy = omyt*omyt; */
      v_vy = _mm512_mul_ps(f,f);
/*    vz = omzt*omzt; */
      v_vz = _mm512_mul_ps(g,g);
/*    omt = omxt*omxt + omyt*omyt + omzt*omzt; */
      v_at = _mm512_add_ps(_mm512_add_ps(v_vx,v_vy),v_vz);
/*    anorm = 2.0f/(1.0f + omt); */
      d = _mm512_div_ps(v_two,_mm512_add_ps(v_one,v_at));
/*    omt = 0.5f*(1.0f - omt); */
      h = _mm512_mul_ps(v_half,_mm512_sub_ps(v_one,v_at));
/*    vx = (omt + vx)*acx; */
      v_vx = _mm512_mul_ps(_mm512_add_ps(h,v_vx),a);
/*    vy = (omt + vy)*acy; */
      v_vy = _mm512_mul_ps(_mm512_add_ps(h,v_vy),b);
/*    vz = (omt + vz)*acz; */
      v_vz = _mm512_mul_ps(_mm512_add_ps(h,v_vz),c);
/*    omt = omxt*omyt; */
      h = _mm512_mul_ps(e,f);
/*    vx = vx + (omzt + omt)*acy; */
      v_vx = _mm512_fmadd_ps(_mm512_add_ps(h,g),b,v_vx);
/*    vy = vy + (omt - omzt)*acx; */
      v_vy = _mm512_fmadd_ps(_mm512_sub_ps(h,g),a,v_vy);
/*    omt = omxt*omzt;  */
      h = _mm512_mul_ps(e,g);
/*    vx = vx + (omt - omyt)*acz; */
      v_vx = _mm512_fmadd_ps(_mm512_sub_ps(h,f),c,v_vx);
/*    vz = vz + (omt + omyt)*acx; */
      v_vz = _mm512_fmadd_ps(_mm512_add_ps(h,f),a,v_vz);
/*    omt = omyt*omzt; */
      h = _mm512_mul_ps(f,g);
/*    vy = vy + (omt + omxt)*acz; */
      v_vy = _mm512_fmadd_ps(_mm512_add_ps(h,e),c,v_vy);
/*    vz = vz + (omt - omxt)*acy; */
      v_vz = _mm512_fmadd_ps(_mm512_sub_ps(h,e),b,v_vz);
/* new momentum  */
/*    vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm; */
/*    vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm; */
/*    vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm; */
      v_vx = _mm512_fmadd_ps(v_vx,d,v_dx);
      v_vy = _mm512_fmadd_ps(v_vy,d,v_dy);
      v_vz = _mm512_fmadd_ps(v_vz,d,v_dz);
/* update inverse gamma */
/*    p2 = vx*vx + vy*vy + vz*vz; */
      v_at = _mm512_fmadd_ps(v_vy,v_vy,_mm512_mul_ps(v_vx,v_vx));
      v_at = _mm512_fmadd_ps(v_vz,v_vz,v_at);
/*    dtg = dtc/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_at = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*    v_at = _mm512_mul_ps(v_dtc,v_at);                            */
/* full accuracy calculation */
      v_at = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
      v_at = _mm512_div_ps(v_dtc,v_at);
/* full accuracy calculation with SVML */
/*    v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/*    v_at = _mm512_div_ps(v_dtc,v_at);                              */
/* new position */
/*    dx = x + vx*dtg; */
/*    dy = y + vy*dtg; */
/*    dz = z + vz*dtg; */
      v_dx = _mm512_fmadd_ps(v_vx,v_at,v_x);
      v_dy = _mm512_fmadd_ps(v_vy,v_at,v_y);
      v_dz = _mm512_fmadd_ps(v_vz,v_at,v_z);
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         v_dx = _mm512_mask_add_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         v_dx = _mm512_mask_sub_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dy < edgely) dy += edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         v_dy = _mm512_mask_add_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dy >= edgery) dy -= edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         v_dy = _mm512_mask_sub_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if ((dz < edgelz) || (dz >= edgerz)) { */
/*          dz = z;                             */
/*          vz = -vz;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
               _MM_CMPINT_GE));
         v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
         v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* set new position */
/*    part[j] = dx;       */
/*    part[j+npe] = dy;   */
/*    part[j+2*npe] = dz; */
      _mm512_store_ps(&part[j],v_dx);
      _mm512_store_ps(&part[j+npe],v_dy);
      _mm512_store_ps(&part[j+2*npe],v_dz);
/* set new momentum */
/*    part[j+3*npe] = vx; */
/*    part[j+4*npe] = vy; */
/*    part[j+5*npe] = vz; */
      _mm512_store_ps(&part[j+3*npe],v_vx);
      _mm512_store_ps(&part[j+4*npe],v_vy);
      _mm512_store_ps(&part[j+5*npe],v_vz);
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nm = 4*(nn + nxv*mm + nxyv*ll);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amz = 1.0f - dzp;
      amy = dxp*amy;
/* find electric field */
      nn = nm;
      dx = amx*fxyz[nn] + amy*fxyz[nn+4];
      dy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      dz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = amz*(dx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = amz*(dy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = amz*(dz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
      nn += 4*nxyv;
      acx = amx*fxyz[nn] + amy*fxyz[nn+4];
      acy = amx*fxyz[nn+1] + amy*fxyz[nn+1+4];
      acz = amx*fxyz[nn+2] + amy*fxyz[nn+2+4];
      mm = nn + 4*nxv;
      dx = dx + dzp*(acx + dyp*fxyz[mm] + dx1*fxyz[mm+4]);
      dy = dy + dzp*(acy + dyp*fxyz[mm+1] + dx1*fxyz[mm+1+4]);
      dz = dz + dzp*(acz + dyp*fxyz[mm+2] + dx1*fxyz[mm+2+4]);
/* find magnetic field */
      nn = nm;
      ox = amx*bxyz[nn] + amy*bxyz[nn+4];
      oy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];
      oz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];
      mm = nn + 4*nxv;
      ox = amz*(ox + dyp*bxyz[mm] + dx1*bxyz[mm+4]);
      oy = amz*(oy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);
      oz = amz*(oz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);
      nn += 4*nxyv;
      acx = amx*bxyz[nn] + amy*bxyz[nn+4];
      acy = amx*bxyz[nn+1] + amy*bxyz[nn+1+4];
      acz = amx*bxyz[nn+2] + amy*bxyz[nn+2+4];
      mm = nn + 4*nxv;
      ox = ox + dzp*(acx + dyp*bxyz[mm] + dx1*bxyz[mm+4]);
      oy = oy + dzp*(acy + dyp*bxyz[mm+1] + dx1*bxyz[mm+1+4]);
      oz = oz + dzp*(acz + dyp*bxyz[mm+2] + dx1*bxyz[mm+2+4]);
/* calculate half impulse */
      dx *= qtmh;
      dy *= qtmh;
      dz *= qtmh;
/* half acceleration */
      acx = part[j+3*npe] + dx;
      acy = part[j+4*npe] + dy;
      acz = part[j+5*npe] + dz;
/* find inverse gamma */
      p2 = acx*acx + acy*acy + acz*acz;
      gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* renormalize magnetic field */
      qtmg = qtmh*gami;
/* time-centered kinetic energy */
      sum1 += gami*p2/(1.0f + gami);
/* calculate cyclotron frequency */
      omxt = qtmg*ox;
      omyt = qtmg*oy;
      omzt = qtmg*oz;
/* calculate rotation matrix */
      omt = omxt*omxt + omyt*omyt + omzt*omzt;
      anorm = 2.0f/(1.0f + omt);
      omt = 0.5f*(1.0f - omt);
      rot4 = omxt*omyt;
      rot7 = omxt*omzt;
      rot8 = omyt*omzt;
      rot1 = omt + omxt*omxt;
      rot5 = omt + omyt*omyt;
      rot9 = omt + omzt*omzt;
      rot2 = omzt + rot4;
      rot4 -= omzt;
      rot3 = -omyt + rot7;
      rot7 += omyt;
      rot6 = omxt + rot8;
      rot8 -= omxt;
/* new momentum  */
      vx = dx + (rot1*acx + rot2*acy + rot3*acz)*anorm;
      vy = dy + (rot4*acx + rot5*acy + rot6*acz)*anorm;
      vz = dz + (rot7*acx + rot8*acy + rot9*acz)*anorm;
/* update inverse gamma */
      p2 = vx*vx + vy*vy + vz*vz;
      dtg = dtc/sqrtf(1.0f + p2*ci2);
/* new position */
      dx = x + vx*dtg;
      dy = y + vy*dtg;
      dz = z + vz*dtg;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            vy = -vy;
         }
         if ((dz < edgelz) || (dz >= edgerz)) {
            dz = z;
            vz = -vz;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            vy = -vy;
         }
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
      part[j+2*npe] = dz;
/* set new momentum */
      part[j+3*npe] = vx;
      part[j+4*npe] = vy;
      part[j+5*npe] = vz;
   }
/* normalize kinetic energy */
/* *ek += sum1; */
   d0 = _mm512_reduce_add_pd(v_sum1);
   *ek += sum1 + d0;
   return;
}

/*--------------------------------------------------------------------*/
void ckncgpost3lt(float part[], float q[], float qm, int nop, int npe,
                  int idimp, int nxv, int nyv, int nzv) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vector version using guard cells
   33 flops/particle, 11 loads, 8 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
   q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
   q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
   q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
   q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
   q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
   q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
   q(n+1,m+1,l+1)=qm*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   q[l][k][j] = charge density at grid point j,k,l
   qm = charge on particle, in units of e
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, nps, nn, mm, ll, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_nxv, v_nxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_as, v_at;
   __m512 a, b, c, d, e, f, g, h, qp, qr;
   __mmask16 msk, msks;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   v_nxv = _mm512_set1_epi32(nxv);
   v_nxyv = _mm512_set1_epi32(nxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = nn + nxv*mm + nxyv*ll; */
      v_it = _mm512_mullo_epi32(v_nxyv,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv,v_mm));
      v_nn = _mm512_add_epi32(v_nn,v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    a = amx*amz; */
/*    b = amy*amz; */
/*    d = dyp*amz; */
/*    d = dx1*amz; */
      a = _mm512_mul_ps(v_amx,v_amz);
      b = _mm512_mul_ps(v_amy,v_amz);
      c = _mm512_mul_ps(v_dyp,v_amz);
      d = _mm512_mul_ps(v_dx1,v_amz);
/*    e = amx*dzp; */
/*    f = amy*dzp; */
/*    g = dyp*dzp; */
/*    h = dx1*dzp; */
      e = _mm512_mul_ps(v_amx,v_dzp);
      f = _mm512_mul_ps(v_amy,v_dzp);
      g = _mm512_mul_ps(v_dyp,v_dzp);
      h = _mm512_mul_ps(v_dx1,v_dzp);
      _mm512_store_epi32(kk,v_nn);
/* deposit charge */
/*    x = q[nn] + amx*amz;       */
/*    y = q[nn+1] + amy*amz;     */
/*    z = q[nn+nxv] + dyp*amz;   */
/*    w = q[nn+1+nxv] + dx1*amz; */
/*    q[nn] = x;                 */
/*    q[nn+1] = y;               */
/*    q[nn+nxv] = z;             */
/*    q[nn+1+nxv] = w;           */
/*    mm = nn + nxyv;            */
/*    x = q[mm] + amx*dzp;       */
/*    y = q[mm+1] + amy*dzp;     */
/*    z = q[mm+nxv] + dyp*dzp;   */
/*    w = q[mm+1+nxv] + dx1*dzp; */
/*    q[mm] = x;                 */
/*    q[mm+1] = y;               */
/*    q[mm+nxv] = z;             */
/*    q[mm+1+nxv] = w;           */
/* deposit charge for two particles at a time */
      for (i = 0; i < 8; i++) {
/* first particle */
         mm = kk[2*i];
         msk = _mm512_int2mask(3<<(2*i));
         msks = _mm512_int2mask(2<<(2*i));
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)a,msks,
                (__m512i)b,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)c,msks,
                (__m512i)d,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
         mm = mm + nxyv;
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)e,msks,
                (__m512i)f,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)g,msks,
                (__m512i)h,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
/* second particle */
         mm = kk[2*i+1];
         msks = _mm512_int2mask(1<<(2*i));
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)b,msks,
                (__m512i)a,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)d,msks,
                (__m512i)c,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
         mm = mm + nxyv;
         qp = _mm512_mask_loadunpacklo_ps(qp,msk,&q[mm]);
         qp = _mm512_mask_loadunpackhi_ps(qp,msk,&q[mm+16]);
         v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)f,msks,
                (__m512i)e,177);
         qp = _mm512_mask_add_ps(qp,msk,qp,v_at);
         _mm512_mask_packstorelo_ps(&q[mm],msk,qp);
         _mm512_mask_packstorehi_ps(&q[mm+16],msk,qp);
         ll = mm + nxv;
         qr = _mm512_mask_loadunpacklo_ps(qr,msk,&q[ll]);
         qr = _mm512_mask_loadunpackhi_ps(qr,msk,&q[ll+16]);
         v_as = (__m512)_mm512_mask_shuffle_epi32((__m512i)h,msks,
                (__m512i)g,177);
         qr = _mm512_mask_add_ps(qr,msk,qr,v_as);
         _mm512_mask_packstorelo_ps(&q[ll],msk,qr);
         _mm512_mask_packstorehi_ps(&q[ll+16],msk,qr);
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = nn + nxv*mm + nxyv*ll;
      amx = qm - dxp;
      amy = 1.0f - dyp;
      amz = 1.0f - dzp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amy = dxp*amy;
/* deposit charge */
      x = q[nn] + amx*amz;
      y = q[nn+1] + amy*amz;
      z = q[nn+nxv] + dyp*amz;
      w = q[nn+1+nxv] + dx1*amz;
      q[nn] = x;
      q[nn+1] = y;
      q[nn+nxv] = z;
      q[nn+1+nxv] = w;
      mm = nn + nxyv;
      x = q[mm] + amx*dzp;
      y = q[mm+1] + amy*dzp;
      z = q[mm+nxv] + dyp*dzp;
      w = q[mm+1+nxv] + dx1*dzp;
      q[mm] = x;
      q[mm+1] = y;
      q[mm+nxv] = z;
      q[mm+1+nxv] = w;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gpost3lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv, int nzv) {
/* for 3d code, this subroutine calculates particle charge density
   using first-order linear interpolation, periodic boundaries
   vector version using guard cells
   33 flops/particle, 11 loads, 8 stores
   input: all, output: q
   charge density is approximated by values at the nearest grid points
   q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
   q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
   q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
   q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
   q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
   q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
   q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
   q(n+1,m+1,l+1)=qm*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   q[l][k][j] = charge density at grid point j,k,l
   qm = charge on particle, in units of e
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nxv = first dimension of charge array, must be >= nx+1
   nyv = second dimension of charge array, must be >= ny+1
   nzv = third dimension of charge array, must be >= nz+1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, nps, nn, mm, ll, nxyv;
   float x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz;
   __m512i v_nxv, v_nxyv;
   __m512i v_nn, v_mm, v_ll, v_it;
   __m512 v_qm, v_one;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1;
   __attribute__((aligned(64))) unsigned int kk[16];
   typedef union vfloat {float v[16]; __m512 v16;} vf;
   vf vv[8];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   v_nxv = _mm512_set1_epi32(nxv);
   v_nxyv = _mm512_set1_epi32(nxyv);
   v_qm = _mm512_set1_ps(qm);
   v_one = _mm512_set1_ps(1.0f);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm; */
/*    dzp = z - (float) ll; */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = nn + nxv*mm + nxyv*ll; */
      v_it = _mm512_mullo_epi32(v_nxyv,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv,v_mm));
      v_nn = _mm512_add_epi32(v_nn,v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    x = amx*amz; */
/*    y = amy*amz; */
/*    z = dyp*amz; */
/*    w = dx1*amz; */
      vv[0].v16 = _mm512_mul_ps(v_amx,v_amz);
      vv[1].v16 = _mm512_mul_ps(v_amy,v_amz);
      vv[2].v16 = _mm512_mul_ps(v_dyp,v_amz);
      vv[3].v16 = _mm512_mul_ps(v_dx1,v_amz);
      vv[4].v16 = _mm512_mul_ps(v_amx,v_dzp);
      vv[5].v16 = _mm512_mul_ps(v_amy,v_dzp);
      vv[6].v16 = _mm512_mul_ps(v_dyp,v_dzp);
      vv[7].v16 = _mm512_mul_ps(v_dx1,v_dzp);
      _mm512_store_epi32(kk,v_nn);
/* deposit charge */
/*    x = q[nn] + amx*amz;       */
/*    y = q[nn+1] + amy*amz;     */
/*    z = q[nn+nxv] + dyp*amz;   */
/*    w = q[nn+1+nxv] + dx1*amz; */
/*    q[nn] = x;                 */
/*    q[nn+1] = y;               */
/*    q[nn+nxv] = z;             */
/*    q[nn+1+nxv] = w;           */
/*    mm = nn + nxyv;            */
/*    x = q[mm] + amx*dzp;       */
/*    y = q[mm+1] + amy*dzp;     */
/*    z = q[mm+nxv] + dyp*dzp;   */
/*    w = q[mm+1+nxv] + dx1*dzp; */
/*    q[mm] = x;                 */
/*    q[mm+1] = y;               */
/*    q[mm+nxv] = z;             */
/*    q[mm+1+nxv] = w;           */
      for (i = 0; i < 16; i++) {
         nn = kk[i];
         x = q[nn] + vv[0].v[i];
         y = q[nn+1] + vv[1].v[i];
         z = q[nn+nxv] + vv[2].v[i];
         w = q[nn+1+nxv] + vv[3].v[i];
         q[nn] = x;
         q[nn+1] = y;
         q[nn+nxv] = z;
         q[nn+1+nxv] = w;
         mm = nn + nxyv;
         x = q[mm] + vv[4].v[i];
         y = q[mm+1] + vv[5].v[i];
         z = q[mm+nxv] + vv[6].v[i];
         w = q[mm+1+nxv] + vv[7].v[i];
         q[mm] = x;
         q[mm+1] = y;
         q[mm+nxv] = z;
         q[mm+1+nxv] = w;
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = nn + nxv*mm + nxyv*ll;
      amx = qm - dxp;
      amy = 1.0f - dyp;
      amz = 1.0f - dzp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amy = dxp*amy;
/* deposit charge */
      x = q[nn] + amx*amz;
      y = q[nn+1] + amy*amz;
      z = q[nn+nxv] + dyp*amz;
      w = q[nn+1+nxv] + dx1*amz;
      q[nn] = x;
      q[nn+1] = y;
      q[nn+nxv] = z;
      q[nn+1+nxv] = w;
      mm = nn + nxyv;
      x = q[mm] + amx*dzp;
      y = q[mm+1] + amy*dzp;
      z = q[mm+nxv] + dyp*dzp;
      w = q[mm+1+nxv] + dx1*dzp;
      q[mm] = x;
      q[mm+1] = y;
      q[mm+nxv] = z;
      q[mm+1+nxv] = w;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncgjpost3lt(float part[], float cu[], float qm, float dt,
                   int nop, int npe, int idimp, int nx, int ny, int nz,
                   int nxv, int nyv, int nzv, int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation
   in addition, particle positions are advanced a half time-step
   vector version using guard cells
   69 flops/particle, 30 loads, 27 stores
   input: all, output: part, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   and qci = qm*vi, where i = x,y,z
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   part[3][n] = x velocity of particle n
   part[4][n] = y velocity of particle n
   part[5][n] = z velocity of particle n
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nx/ny/nz = system length in x/y/z direction
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   nzv = fourth dimension of current array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
   int i, j, nps, nn, mm, ll, ii, nxyv;
   float edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float x, y, z;
   __m512i v_nxv4, v_nxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_perm;
   __m512 v_qm, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   v_nxv4 = _mm512_set1_epi32(4*nxv);
   v_nxyv4 = _mm512_set1_epi32(4*nxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm;      */
/*    dzp = z - (float) ll;      */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/*    nn = 4*(nn + nxv*mm + nxyv*ll); */
      v_it = _mm512_mullo_epi32(v_nxyv4,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv4,v_mm));
      v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    a = amx*amz; */
/*    b = amy*amz; */
/*    c = dyp*amz; */
/*    d = dx1*amz; */
      a = _mm512_mul_ps(v_amx,v_amz);
      b = _mm512_mul_ps(v_amy,v_amz);
      c = _mm512_mul_ps(v_dyp,v_amz);
      d = _mm512_mul_ps(v_dx1,v_amz);
/*    e = amx*dzp; */
/*    f = amy*dzp; */
/*    g = dyp*dzp; */
/*    h = dx1*dzp; */
      e = _mm512_mul_ps(v_amx,v_dzp);
      f = _mm512_mul_ps(v_amy,v_dzp);
      g = _mm512_mul_ps(v_dyp,v_dzp);
      h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*    vx = part[j+3*npe]; */
/*    vy = part[j+4*npe]; */
/*    vz = part[j+5*npe]; */
      v_vx = _mm512_load_ps(&part[j+3*npe]);
      v_vy = _mm512_load_ps(&part[j+4*npe]);
      v_vz = _mm512_load_ps(&part[j+5*npe]);
      v_ll = _mm512_add_epi32(v_nn,v_nxyv4);
/* deposit charge for one particle at a time */
      for (i = 0; i < 16; i++) {
         ii = i >> 2;
         if (i==(ii<<2)) {
            switch (ii)
            {
               case 0:
/* replicate velocities of first group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,0);
                  q = _mm512_permute4f128_ps(v_vy,0);
                  r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(a,
                      _mm512_int2mask(61680),b,177);
                  t = _mm512_mask_permute4f128_ps(c,
                      _mm512_int2mask(61680),d,177);
                  u = _mm512_mask_permute4f128_ps(e,
                      _mm512_int2mask(61680),f,177);
                  v = _mm512_mask_permute4f128_ps(g,
                      _mm512_int2mask(61680),h,177);
                  break;
               case 1:
/* replicate velocities of second group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,85);
                  q = _mm512_permute4f128_ps(v_vy,85);
                  r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(b,
                      _mm512_int2mask(3855),a,177);
                  t = _mm512_mask_permute4f128_ps(d,
                      _mm512_int2mask(3855),c,177);
                  u = _mm512_mask_permute4f128_ps(f,
                      _mm512_int2mask(3855),e,177);
                  v = _mm512_mask_permute4f128_ps(h,
                      _mm512_int2mask(3855),g,177);
                  break;
               case 2:
/* replicate velocities of third group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,170);
                  q = _mm512_permute4f128_ps(v_vy,170);
                  r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(a,
                      _mm512_int2mask(61680),b,177);
                  s = _mm512_permute4f128_ps(s,78);
                  t = _mm512_mask_permute4f128_ps(c,
                      _mm512_int2mask(61680),d,177);
                  t = _mm512_permute4f128_ps(t,78);
                  u = _mm512_mask_permute4f128_ps(e,
                      _mm512_int2mask(61680),f,177);
                  u = _mm512_permute4f128_ps(u,78);
                  v = _mm512_mask_permute4f128_ps(g,
                      _mm512_int2mask(61680),h,177);
                  v = _mm512_permute4f128_ps(v,78);
                  break;
               case 3:
/* replicate velocities of fourth group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,255);
                  q = _mm512_permute4f128_ps(v_vy,255);
                  r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(b,
                      _mm512_int2mask(3855),a,177);
                  s = _mm512_permute4f128_ps(s,78);
                  t = _mm512_mask_permute4f128_ps(d,
                      _mm512_int2mask(3855),c,177);
                  t = _mm512_permute4f128_ps(t,78);
                  u = _mm512_mask_permute4f128_ps(f,
                      _mm512_int2mask(3855),e,177);
                  u = _mm512_permute4f128_ps(u,78);
                  v = _mm512_mask_permute4f128_ps(h,
                      _mm512_int2mask(3855),g,177);
                  v = _mm512_permute4f128_ps(v,78);
                  break;
            }
         }
         v_it = _mm512_setzero_epi32();
         switch (i-(ii<<2))
         {
/* first particle */
            case 0:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                      _mm512_int2mask(170),(__m512i)q,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,0);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,0);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,0);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,0);
               break;
/* second particle */
            case 1:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                      _mm512_int2mask(85),(__m512i)p,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,85);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,85);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,85);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,85);
               break;
/* third particle */
            case 2:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                      _mm512_int2mask(170),(__m512i)q,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                      _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,170);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,170);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,170);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,170);
               break;
/* fourth particle */
            case 3:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                      _mm512_int2mask(85),(__m512i)p,177);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,255);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,255);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,255);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,255);
               break;
         }
         _mm512_store_epi32(kk,v_nn);
/* load cu[nn:nn+3] and cu[nn+4:nn+7] field components */
/*    dx = amx*amz;        */
/*    dy = amy*amz;        */
/*    cu[nn] += vx*dx;     */
/*    cu[nn+1] += vy*dx;   */
/*    cu[nn+2] += vz*dx;   */
/*    dx = dyp*amz;        */
/*    cu[nn+4] += vx*dy;   */
/*    cu[nn+1+4] += vy*dy; */
/*    cu[nn+2+4] += vz*dy; */
         mm = kk[i];
         cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
              &cu[mm]);
         cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
              &cu[mm+16]);
         cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cp);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cp);
/*    mm = nn + 4*nxv;                                 */
/* load cu[mm:mm+3] and cu[mm+4:mm+7] field components */
/*    dx = dyp*amz;        */
/*    dy = dx1*amz;        */
/*    cu[mm] += vx*dx;     */
/*    cu[mm+1] += vy*dx;   */
/*    cu[mm+2] += vz*dx;   */
/*    cu[mm+4] += vx*dy;   */
/*    cu[mm+1+4] += vy*dy; */
/*    cu[mm+2+4] += vz*dy; */
         mm = kk[i] + 4*nxv;
         cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
              &cu[mm]);
         cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
              &cu[mm+16]);
         cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cr);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cr);
         _mm512_store_epi32(kk,v_ll);
/*    nn += 4*nxyv;                                    */
/* load cu[nn:nn+3] and cu[nn+4:nn+7] field components */
/*    nn += 4*nxyv; */
/*    dx = amx*dzp;        */
/*    dy = amy*dzp;        */
/*    cu[nn] += vx*dx;     */
/*    cu[nn+1] += vy*dx;   */
/*    cu[nn+2] += vz*dx;   */
/*    cu[nn+4] += vx*dy;   */
/*    cu[nn+1+4] += vy*dy; */
/*    cu[nn+2+4] += vz*dy; */
         mm = kk[i];
         cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
              &cu[mm]);
         cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
              &cu[mm+16]);
         cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cp);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cp);
/*    mm = nn + 4*nxv;                                 */
/* load cu[mm:mm+3] and cu[mm+4:mm+7] field components */
/*    dx = dyp*dzp;        */
/*    dy = dx1*dzp;        */
/*    cu[mm] += vx*dx;     */
/*    cu[mm+1] += vy*dx;   */
/*    cu[mm+2] += vz*dx;   */
/*    cu[mm+4] += vx*dy;   */
/*    cu[mm+1+4] += vy*dy; */
/*    cu[mm+2+4] += vz*dy; */
         mm = kk[i] + 4*nxv;
         cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
              &cu[mm]);
         cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
              &cu[mm+16]);
         cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cr);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cr);
      }
/* advance position half a time-step */
/*    dx = x + vx*dt; */
/*    dy = y + vy*dt; */
/*    dz = z + vz*dt; */
      v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
      v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
      v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         v_dx = _mm512_mask_add_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         v_dx = _mm512_mask_sub_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dy < edgely) dy += edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         v_dy = _mm512_mask_add_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dy >= edgery) dy -= edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         v_dy = _mm512_mask_sub_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+3*npe] = -vx;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+3*npe],v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+4*npe] = -vy;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+4*npe],v_vy);
/*       if ((dz < edgelz) || (dz >= edgerz)) { */
/*          dz = z;                             */
/*          part[j+5*npe] = -vz;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
               _MM_CMPINT_GE));
         v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
         v_vz = _mm512_mask_sub_ps(v_vz,msk,v_zero,v_vz);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+5*npe],v_vz);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+3*npe] = -vx;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_vx = _mm512_mask_sub_ps(v_vx,msk,v_zero,v_vx);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+3*npe],v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+4*npe] = -vy;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_vy = _mm512_mask_sub_ps(v_vy,msk,v_zero,v_vy);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+4*npe],v_vy);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* set new position */
/*    part[j] = dx;       */
/*    part[j+npe] = dy;   */
/*    part[j+2*npe] = dz; */
      _mm512_store_ps(&part[j],v_dx);
      _mm512_store_ps(&part[j+npe],v_dy);
      _mm512_store_ps(&part[j+2*npe],v_dz);
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
      nn = 4*(nn + nxv*mm + nxyv*ll);
      amx = qm - dxp;
      amy = 1.0f - dyp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amz = 1.0f - dzp;
      amy = dxp*amy;
/* deposit current */
      dx = amx*amz;
      dy = amy*amz;
      vx = part[j+3*npe];
      vy = part[j+4*npe];
      vz = part[j+5*npe];
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = dyp*amz;
      cu[nn+4] += vx*dy;
      cu[nn+1+4] += vy*dy;
      cu[nn+2+4] += vz*dy;
      dy = dx1*amz;
      mm = nn + 4*nxv;
      cu[mm] += vx*dx;
      cu[mm+1] += vy*dx;
      cu[mm+2] += vz*dx;
      dx = amx*dzp;
      cu[mm+4] += vx*dy;
      cu[mm+1+4] += vy*dy;
      cu[mm+2+4] += vz*dy;
      dy = amy*dzp;
      nn += 4*nxyv;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = dyp*dzp;
      cu[nn+4] += vx*dy;
      cu[nn+1+4] += vy*dy;
      cu[nn+2+4] += vz*dy;
      dy = dx1*dzp;
      mm = nn + 4*nxv;
      cu[mm] += vx*dx;
      cu[mm+1] += vy*dx;
      cu[mm+2] += vz*dx;
      cu[mm+4] += vx*dy;
      cu[mm+1+4] += vy*dy;
      cu[mm+2+4] += vz*dy;
/* advance position half a time-step */
      dx = x + vx*dt;
      dy = y + vy*dt;
      dz = z + vz*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+3*npe] = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+4*npe] = -vy;
         }
         if ((dz < edgelz) || (dz >= edgerz)) {
            dz = z;
            part[j+5*npe] = -vz;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+3*npe] = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+4*npe] = -vy;
         }
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
      part[j+2*npe] = dz;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrjpost3lt(float part[], float cu[], float qm, float dt,
                    float ci, int nop, int npe, int idimp, int nx,
                    int ny, int nz, int nxv, int nyv, int nzv,
                    int ipbc) {
/* for 3d code, this subroutine calculates particle current density
   using first-order linear interpolation for relativistic particles
   in addition, particle positions are advanced a half time-step
   vector version using guard cells
   79 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
   input: all, output: part, cu
   current density is approximated by values at the nearest grid points
   cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
   cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
   cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
   cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
   cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
   cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
   cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
   cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
   where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
   and qci = qm*pi*gami, where i = x,y,z
   where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   part[3][n] = x momentum of particle n
   part[4][n] = y momentum of particle n
   part[5][n] = z momentum of particle n
   cu[l][k][j][i] = ith component of current density at grid point j,k,l
   qm = charge on particle, in units of e
   dt = time interval between successive calculations
   ci = reciprocal of velocity of light
   nop = number of particles
   npe = first dimension of particle array
   idimp = size of phase space = 6
   nx/ny/nz = system length in x/y/z direction
   nxv = second dimension of current array, must be >= nx+1
   nyv = third dimension of current array, must be >= ny+1
   nzv = fourth dimension of current array, must be >= nz+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
   cu needs to have 4 components, although one is not used
local data                                                            */
   int i, j, nps, nn, mm, ll, ii, nxyv;
   float ci2, edgelx, edgely, edgelz, edgerx, edgery, edgerz;
   float dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz, vx, vy, vz;
   float p2, gami;
   float x, y, z, ux, uy, uz;
   __m512i v_nxv4, v_nxyv4;
   __m512i v_nn, v_mm, v_ll, v_it, v_perm;
   __m512 v_qm, v_ci2, v_dt, v_one, v_zero;
   __m512 v_x, v_y, v_z, v_dxp, v_dyp, v_dzp, v_amx, v_amy, v_amz;
   __m512 v_dx1, v_gami, v_at, v_dx, v_dy, v_dz, v_vx, v_vy, v_vz;
   __m512 v_ux, v_uy, v_uz;
   __m512 v_edgelx, v_edgely, v_edgelz, v_edgerx, v_edgery, v_edgerz;
   __m512 a, b, c, d, e, f, g, h, p, q, r, s, t, u, v, ws, wt, wu, wv;
   __m512 cp, cr;
   __mmask16 msk;
   __attribute__((aligned(64))) unsigned int kk[16];
   nxyv = nxv*nyv;
   nps = 16*(nop/16);
   ci2 = ci*ci;
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgelz = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   edgerz = (float) nz;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgelz = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
      edgerz = (float) (nz-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   v_nxv4 = _mm512_set1_epi32(4*nxv);
   v_nxyv4 = _mm512_set1_epi32(4*nxyv);
   v_perm = _mm512_set_epi32(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);
   v_qm = _mm512_set1_ps(qm);
   v_ci2 = _mm512_set1_ps(ci2);
   v_dt = _mm512_set1_ps(dt);
   v_one = _mm512_set1_ps(1.0f);
   v_zero = _mm512_setzero_ps();
   v_edgelx = _mm512_set1_ps(edgelx);
   v_edgely = _mm512_set1_ps(edgely);
   v_edgelz = _mm512_set1_ps(edgelz);
   v_edgerx = _mm512_set1_ps(edgerx);
   v_edgery = _mm512_set1_ps(edgery);
   v_edgerz = _mm512_set1_ps(edgerz);
/* vector loop over particles in blocks of 16 */
   for (j = 0; j < nps; j+=16) {
/* find interpolation weights */
/*    x = part[j];       */
/*    y = part[j+npe];   */
/*    z = part[j+2*npe]; */
      v_x = _mm512_load_ps(&part[j]);
      v_y = _mm512_load_ps(&part[j+npe]);
      v_z = _mm512_load_ps(&part[j+2*npe]);
/*    nn = x; */
/*    mm = y; */
/*    ll = z; */
      v_nn = _mm512_cvtfxpnt_round_adjustps_epi32(v_x,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_mm = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_ll = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
             _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    dxp = qm*(x - (float) nn); */
/*    dyp = y - (float) mm;      */
/*    dzp = z - (float) ll;      */
      v_dxp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_nn,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dxp = _mm512_mul_ps(v_qm,_mm512_sub_ps(v_x,v_dxp));
      v_dyp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_mm,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dyp = _mm512_sub_ps(v_y,v_dyp);
      v_dzp = _mm512_cvtfxpnt_round_adjustepi32_ps(v_ll,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dzp = _mm512_sub_ps(v_z,v_dzp);
/* find inverse gamma */
/*    ux = part[j+3*npe]; */
/*    uy = part[j+4*npe]; */
/*    uz = part[j+5*npe]; */
      v_ux = _mm512_load_ps(&part[j+3*npe]);
      v_uy = _mm512_load_ps(&part[j+4*npe]);
      v_uz = _mm512_load_ps(&part[j+5*npe]);
/*    p2 = ux*ux + uy*uy + uz*uz;       */
      v_at = _mm512_fmadd_ps(v_uy,v_uy,_mm512_mul_ps(v_ux,v_ux));
      v_at = _mm512_fmadd_ps(v_uz,v_uz,v_at);
/*    gami = 1.0f/sqrtf(1.0f + p2*ci2); */
/* approximate calculation */
/*    v_gami = _mm512_rsqrt23_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* full accuracy calculation */
      v_gami = _mm512_sqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one));
      v_gami = _mm512_div_ps(v_one,v_gami);
/* full accuracy calculation with SVML */
/*    v_gami = _mm512_invsqrt_ps(_mm512_fmadd_ps(v_at,v_ci2,v_one)); */
/* calculate weights */
/*    nn = 4*(nn + nxv*mm + nxyv*ll); */
      v_it = _mm512_mullo_epi32(v_nxyv4,v_ll);
      v_it = _mm512_add_epi32(v_it,_mm512_mullo_epi32(v_nxv4,v_mm));
      v_nn = _mm512_add_epi32(_mm512_slli_epi32(v_nn,2),v_it);
/*    amx = qm - dxp;   */
/*    amy = 1.0f - dyp; */
/*    amz = 1.0f - dzp; */
      v_amx = _mm512_sub_ps(v_qm,v_dxp);
      v_amy = _mm512_sub_ps(v_one,v_dyp);
      v_amz = _mm512_sub_ps(v_one,v_dzp);
/*    dx1 = dxp*dyp; */
/*    dyp = amx*dyp; */
/*    amx = amx*amy; */
/*    amy = dxp*amy; */
      v_dx1 = _mm512_mul_ps(v_dxp,v_dyp);
      v_dyp = _mm512_mul_ps(v_amx,v_dyp);
      v_amx = _mm512_mul_ps(v_amx,v_amy);
      v_amy = _mm512_mul_ps(v_dxp,v_amy);
/*    a = amx*amz; */
/*    b = amy*amz; */
/*    c = dyp*amz; */
/*    d = dx1*amz; */
      a = _mm512_mul_ps(v_amx,v_amz);
      b = _mm512_mul_ps(v_amy,v_amz);
      c = _mm512_mul_ps(v_dyp,v_amz);
      d = _mm512_mul_ps(v_dx1,v_amz);
/*    e = amx*dzp; */
/*    f = amy*dzp; */
/*    g = dyp*dzp; */
/*    h = dx1*dzp; */
      e = _mm512_mul_ps(v_amx,v_dzp);
      f = _mm512_mul_ps(v_amy,v_dzp);
      g = _mm512_mul_ps(v_dyp,v_dzp);
      h = _mm512_mul_ps(v_dx1,v_dzp);
/* deposit current */
/*    vx = ux*gami; */
/*    vy = uy*gami; */
/*    vz = uz*gami; */
      v_vx = _mm512_mul_ps(v_ux,v_gami);
      v_vy = _mm512_mul_ps(v_uy,v_gami);
      v_vz = _mm512_mul_ps(v_uz,v_gami);
      v_ll = _mm512_add_epi32(v_nn,v_nxyv4);
/* deposit charge for one particle at a time */
      for (i = 0; i < 16; i++) {
         ii = i >> 2;
         if (i==(ii<<2)) {
            switch (ii)
            {
               case 0:
/* replicate velocities of first group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,0);
                  q = _mm512_permute4f128_ps(v_vy,0);
                  r = _mm512_permute4f128_ps(v_vz,0);
/* regroup weights for first group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(a,
                      _mm512_int2mask(61680),b,177);
                  t = _mm512_mask_permute4f128_ps(c,
                      _mm512_int2mask(61680),d,177);
                  u = _mm512_mask_permute4f128_ps(e,
                      _mm512_int2mask(61680),f,177);
                  v = _mm512_mask_permute4f128_ps(g,
                      _mm512_int2mask(61680),h,177);
                  break;
               case 1:
/* replicate velocities of second group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,85);
                  q = _mm512_permute4f128_ps(v_vy,85);
                  r = _mm512_permute4f128_ps(v_vz,85);
/* regroup weights for second group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(b,
                      _mm512_int2mask(3855),a,177);
                  t = _mm512_mask_permute4f128_ps(d,
                      _mm512_int2mask(3855),c,177);
                  u = _mm512_mask_permute4f128_ps(f,
                      _mm512_int2mask(3855),e,177);
                  v = _mm512_mask_permute4f128_ps(h,
                      _mm512_int2mask(3855),g,177);
                  break;
               case 2:
/* replicate velocities of third group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,170);
                  q = _mm512_permute4f128_ps(v_vy,170);
                  r = _mm512_permute4f128_ps(v_vz,170);
/* regroup weights for third group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(a,
                      _mm512_int2mask(61680),b,177);
                  s = _mm512_permute4f128_ps(s,78);
                  t = _mm512_mask_permute4f128_ps(c,
                      _mm512_int2mask(61680),d,177);
                  t = _mm512_permute4f128_ps(t,78);
                  u = _mm512_mask_permute4f128_ps(e,
                      _mm512_int2mask(61680),f,177);
                  u = _mm512_permute4f128_ps(u,78);
                  v = _mm512_mask_permute4f128_ps(g,
                      _mm512_int2mask(61680),h,177);
                  v = _mm512_permute4f128_ps(v,78);
                  break;
               case 3:
/* replicate velocities of fourth group of 4 particles */
                  p = _mm512_permute4f128_ps(v_vx,255);
                  q = _mm512_permute4f128_ps(v_vy,255);
                  r = _mm512_permute4f128_ps(v_vz,255);
/* regroup weights for fourth group of 4 particles */
                  s = _mm512_mask_permute4f128_ps(b,
                      _mm512_int2mask(3855),a,177);
                  s = _mm512_permute4f128_ps(s,78);
                  t = _mm512_mask_permute4f128_ps(d,
                      _mm512_int2mask(3855),c,177);
                  t = _mm512_permute4f128_ps(t,78);
                  u = _mm512_mask_permute4f128_ps(f,
                      _mm512_int2mask(3855),e,177);
                  u = _mm512_permute4f128_ps(u,78);
                  v = _mm512_mask_permute4f128_ps(h,
                      _mm512_int2mask(3855),g,177);
                  v = _mm512_permute4f128_ps(v,78);
                  break;
            }
         }
         v_it = _mm512_setzero_epi32();
         switch (i-(ii<<2))
         {
/* first particle */
            case 0:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                      _mm512_int2mask(170),(__m512i)q,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,78);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,0);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,0);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,0);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,0);
               break;
/* second particle */
            case 1:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                      _mm512_int2mask(85),(__m512i)p,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,24);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,85);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,85);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,85);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,85);
               break;
/* third particle */
            case 2:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)p,
                      _mm512_int2mask(170),(__m512i)q,177);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)r,
                      _mm512_int2mask(51),(__m512i)v_at,78);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,170);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,170);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,170);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,170);
               break;
/* fourth particle */
            case 3:
/* reorder velocity components */
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)q,
                      _mm512_int2mask(85),(__m512i)p,177);
               v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_at,78);
               v_at = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at,
                      _mm512_int2mask(68),(__m512i)r,177);
/* reorder weights */
               ws = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)s,255);
               wt = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)t,255);
               wu = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)u,255);
               wv = (__m512)_mm512_mask_shuffle_epi32(v_it,
                    _mm512_int2mask(119),(__m512i)v,255);
               break;
         }
         _mm512_store_epi32(kk,v_nn);
/* load cu[nn:nn+3] and cu[nn+4:nn+7] field components */
/*    dx = amx*amz;        */
/*    dy = amy*amz;        */
/*    cu[nn] += vx*dx;     */
/*    cu[nn+1] += vy*dx;   */
/*    cu[nn+2] += vz*dx;   */
/*    dx = dyp*amz;        */
/*    cu[nn+4] += vx*dy;   */
/*    cu[nn+1+4] += vy*dy; */
/*    cu[nn+2+4] += vz*dy; */
         mm = kk[i];
         cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
              &cu[mm]);
         cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
              &cu[mm+16]);
         cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),ws,cp);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cp);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cp);
/*    mm = nn + 4*nxv;                                 */
/* load cu[mm:mm+3] and cu[mm+4:mm+7] field components */
/*    dx = dyp*amz;        */
/*    dy = dx1*amz;        */
/*    cu[mm] += vx*dx;     */
/*    cu[mm+1] += vy*dx;   */
/*    cu[mm+2] += vz*dx;   */
/*    cu[mm+4] += vx*dy;   */
/*    cu[mm+1+4] += vy*dy; */
/*    cu[mm+2+4] += vz*dy; */
         mm = kk[i] + 4*nxv;
         cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
              &cu[mm]);
         cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
              &cu[mm+16]);
         cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wt,cr);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cr);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cr);
         _mm512_store_epi32(kk,v_ll);
/*    nn += 4*nxyv;                                    */
/* load cu[nn:nn+3] and cu[nn+4:nn+7] field components */
/*    nn += 4*nxyv; */
/*    dx = amx*dzp;        */
/*    dy = amy*dzp;        */
/*    cu[nn] += vx*dx;     */
/*    cu[nn+1] += vy*dx;   */
/*    cu[nn+2] += vz*dx;   */
/*    cu[nn+4] += vx*dy;   */
/*    cu[nn+1+4] += vy*dy; */
/*    cu[nn+2+4] += vz*dy; */
         mm = kk[i];
         cp = _mm512_mask_loadunpacklo_ps(cp,_mm512_int2mask(255),
              &cu[mm]);
         cp = _mm512_mask_loadunpackhi_ps(cp,_mm512_int2mask(255),
              &cu[mm+16]);
         cp = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wu,cp);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cp);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cp);
/*    mm = nn + 4*nxv;                                 */
/* load cu[mm:mm+3] and cu[mm+4:mm+7] field components */
/*    dx = dyp*dzp;        */
/*    dy = dx1*dzp;        */
/*    cu[mm] += vx*dx;     */
/*    cu[mm+1] += vy*dx;   */
/*    cu[mm+2] += vz*dx;   */
/*    cu[mm+4] += vx*dy;   */
/*    cu[mm+1+4] += vy*dy; */
/*    cu[mm+2+4] += vz*dy; */
         mm = kk[i] + 4*nxv;
         cr = _mm512_mask_loadunpacklo_ps(cr,_mm512_int2mask(255),
              &cu[mm]);
         cr = _mm512_mask_loadunpackhi_ps(cr,_mm512_int2mask(255),
              &cu[mm+16]);
         cr = _mm512_mask_fmadd_ps(v_at,_mm512_int2mask(119),wv,cr);
         _mm512_mask_packstorelo_ps(&cu[mm],_mm512_int2mask(255),cr);
         _mm512_mask_packstorehi_ps(&cu[mm+16],_mm512_int2mask(255),cr);
      }
/* advance position half a time-step */
/*    dx = x + vx*dt; */
/*    dy = y + vy*dt; */
/*    dz = z + vz*dt; */
      v_dx = _mm512_fmadd_ps(v_vx,v_dt,v_x);
      v_dy = _mm512_fmadd_ps(v_vy,v_dt,v_y);
      v_dz = _mm512_fmadd_ps(v_vz,v_dt,v_z);
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         v_dx = _mm512_mask_add_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgerx,_MM_CMPINT_GE);
         v_dx = _mm512_mask_sub_ps(v_dx,msk,v_dx,v_edgerx);
/*       if (dy < edgely) dy += edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         v_dy = _mm512_mask_add_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dy >= edgery) dy -= edgery; */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgery,_MM_CMPINT_GE);
         v_dy = _mm512_mask_sub_ps(v_dy,msk,v_dy,v_edgery);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+3*npe] = -ux;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+3*npe],v_ux);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+4*npe] = -uy;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+4*npe],v_uy);
/*       if ((dz < edgelz) || (dz >= edgerz)) { */
/*          dz = z;                             */
/*          part[j+5*npe] = -uz;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dz,v_edgerz,
               _MM_CMPINT_GE));
         v_dz = _mm512_mask_blend_ps(msk,v_dz,v_z);
         v_uz = _mm512_mask_sub_ps(v_uz,msk,v_zero,v_uz);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+5*npe],v_uz);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          part[j+3*npe] = -ux;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dx,v_edgelx,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dx,v_edgerx,
               _MM_CMPINT_GE));
         v_dx = _mm512_mask_blend_ps(msk,v_dx,v_x);
         v_ux = _mm512_mask_sub_ps(v_ux,msk,v_zero,v_ux);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+3*npe],v_ux);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          part[j+4*npe] = -uy;                */
/*       }                                      */
         msk = _mm512_cmp_ps_mask(v_dy,v_edgely,_MM_CMPINT_LT);
         msk = _mm512_kor(msk,_mm512_cmp_ps_mask(v_dy,v_edgery,
               _MM_CMPINT_GE));
         v_dy = _mm512_mask_blend_ps(msk,v_dy,v_y);
         v_uy = _mm512_mask_sub_ps(v_uy,msk,v_zero,v_uy);
/* write output if test result is true for any particle */
         if (msk)
            _mm512_store_ps(&part[j+4*npe],v_uy);
/*       if (dz < edgelz) dz += edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgelz,_MM_CMPINT_LT);
         v_dz = _mm512_mask_add_ps(v_dz,msk,v_dz,v_edgerz);
/*       if (dz >= edgerz) dz -= edgerz; */
         msk = _mm512_cmp_ps_mask(v_dz,v_edgerz,_MM_CMPINT_GE);
         v_dz = _mm512_mask_sub_ps(v_dz,msk,v_dz,v_edgerz);
      }
/* set new position */
/*    part[j] = dx;       */
/*    part[j+npe] = dy;   */
/*    part[j+2*npe] = dz; */
      _mm512_store_ps(&part[j],v_dx);
      _mm512_store_ps(&part[j+npe],v_dy);
      _mm512_store_ps(&part[j+2*npe],v_dz);
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      z = part[j+2*npe];
      nn = x;
      mm = y;
      ll = z;
      dxp = qm*(x - (float) nn);
      dyp = y - (float) mm;
      dzp = z - (float) ll;
/* find inverse gamma */
      ux = part[j+3*npe];
      uy = part[j+4*npe];
      uz = part[j+5*npe];
      p2 = ux*ux + uy*uy + uz*uz;
      gami = 1.0f/sqrtf(1.0f + p2*ci2);
/* calculate weights */
      nn = 4*(nn + nxv*mm + nxyv*ll);
      amx = qm - dxp;
      amy = 1.0f - dyp;
      dx1 = dxp*dyp;
      dyp = amx*dyp;
      amx = amx*amy;
      amz = 1.0f - dzp;
      amy = dxp*amy;
/* deposit current */
      dx = amx*amz;
      dy = amy*amz;
      vx = ux*gami;
      vy = uy*gami;
      vz = uz*gami;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = dyp*amz;
      cu[nn+4] += vx*dy;
      cu[nn+1+4] += vy*dy;
      cu[nn+2+4] += vz*dy;
      dy = dx1*amz;
      mm = nn + 4*nxv;
      cu[mm] += vx*dx;
      cu[mm+1] += vy*dx;
      cu[mm+2] += vz*dx;
      dx = amx*dzp;
      cu[mm+4] += vx*dy;
      cu[mm+1+4] += vy*dy;
      cu[mm+2+4] += vz*dy;
      dy = amy*dzp;
      nn += 4*nxyv;
      cu[nn] += vx*dx;
      cu[nn+1] += vy*dx;
      cu[nn+2] += vz*dx;
      dx = dyp*dzp;
      cu[nn+4] += vx*dy;
      cu[nn+1+4] += vy*dy;
      cu[nn+2+4] += vz*dy;
      dy = dx1*dzp;
      mm = nn + 4*nxv;
      cu[mm] += vx*dx;
      cu[mm+1] += vy*dx;
      cu[mm+2] += vz*dx;
      cu[mm+4] += vx*dy;
      cu[mm+1+4] += vy*dy;
      cu[mm+2+4] += vz*dy;
/* advance position half a time-step */
      dx = x + vx*dt;
      dy = y + vy*dt;
      dz = z + vz*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+3*npe] = -ux;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+4*npe] = -uy;
         }
         if ((dz < edgelz) || (dz >= edgerz)) {
            dz = z;
            part[j+5*npe] = -uz;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            part[j+3*npe] = -ux;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            part[j+4*npe] = -uy;
         }
         if (dz < edgelz) dz += edgerz;
         if (dz >= edgerz) dz -= edgerz;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
      part[j+2*npe] = dz;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncdsortp3yzlt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1, int nyz1) {
/* this subroutine sorts particles by y,z grid
   linear interpolation
   part = particle array
   part[1][n] = position y of particle n
   part[2][n] = position z of particle n
   npic = address offset for reordering particles
   idimp = size of phase space = 6
   npe = first dimension of particle array
   ny1 = system length in y direction + 1
   nyz1 = ny1*nz1, where nz1 = system length in z direction + 1
   requires KNC, part needs to be 64 byte aligned
   npe needs to be a multiple of 16
local data                                                            */
   int i, j, k, m, l, nps, ip;
/* int isum, ist; */
   __m512i v_ny1, v_m, v_l, v_it;
   __m512 v_y, v_z, v_at;
   __attribute__((aligned(64))) unsigned int pp[16];
   typedef union vint {int ll[16]; __m512i v_l;} vi;
   vi ii;
   nps = 16*(nop/16);
   v_ny1 = _mm512_set1_epi32(ny1);
/* clear counter array */
/* for (k = 0; k < nyz1; k++) { */
/*    npic[k] = 0;              */
/* }                            */
   memset((void *)npic,0,nyz1*sizeof(int));
/* find how many particles in each grid */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=16) {
/*    m = parta[j+npe];   */
/*    l = parta[j+2*npe]; */
      v_y = _mm512_load_ps(&parta[j+npe]);
      v_z = _mm512_load_ps(&parta[j+2*npe]);
      v_m = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_l = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    l = m + ny1*l; */
      ii.v_l =_mm512_fmadd_epi32(v_ny1,v_l,v_m);
      for (k = 0; k < 16; k++) {
         l = ii.ll[k];
         npic[l] += 1;
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
      m = parta[j+npe];
      l = parta[j+2*npe];
      l = m + ny1*l;
      npic[l] += 1;
   }
/* find address offset */
/* isum = 0;                    */
/* for (k = 0; k < nyz1; k++) { */
/*    ist = npic[k];            */
/*    npic[k] = isum;           */
/*    isum += ist;              */
/* }                            */
   ckncxiscan2(npic,nyz1);
/* find addresses of particles at each grid and reorder particles */
/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=16) {
/*    m = parta[j+npe];   */
/*    l = parta[j+2*npe]; */
      v_y = _mm512_load_ps(&parta[j+npe]);
      v_z = _mm512_load_ps(&parta[j+2*npe]);
      v_m = _mm512_cvtfxpnt_round_adjustps_epi32(v_y,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
      v_l = _mm512_cvtfxpnt_round_adjustps_epi32(v_z,
            _MM_ROUND_MODE_DOWN,_MM_EXPADJ_NONE);
/*    l = m + ny1*l; */
      ii.v_l = _mm512_fmadd_epi32(v_ny1,v_l,v_m);
      for (k = 0; k < 16; k++) {
         l = ii.ll[k];
         ip = npic[l];
         npic[l] = ip + 1;
         pp[k] = ip;
      }
      v_it = _mm512_load_epi32(pp);
      for (i = 0; i < idimp; i++) {
/*          partb[ip+npe*i] = parta[j+npe*i]; */
            v_at = _mm512_load_ps(&parta[j+npe*i]);
            _mm512_i32scatter_ps(&partb[npe*i],v_it,v_at,4);
      }
   }
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
      m = parta[j+npe];
      l = parta[j+2*npe];
      l = m + ny1*l;
      ip = npic[l];
      for (i = 0; i < idimp; i++) {
         partb[ip+npe*i] = parta[j+npe*i];
      }
      npic[l] = ip + 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cknccguard3l(float fxyz[], int nx, int ny, int nz, int nxe,
                  int nye, int nze) {
/* replicate extended periodic vector field fxyz
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, fxyz needs to be 64 byte aligned
   nxe needs to be a multiple of 4
local data                                                 */
   int j, k, l, nxs, nxyen, ll;
   nxs = 4*(nx/4);
   nxyen = 4*nxe*nye;
/* copy edges of extended field */
   for (l = 0; l < nz; l++) {
      ll = nxyen*l;
      for (k = 0; k < ny; k++) {
         fxyz[4*nx+4*nxe*k+ll] = fxyz[4*nxe*k+ll];
         fxyz[1+4*nx+4*nxe*k+ll] = fxyz[1+4*nxe*k+ll];
         fxyz[2+4*nx+4*nxe*k+ll] = fxyz[2+4*nxe*k+ll];
      }
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[4*j+4*nxe*ny+ll],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j+ll]));
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         fxyz[4*j+4*nxe*ny+ll] = fxyz[4*j+ll];
         fxyz[1+4*j+4*nxe*ny+ll] = fxyz[1+4*j+ll];
         fxyz[2+4*j+4*nxe*ny+ll] = fxyz[2+4*j+ll];
      }
      fxyz[4*nx+4*nxe*ny+ll] = fxyz[ll];
      fxyz[1+4*nx+4*nxe*ny+ll] = fxyz[1+ll];
      fxyz[2+4*nx+4*nxe*ny+ll] = fxyz[2+ll];
   }
   for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[4*j+4*nxe*k+nxyen*nz],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j+4*nxe*k]));
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         fxyz[4*j+4*nxe*k+nxyen*nz] = fxyz[4*j+4*nxe*k];
         fxyz[1+4*j+4*nxe*k+nxyen*nz] = fxyz[1+4*j+4*nxe*k];
         fxyz[2+4*j+4*nxe*k+nxyen*nz] = fxyz[2+4*j+4*nxe*k];
      }
      fxyz[4*nx+4*nxe*k+nxyen*nz] = fxyz[4*nxe*k];
      fxyz[1+4*nx+4*nxe*k+nxyen*nz] = fxyz[1+4*nxe*k];
      fxyz[2+4*nx+4*nxe*k+nxyen*nz] = fxyz[2+4*nxe*k];
   }
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         _mm512_mask_store_ps(&fxyz[4*j+4*nxe*ny+nxyen*nz],
            _mm512_int2mask(30583),_mm512_load_ps(&fxyz[4*j]));
   }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
      fxyz[4*j+4*nxe*ny+nxyen*nz] = fxyz[4*j];
      fxyz[1+4*j+4*nxe*ny+nxyen*nz] = fxyz[1+4*j];
      fxyz[2+4*j+4*nxe*ny+nxyen*nz] = fxyz[2+4*j];
   }
   fxyz[4*nx+4*nxe*ny+nxyen*nz] = fxyz[0];
   fxyz[1+4*nx+4*nxe*ny+nxyen*nz] = fxyz[1];
   fxyz[2+4*nx+4*nxe*ny+nxyen*nz] = fxyz[2];
   return;
}

/*--------------------------------------------------------------------*/
void ckncacguard3l(float cu[], int nx, int ny, int nz, int nxe, int nye,
                   int nze) {
/* accumulate extended periodic field cu
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, fxyz needs to be 64 byte aligned
   nxe needs to be a multiple of 4
local data                                                 */
   int j, k, l, nxs, nxyen, ll;
   __m512 v_cu, v_zero;
   nxs = 4*(nx/4);
   nxyen = 4*nxe*nye;
   v_zero = _mm512_set1_ps(0.0f);
/* accumulate edges of extended field */
   for (l = 0; l < nz; l++) {
      ll = nxyen*l;
      for (k = 0; k < ny; k++) {
         cu[4*nxe*k+ll] += cu[4*nx+4*nxe*k+ll];
         cu[1+4*nxe*k+ll] += cu[1+4*nx+4*nxe*k+ll];
         cu[2+4*nxe*k+ll] += cu[2+4*nx+4*nxe*k+ll];
         cu[4*nx+4*nxe*k+ll] = 0.0;
         cu[1+4*nx+4*nxe*k+ll] = 0.0;
         cu[2+4*nx+4*nxe*k+ll] = 0.0;
      }
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         v_cu = _mm512_load_ps(&cu[4*j+4*nxe*ny+ll]);
         v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j+ll]),v_cu);
         _mm512_store_ps(&cu[4*j+ll],v_cu);
         _mm512_store_ps(&cu[4*j+4*nxe*ny+ll],v_zero);
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         cu[4*j+ll] += cu[4*j+4*nxe*ny+ll];
         cu[1+4*j+ll] += cu[1+4*j+4*nxe*ny+ll];
         cu[2+4*j+ll] += cu[2+4*j+4*nxe*ny+ll];
         cu[4*j+4*nxe*ny+ll] = 0.0;
         cu[1+4*j+4*nxe*ny+ll] = 0.0;
         cu[2+4*j+4*nxe*ny+ll] = 0.0;
      }
      cu[ll] += cu[4*nx+4*nxe*ny+ll];
      cu[1+ll] += cu[1+4*nx+4*nxe*ny+ll];
      cu[2+ll] += cu[2+4*nx+4*nxe*ny+ll];
      cu[4*nx+4*nxe*ny+ll] = 0.0;
      cu[1+4*nx+4*nxe*ny+ll] = 0.0;
      cu[2+4*nx+4*nxe*ny+ll] = 0.0;
   }
   for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 4 */
      for (j = 0; j < nxs; j+=4) {
         v_cu = _mm512_load_ps(&cu[4*j+4*nxe*k+nxyen*nz]);
         v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j+4*nxe*k]),v_cu);
         _mm512_store_ps(&cu[4*j+4*nxe*k],v_cu);
         _mm512_store_ps(&cu[4*j+4*nxe*k+nxyen*nz],v_zero);
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         cu[4*j+4*nxe*k] += cu[4*j+4*nxe*k+nxyen*nz];
         cu[1+4*j+4*nxe*k] += cu[1+4*j+4*nxe*k+nxyen*nz];
         cu[2+4*j+4*nxe*k] += cu[2+4*j+4*nxe*k+nxyen*nz];
         cu[4*j+4*nxe*k+nxyen*nz] = 0.0;
         cu[1+4*j+4*nxe*k+nxyen*nz] = 0.0;
         cu[2+4*j+4*nxe*k+nxyen*nz] = 0.0;
      }
      cu[4*nxe*k] += cu[4*nx+4*nxe*k+nxyen*nz];
      cu[1+4*nxe*k] += cu[1+4*nx+4*nxe*k+nxyen*nz];
      cu[2+4*nxe*k] += cu[2+4*nx+4*nxe*k+nxyen*nz];
      cu[4*nx+4*nxe*k+nxyen*nz] = 0.0;
      cu[1+4*nx+4*nxe*k+nxyen*nz] = 0.0;
      cu[2+4*nx+4*nxe*k+nxyen*nz] = 0.0;
   }
/* vector loop over elements in blocks of 4 */
   for (j = 0; j < nxs; j+=4) {
      v_cu = _mm512_load_ps(&cu[4*j+4*nxe*ny+nxyen*nz]);
      v_cu = _mm512_add_ps(_mm512_load_ps(&cu[4*j]),v_cu);
      _mm512_store_ps(&cu[4*j],v_cu);
      _mm512_store_ps(&cu[4*j+4*nxe*ny+nxyen*nz],v_zero);
   }
/* loop over remaining elements */
   for (j = nxs; j < nx; j++) {
      cu[4*j] += cu[4*j+4*nxe*ny+nxyen*nz];
      cu[1+4*j] += cu[1+4*j+4*nxe*ny+nxyen*nz];
      cu[2+4*j] += cu[2+4*j+4*nxe*ny+nxyen*nz];
      cu[4*j+4*nxe*ny+nxyen*nz] = 0.0;
      cu[1+4*j+4*nxe*ny+nxyen*nz] = 0.0;
      cu[2+4*j+4*nxe*ny+nxyen*nz] = 0.0;
   }
   cu[0] += cu[4*nx+4*nxe*ny+nxyen*nz];
   cu[1] += cu[1+4*nx+4*nxe*ny+nxyen*nz];
   cu[2] += cu[2+4*nx+4*nxe*ny+nxyen*nz];
   cu[4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   cu[1+4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   cu[2+4*nx+4*nxe*ny+nxyen*nz] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void ckncaguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
                  int nze) {
/* accumulate extended periodic scalar field q
   linear interpolation
   nx/ny/nz = system length in x/y direction
   nxe = first dimension of field arrays, must be >= nx+1
   nye = second dimension of field arrays, must be >= ny+1
   nze = third dimension of field arrays, must be >= nz+1
   requires KNC, q needs to be 64 byte aligned
   nxe needs to be a multiple of 16
local data                                                 */
   int j, k, l, nxs, nxye, ll;
   __m512 v_q;
   nxs = 16*(nx/16);
   nxye = nxe*nye;
/* accumulate edges of extended field */
   for (l = 0; l < nz; l++) {
      ll = nxye*l;
      for (k = 0; k < ny; k++) {
         q[nxe*k+ll] += q[nx+nxe*k+ll];
         q[nx+nxe*k+ll] = 0.0;
      }
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < nxs; j+=16) {
         v_q = _mm512_load_ps(&q[j+nxe*ny+ll]);
         v_q = _mm512_add_ps(_mm512_load_ps(&q[j+ll]),v_q);
         _mm512_store_ps(&q[j+ll],v_q);
         _mm512_store_ps(&q[j+nxe*ny+ll],_mm512_setzero_ps());
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         q[j+ll] += q[j+nxe*ny+ll];
         q[j+nxe*ny+ll] = 0.0;
      }
      q[ll] += q[nx+nxe*ny+ll];
      q[nx+nxe*ny+ll] = 0.0;
   }
   for (k = 0; k < ny; k++) {
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < nxs; j+=16) {
         v_q = _mm512_load_ps(&q[j+nxe*k+nxye*nz]);
         v_q = _mm512_add_ps(_mm512_load_ps(&q[j+nxe*k]),v_q);
         _mm512_store_ps(&q[j+nxe*k],v_q);
         _mm512_store_ps(&q[j+nxe*k+nxye*nz],_mm512_setzero_ps());
      }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
         q[j+nxe*k] += q[j+nxe*k+nxye*nz];
         q[j+nxe*k+nxye*nz] = 0.0;
      }
      q[nxe*k] += q[nx+nxe*k+nxye*nz];
      q[nx+nxe*k+nxye*nz] = 0.0;
   }
/* vector loop over elements in blocks of 16 */
      for (j = 0; j < nxs; j+=16) {
         v_q = _mm512_load_ps(&q[j+nxe*ny+nxye*nz]);
         v_q = _mm512_add_ps(_mm512_load_ps(&q[j]),v_q);
         _mm512_store_ps(&q[j],v_q);
         _mm512_store_ps(&q[j+nxe*ny+nxye*nz],_mm512_setzero_ps());
   }
/* loop over remaining elements */
      for (j = nxs; j < nx; j++) {
      q[j] += q[j+nxe*ny+nxye*nz];
      q[j+nxe*ny+nxye*nz] = 0.0;
   }
   q[0] += q[nx+nxe*ny+nxye*nz];
   q[nx+nxe*ny+nxye*nz] = 0.0;
   return;
}

/*--------------------------------------------------------------------*/
void ckncpois33(float complex q[], float complex fxyz[], int isign,
                float complex ffc[], float ax, float ay, float az,
                float affp, float *we, int nx, int ny, int nz, int nxvh,
                int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   force/charge (or convolution of electric field over particle shape)
   with periodic boundary conditions.
   for isign = 0, output: ffc
   input: isign,ax,ay,az,affp,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   for isign = -1, output: fxyz, we
   input: q,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   approximate flop count is:
   59*nxc*nyc*nzc + 26*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   if isign = 0, form factor array is prepared
   if isign is not equal to 0, force/charge is calculated
   equation used is:
   fx[kz][ky][kx] = -sqrt(-1)*kx*g[kz][ky][kx]*s[kz][ky][kx],
   fy[kz][ky][kx] = -sqrt(-1)*ky*g[kz][ky][kx]*s[kz][ky][kx],
   fz[kz][ky][kx] = -sqrt(-1)*kz*g[kz][ky][kx]*s[kz][ky][kx],
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s[kz][ky][kx],
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
   fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
   fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
   fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
   q[l][k][j] = complex charge density for fourier mode (j,k,l)
   fxyz[l][k][j][0] = x component of complex force/charge
   fxyz[l][k][j][1] = y component of complex force/charge
   fxyz[l][k][j][2] = z component of complex force/charge
   all for fourier mode (j,k,l)
   cimag(ffc[l][k][j]) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   creal(ffc[l][k][j]) = potential green's function g
   for fourier mode (j,k,l)
   ax/ay/az = half-width of particle in x/y/z direction
   affp = normalization constant = nx*ny*nz/np,
   where np=number of particles
   electric field energy is also calculated, using
   we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
      |q[kz][ky][kx]*s[kz][ky][kx]|**2)
   nx/ny/nz = system length in x/y/z direction
   nxvh = first dimension of field arrays, must be >= nxh
   nyv = second dimension of field arrays, must be >= ny
   nzv = third dimension of field arrays, must be >= nz
   nxhd = first dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
   requires KNC, q, fxy, ffc need to be 64 byte aligned
   nxhd, nxvh need to be a multiple of 8
   fxyz needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh; 
   float dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6;
   float complex zero, zt1, zt2;
   double wp, d0;
   __m512i v_j, v_it, v_perm;
   __m512 v_dnx, v_dny, v_dnz, v_dky, v_dkz, v_at1, v_at2, v_at3, v_at4;
   __m512 v_zero, v_zt1, v_zt2, v_zt3, v_zt4;
   __m512 a, b, c, d, e, f, g, h;
   __m512d v_wp, v_d;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 8*(nxh/8);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   v_j = _mm512_set_epi32(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_perm = _mm512_set_epi32(15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
   if (isign != 0)
      goto L40;
/* prepare form factor array */
   for (l = 0; l < nzh; l++) {
      dkz = dnz*(float) l;
      ll = nxyhd*l;
      at1 = dkz*dkz;
      at2 = pow((dkz*az),2);
      for (k = 0; k < nyh; k++) {
         dky = dny*(float) k;
         kk = nxhd*k;
         at3 = dky*dky + at1;
         at4 = pow((dky*ay),2) + at2;
         for (j = 0; j < nxh; j++) {
            dkx = dnx*(float) j;
            at5 = dkx*dkx + at3;
            at6 = exp(-0.5*(pow((dkx*ax),2) + at4));
            if (at5==0.0) {
               ffc[j+kk+ll] = affp + 1.0*_Complex_I;
            }
            else {
               ffc[j+kk+ll] = (affp*at6/at5) + at6*_Complex_I;
            }
         }
      }
   }
   return;
/* calculate force/charge and sum field energy */
L40: wp = 0.0;
   v_wp = _mm512_setzero_pd();
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
   for (l = 1; l < nzh; l++) {
      dkz = dnz*(float) l;
      v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
      ll = nxyhd*l;
      lj = nxvyh*l;
      l1 = nxvyh*nz - lj;
      for (k = 1; k < nyh; k++) {
         dky = dny*(float) k;
         v_it = _mm512_set1_epi32(k);
         v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dky = _mm512_mul_ps(v_dny,v_dky);
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhs; j+=8) {
/*          at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]); */
            v_at1 = _mm512_load_ps((float *)&ffc[j+kk+ll]);
            v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
            v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*          at2 = at1*dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*          at3 = dky*at1; */
            v_at3 = _mm512_mul_ps(v_dky,v_at1);
/*          at4 = dkz*at1; */
            v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*          zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I; */
            v_zt1 = _mm512_load_ps((float *)&q[j+kj+lj]);
            v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                    v_zero,v_zt1);
            v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*          zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I; */
            v_zt2 = _mm512_load_ps((float *)&q[j+k1+lj]);
            v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                    v_zero,v_zt2);
            v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
               v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
            }
/*          fxyz[4*(j+kj+lj)] = at2*zt1;   */
/*          fxyz[1+4*(j+kj+lj)] = at3*zt1; */
/*          fxyz[2+4*(j+kj+lj)] = at4*zt1; */
            a = _mm512_mul_ps(v_at2,v_zt1);
            b = _mm512_mul_ps(v_at3,v_zt1);
            c = _mm512_mul_ps(v_at4,v_zt1);
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                 c,78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                a,78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                 v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                 b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+kj+lj)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+kj+lj)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+kj+lj)],d);
/*          fxyz[4*(j+k1+lj)] = at2*zt2;    */
/*          fxyz[1+4*(j+k1+lj)] = -at3*zt2; */
/*          fxyz[2+4*(j+k1+lj)] = at4*zt2;  */
            a = _mm512_mul_ps(v_at2,v_zt2);
            b = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at3,v_zt2));
            c = _mm512_mul_ps(v_at4,v_zt2);
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                 c,78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                a,78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                 v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                 b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],d);
/*          wp += at1*(q[j+kj+lj]*conjf(q[j+kj+lj]) */
/*             + q[j+k1+lj]*conjf(q[j+k1+lj]));     */
            v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
            v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
            v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/*          zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I; */
            v_zt1 = _mm512_load_ps((float *)&q[j+kj+l1]);
            v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                    v_zero,v_zt1);
            v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*          zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I; */
            v_zt2 = _mm512_load_ps((float *)&q[j+k1+l1]);
            v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                    v_zero,v_zt2);
            v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
               v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
            }
/*          fxyz[4*(j+kj+l1)] = at2*zt1;    */
/*          fxyz[1+4*(j+kj+l1)] = at3*zt1;  */
/*          fxyz[2+4*(j+kj+l1)] = -at4*zt1; */
            a = _mm512_mul_ps(v_at2,v_zt1);
            b = _mm512_mul_ps(v_at3,v_zt1);
            c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt1));
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                 c,78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                a,78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                 v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                 b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+kj+l1)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+kj+l1)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+kj+l1)],d);
/*          fxyz[4*(j+k1+l1)] = at2*zt2;    */
/*          fxyz[1+4*(j+k1+l1)] = -at3*zt2; */
/*          fxyz[2+4*(j+k1+l1)] = -at4*zt2; */
            a = _mm512_mul_ps(v_at2,v_zt2);
            b = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at3,v_zt2));
            c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt2));
/* perform 4x16 transpose for fxyz field components */
            e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),
                 c,78);
            f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),
                a,78);
            g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
                 v_zero,78);
            h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),
                 b,78);
            a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
                177);
            b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,
                177);
            c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
                177);
            d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,
                177);
            a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
            b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
            c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
            d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],a);
            _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],b);
            _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],c);
            _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],d);
/*          wp += at1*(q[j+kj+l1]*conjf(q[j+kj+l1]) */
/*             + q[j+k1+l1]*conjf(q[j+k1+l1]));     */
            v_zt4 = _mm512_mul_ps(v_zt1,v_zt1);
            v_zt4 = _mm512_add_ps(v_zt4,_mm512_mul_ps(v_zt2,v_zt2));
            v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_at1,v_zt4));
/* convert to double precision before accumulating */
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = crealf(ffc[j+kk+ll])*cimagf(ffc[j+kk+ll]);
            at2 = at1*dnx*(float) j;
            at3 = dky*at1;
            at4 = dkz*at1;
            zt1 = cimagf(q[j+kj+lj]) - crealf(q[j+kj+lj])*_Complex_I;
            zt2 = cimagf(q[j+k1+lj]) - crealf(q[j+k1+lj])*_Complex_I;
            fxyz[4*(j+kj+lj)] = at2*zt1;
            fxyz[1+4*(j+kj+lj)] = at3*zt1;
            fxyz[2+4*(j+kj+lj)] = at4*zt1;
            fxyz[4*(j+k1+lj)] = at2*zt2;
            fxyz[1+4*(j+k1+lj)] = -at3*zt2;
            fxyz[2+4*(j+k1+lj)] = at4*zt2;
            zt1 = cimagf(q[j+kj+l1]) - crealf(q[j+kj+l1])*_Complex_I;
            zt2 = cimagf(q[j+k1+l1]) - crealf(q[j+k1+l1])*_Complex_I;
            fxyz[4*(j+kj+l1)] = at2*zt1;
            fxyz[1+4*(j+kj+l1)] = at3*zt1;
            fxyz[2+4*(j+kj+l1)] = -at4*zt1;
            fxyz[4*(j+k1+l1)] = at2*zt2;
            fxyz[1+4*(j+k1+l1)] = -at3*zt2;
            fxyz[2+4*(j+k1+l1)] = -at4*zt2;
            at1 = at1*(q[j+kj+lj]*conjf(q[j+kj+lj])                 
                + q[j+k1+lj]*conjf(q[j+k1+lj])
                + q[j+kj+l1]*conjf(q[j+kj+l1])      
                + q[j+k1+l1]*conjf(q[j+k1+l1]));
            wp += (double) at1;
         }
      }
/* mode numbers kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         at1 = crealf(ffc[kk+ll])*cimagf(ffc[kk+ll]);
         at3 = at1*dny*(float) k;
         at4 = dkz*at1;
         zt1 = cimagf(q[kj+lj]) - crealf(q[kj+lj])*_Complex_I;
         zt2 = cimagf(q[kj+l1]) - crealf(q[kj+l1])*_Complex_I;
         fxyz[4*(kj+lj)] = zero;
         fxyz[1+4*(kj+lj)] = at3*zt1;
         fxyz[2+4*(kj+lj)] = at4*zt1;
         fxyz[4*(k1+lj)] = zero;
         fxyz[1+4*(k1+lj)] = zero;
         fxyz[2+4*(k1+lj)] = zero;
         fxyz[4*(kj+l1)] = zero;
         fxyz[1+4*(kj+l1)] = at3*zt2;
         fxyz[2+4*(kj+l1)] = -at4*zt2;
         fxyz[4*(k1+l1)] = zero;
         fxyz[1+4*(k1+l1)] = zero;
         fxyz[2+4*(k1+l1)] = zero;
         at1 = at1*(q[kj+lj]*conjf(q[kj+lj])
             + q[kj+l1]*conjf(q[kj+l1]));
         wp += (double) at1;
      }
/* mode numbers ky = 0, ny/2 */
      k1 = nxvh*nyh;
/* vector loop over elements in blocks of 8 */
      for (j = 0; j < nxhs; j+=8) {
/*       at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]); */
         v_at1 = _mm512_load_ps((float *)&ffc[j+ll]);
         v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*       at2 = at1*dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*       at4 = dkz*at1; */
         v_at4 = _mm512_mul_ps(v_dkz,v_at1);
/*       zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I; */
         v_zt1 = _mm512_load_ps((float *)&q[j+lj]);
         v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),
                 v_zero,v_zt1);
         v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/*       zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I; */
         v_zt2 = _mm512_load_ps((float *)&q[j+l1]);
         v_zt2 = _mm512_mask_sub_ps(v_zt2,_mm512_int2mask(21845),
                 v_zero,v_zt2);
         v_zt2 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt2,177);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
            v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(3),v_zero);
         }
/*       fxyz[4*(j+lj)] = at2*zt1;   */
/*       fxyz[1+4*(j+lj)] = zero;    */
/*       fxyz[2+4*(j+lj)] = at4*zt1; */
         a = _mm512_mul_ps(v_at2,v_zt1);
         b = v_zero;
         c = _mm512_mul_ps(v_at4,v_zt1);
/* perform 4x16 transpose for fxyz field components */
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
         f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
         g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
             v_zero,78);
         h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),b,
             78);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
             177);
         b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,177);
         c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
             177);
         d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,177);
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         _mm512_store_ps((float *)&fxyz[4*(j+lj)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+lj)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+lj)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+lj)],d);
/*       fxyz[4*(j+k1+lj)] = zero;   */
/*       fxyz[1+4*(j+k1+lj)] = zero; */
/*       fxyz[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1+lj)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1+lj)],v_zero);
/*       fxyz[4*(j+l1)] = at2*zt2;    */
/*       fxyz[1+4*(j+l1)] = zero;     */
/*       fxyz[2+4*(j+l1)] = -at4*zt2; */
         a = _mm512_mul_ps(v_at2,v_zt2);
         b = v_zero;
         c = _mm512_sub_ps(v_zero,_mm512_mul_ps(v_at4,v_zt2));
/* perform 4x16 transpose for fxyz field components */
         e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
         f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
         g = _mm512_mask_permute4f128_ps(b,_mm512_int2mask(65280),
             v_zero,78);
         h = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(255),b,
             78);
         a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),g,
             177);
         b = _mm512_mask_permute4f128_ps(g,_mm512_int2mask(3855),e,177);
         c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),h,
             177);
         d = _mm512_mask_permute4f128_ps(h,_mm512_int2mask(3855),f,177);
         a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
         b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
         c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
         d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],a);
         _mm512_store_ps((float *)&fxyz[8+4*(j+l1)],b);
         _mm512_store_ps((float *)&fxyz[16+4*(j+l1)],c);
         _mm512_store_ps((float *)&fxyz[24+4*(j+l1)],d);
/*       fxyz[4*(j+k1+l1)] = zero;   */
/*       fxyz[1+4*(j+k1+l1)] = zero; */
/*       fxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
         _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*       wp += at1*(q[j+lj]*conjf(q[j+lj]) */
/*          + q[j+l1]*conjf(q[j+l1]));     */
         v_zt3 = _mm512_mul_ps(v_zt1,v_zt1);
         v_zt3 = _mm512_add_ps(v_zt3,_mm512_mul_ps(v_zt2,v_zt2));
         v_zt3 = _mm512_mul_ps(v_at1,v_zt3);
/* convert to double precision before accumulating */
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = crealf(ffc[j+ll])*cimagf(ffc[j+ll]);
         at2 = at1*dnx*(float) j;  
         at4 = dkz*at1;
         zt1 = cimagf(q[j+lj]) - crealf(q[j+lj])*_Complex_I;
         zt2 = cimagf(q[j+l1]) - crealf(q[j+l1])*_Complex_I;
         fxyz[4*(j+lj)] = at2*zt1;
         fxyz[1+4*(j+lj)] = zero;
         fxyz[2+4*(j+lj)] = at4*zt1;
         fxyz[4*(j+k1+lj)] = zero;
         fxyz[1+4*(j+k1+lj)] = zero;
         fxyz[2+4*(j+k1+lj)] = zero;
         fxyz[4*(j+l1)] = at2*zt2;
         fxyz[1+4*(j+l1)] = zero;
         fxyz[2+4*(j+l1)] = -at4*zt2;
         fxyz[4*(j+k1+l1)] = zero;
         fxyz[1+4*(j+k1+l1)] = zero;
         fxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(q[j+lj]*conjf(q[j+lj])                           
             + q[j+l1]*conjf(q[j+l1]));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      at1 = crealf(ffc[ll])*cimagf(ffc[ll]);
      at4 = dkz*at1;
      zt1 = cimagf(q[lj]) - crealf(q[lj])*_Complex_I;
      fxyz[4*lj] = zero;
      fxyz[1+4*lj] = zero;
      fxyz[2+4*lj] = at4*zt1;
      fxyz[4*(k1+lj)] = zero;
      fxyz[1+4*(k1+lj)] = zero;
      fxyz[2+4*(k1+lj)] = zero;
      fxyz[4*l1] = zero;
      fxyz[1+4*l1] = zero;
      fxyz[2+4*l1] = zero;
      fxyz[4*(k1+l1)] = zero;
      fxyz[1+4*(k1+l1)] = zero;
      fxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(q[lj]*conjf(q[lj]));
      wp += (double) at1;
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      for (j = 1; j < nxh; j++) {
         at1 = crealf(ffc[j+kk])*cimagf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         zt1 = cimagf(q[j+kj]) - crealf(q[j+kj])*_Complex_I;
         zt2 = cimagf(q[j+k1]) - crealf(q[j+k1])*_Complex_I;
         fxyz[4*(j+kj)] = at2*zt1;
         fxyz[1+4*(j+kj)] = at3*zt1;
         fxyz[2+4*(j+kj)] = zero;
         fxyz[4*(j+k1)] = at2*zt2;
         fxyz[1+4*(j+k1)] = -at3*zt2;
         fxyz[2+4*(j+k1)] = zero;
         fxyz[4*(j+kj+l1)] = zero;
         fxyz[1+4*(j+kj+l1)] = zero;
         fxyz[2+4*(j+kj+l1)] = zero;
         fxyz[4*(j+k1+l1)] = zero;
         fxyz[1+4*(j+k1+l1)] = zero;
         fxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(q[j+kj]*conjf(q[j+kj]) + q[j+k1]*conjf(q[j+k1]));
         wp += (double) at1;
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = crealf(ffc[kk])*cimagf(ffc[kk]);
      at3 = at1*dny*(float) k;
      zt1 = cimagf(q[kj]) - crealf(q[kj])*_Complex_I;
      fxyz[4*kj] = zero;
      fxyz[1+4*kj] = at3*zt1;
      fxyz[2+4*kj] = zero;
      fxyz[4*k1] = zero;
      fxyz[1+4*k1] = zero;
      fxyz[2+4*k1] = zero;
      fxyz[4*(kj+l1)] = zero;
      fxyz[1+4*(kj+l1)] = zero;
      fxyz[2+4*(kj+l1)] = zero;
      fxyz[4*(k1+l1)] = zero;
      fxyz[1+4*(k1+l1)] = zero;
      fxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(q[kj]*conjf(q[kj]));
      wp += (double) at1;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 8 */
   for (j = 0; j < nxhs; j+=8) {
/*    at1 = crealf(ffc[j])*cimagf(ffc[j]); */
      v_at1 = _mm512_load_ps((float *)&ffc[j]);
      v_at2 = (__m512)_mm512_shuffle_epi32((__m512i)v_at1,177);
      v_at1 = _mm512_mul_ps(v_at1,v_at2);
/*    at2 = at1*dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_at2 = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_at2 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_dnx,v_at2));
/*    zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I; */
      v_zt1 = _mm512_load_ps((float *)&q[j]);
      v_zt1 = _mm512_mask_sub_ps(v_zt1,_mm512_int2mask(21845),v_zero,
              v_zt1);
      v_zt1 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,177);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(3),v_zero);
      }
/*    fxyz[4*j] = at2*zt1; */
/*    fxyz[1+4*j] = zero;  */
/*    fxyz[2+4*j] = zero;  */
      a = _mm512_mul_ps(v_at2,v_zt1);
      b = v_zero;
      c = v_zero;
/* perform 4x16 transpose for fxyz field components */
      e = _mm512_mask_permute4f128_ps(a,_mm512_int2mask(65280),c,78);
      f = _mm512_mask_permute4f128_ps(c,_mm512_int2mask(255),a,78);
      a = _mm512_mask_permute4f128_ps(e,_mm512_int2mask(61680),v_zero,
          177);
      b = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(3855),e,
          177);
      c = _mm512_mask_permute4f128_ps(f,_mm512_int2mask(61680),v_zero,
          177);
      d = _mm512_mask_permute4f128_ps(v_zero,_mm512_int2mask(3855),f,
          177);
      a = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)a);
      b = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)b);
      c = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)c);
      d = (__m512)_mm512_permutevar_epi32(v_perm,(__m512i)d);
      _mm512_store_ps((float *)&fxyz[4*j],a);
      _mm512_store_ps((float *)&fxyz[8+4*j],b);
      _mm512_store_ps((float *)&fxyz[16+4*j],c);
      _mm512_store_ps((float *)&fxyz[24+4*j],d);
/*    fxyz[4*(j+k1)] = zero;   */
/*    fxyz[1+4*(j+k1)] = zero; */
/*    fxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+k1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+k1)],v_zero);
/*    fxyz[4*(j+l1)] = zero;    */
/*    fxyz[1+4*(j+l1)] = zero;  */
/*     fxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+l1)],v_zero);
/*    fxyz[4*(j+k1+l1)] = zero;   */
/*    fxyz[1+4*(j+k1+l1)] = zero; */
/*    fxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[8+4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[16+4*(j+k1+l1)],v_zero);
      _mm512_store_ps((float *)&fxyz[24+4*(j+k1+l1)],v_zero);
/*    wp += at1*(q[j]*conjf(q[j])); */
      v_zt3 = _mm512_mul_ps(v_at1,_mm512_mul_ps(v_zt1,v_zt1));
/* convert to double precision before accumulating */
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt3));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt3,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
   }
/* loop over remaining elements */
   for (j = itn; j < nxh; j++) {
      at1 = crealf(ffc[j])*cimagf(ffc[j]);
      at2 = at1*dnx*(float) j;
      zt1 = cimagf(q[j]) - crealf(q[j])*_Complex_I;
      fxyz[4*j] = at2*zt1;
      fxyz[1+4*j] = zero;
      fxyz[2+4*j] = zero;
      fxyz[4*(j+k1)] = zero;
      fxyz[1+4*(j+k1)] = zero;
      fxyz[2+4*(j+k1)] = zero;
      fxyz[4*(j+l1)] = zero;
      fxyz[1+4*(j+l1)] = zero;
      fxyz[2+4*(j+l1)] = zero;
      fxyz[4*(j+k1+l1)] = zero;
      fxyz[1+4*(j+k1+l1)] = zero;
      fxyz[2+4*(j+k1+l1)] = zero;
      at1 = at1*(q[j]*conjf(q[j]));
      wp += (double) at1;
   }
   fxyz[0] = zero;
   fxyz[1] = zero;
   fxyz[2] = zero;
   fxyz[4*k1] = zero;
   fxyz[1+4*k1] = zero;
   fxyz[2+4*k1] = zero;
   fxyz[4*l1] = zero;
   fxyz[1+4*l1] = zero;
   fxyz[2+4*l1] = zero;
   fxyz[4*(k1+l1)] = zero;
   fxyz[1+4*(k1+l1)] = zero;
   fxyz[2+4*(k1+l1)] = zero;
/* *we = wp*((float) nx)*((float) ny)*((float) nz); */
   d0 = _mm512_reduce_add_pd(v_wp);
   *we = (wp + d0)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void cknccuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
                 int nyv, int nzv) {
/* this subroutine calculates the transverse current in fourier space
   input: all, output: cu
   approximate flop count is:
   100*nxc*nyc*nzc + 36*(nxc*nyc + nxc*nzc + nyc*nzc)
   and (nx/2)*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the transverse current is calculated using the equation:
   cux[kz][ky][kx] = cux[kz][ky][kx]
                   - kx*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
                   +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   cuy([kz][ky][kx] = cuy[kz][ky][kx]
                    - ky*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
                    +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   cuz[kz][ky][kx] = cuz[kz][ky][kx]
                   - kz*(kx*cux[kz][ky][kx]+ky*cuy[kz][ky][kx]
                   +     kz*cuz[kz][ky][kx])/(kx*kx+ky*ky+kz*kz)
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers, except for
   cux(kx=pi) = cuy(kx=pi) = cuz(kx=pi) = 0,
   cux(ky=pi) = cuy(ky=pi) = cux(ky=pi) = 0,
   cux(kz=pi) = cuy(kz=pi) = cuz(kz=pi) = 0,
   cux(kx=0,ky=0,kz=0) = cuy(kx=0,ky=0,kz=0) = cuz(kx=0,ky=0,kz=0) = 0.
   cu[l][k][j][i] = complex current density for fourier mode (j,k,l)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   requires KNC, cu need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kj, lj, nxvyh;
   float dnx, dny, dnz, dkx, dky, dkz, dky2, dkz2, dkyz2, at1;
   float complex zero, zt1;
   __m512i v_j, v_it;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz, v_dkz2, v_dkyz2;
   __m512 v_dk, v_at1, v_zt1, v_zt2, v_zero, v_one, v_at, v_as;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_one = _mm512_set1_ps(1.0f);
/* calculate transverse part of current */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
   for (l = 1; l < nzh; l++) {
      dkz = dnz*(float) l;
      v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
      lj = nxvyh*l;
      l1 = nxvyh*nz - lj;
      dkz2 = dkz*dkz;
      v_dkz2 = _mm512_set1_ps(dkz2);
/* add kz to gradient operator */
      v_dk  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkz);
      for (k = 1; k < nyh; k++) {
         dky = dny*(float) k;
         v_it = _mm512_set1_epi32(k);
         v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dky = _mm512_mul_ps(v_dny,v_dky);
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         dkyz2 = dky*dky + dkz2;
         v_dkyz2 = _mm512_fmadd_ps(v_dky,v_dky,v_dkz2);
/* add ky to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(3084),v_dky);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          dkx = dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*          at1 = 1.0/(dkx*dkx + dkyz2); */
            v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkyz2);
            v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
            v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*          zt1 = at1*(dkx*cu[4*(j+kj+lj)] + dky*cu[1+4*(j+kj+lj)] */
/*                   + dkz*cu[2+4*(j+kj+lj)]);                     */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
            v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
            }
/*          cu[4*(j+kj+lj)] -= dkx*zt1;   */
/*          cu[1+4*(j+kj+lj)] -= dky*zt1; */
/*          cu[2+4*(j+kj+lj)] -= dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+kj+lj)],v_zt2);
/*          zt1 = at1*(dkx*cu[4*(j+k1+lj)] - dky*cu[1+4*(j+k1+lj)] */
/*                   + dkz*cu[2+4*(j+k1+lj)]);                     */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
            v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(3084),v_zero,
                   v_dk);
            v_zt1 = _mm512_mul_ps(v_as,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
            }
/*          cu[4*(j+k1+lj)] -= dkx*zt1;   */
/*          cu[1+4*(j+k1+lj)] += dky*zt1; */
/*          cu[2+4*(j+k1+lj)] -= dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+k1+lj)],v_zt2);
/*          zt1 = at1*(dkx*cu[4*(j+kj+l1)] + dky*cu[1+4*(j+kj+l1)] */
/*                   - dkz*cu[2+4*(j+kj+l1)]);                     */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
            v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(12336),
                   v_zero,v_dk);
            v_zt1 = _mm512_mul_ps(v_as,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
            }
/*          cu[4*(j+kj+l1)] -= dkx*zt1;   */
/*          cu[1+4*(j+kj+l1)] -= dky*zt1; */
/*          cu[2+4*(j+kj+l1)] += dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+kj+l1)],v_zt2);
/*          zt1 = at1*(dkx*cu[4*(j+k1+l1)] - dky*cu[1+4*(j+k1+l1)] */
/*                   - dkz*cu[2+4*(j+k1+l1)]);                     */
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
            v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(15420),
                   v_zero,v_dk);
            v_zt1 = _mm512_mul_ps(v_as,v_zt2);
            v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
            v_zt1 = _mm512_add_ps(v_at,v_zt1);
            v_at = _mm512_permute4f128_ps(v_zt1,177);
            v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
            }
/*          cu[4*(j+k1+l1)] -= dkx*zt1;   */
/*          cu[1+4*(j+k1+l1)] += dky*zt1; */
/*          cu[2+4*(j+k1+l1)] += dkz*zt1; */
            v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
            _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            dkx = dnx*(float) j;
            at1 = 1.0/(dkx*dkx + dkyz2);
            zt1 = at1*(dkx*cu[4*(j+kj+lj)] + dky*cu[1+4*(j+kj+lj)]
                     + dkz*cu[2+4*(j+kj+lj)]);
            cu[4*(j+kj+lj)] -= dkx*zt1;
            cu[1+4*(j+kj+lj)] -= dky*zt1;
            cu[2+4*(j+kj+lj)] -= dkz*zt1;
            zt1 = at1*(dkx*cu[4*(j+k1+lj)] - dky*cu[1+4*(j+k1+lj)]
                     + dkz*cu[2+4*(j+k1+lj)]);
            cu[4*(j+k1+lj)] -= dkx*zt1;
            cu[1+4*(j+k1+lj)] += dky*zt1;
            cu[2+4*(j+k1+lj)] -= dkz*zt1;
            zt1 = at1*(dkx*cu[4*(j+kj+l1)] + dky*cu[1+4*(j+kj+l1)]
                     - dkz*cu[2+4*(j+kj+l1)]);
            cu[4*(j+kj+l1)] -= dkx*zt1;
            cu[1+4*(j+kj+l1)] -= dky*zt1;
            cu[2+4*(j+kj+l1)] += dkz*zt1;
            zt1 = at1*(dkx*cu[4*(j+k1+l1)] - dky*cu[1+4*(j+k1+l1)]
                     - dkz*cu[2+4*(j+k1+l1)]);
            cu[4*(j+k1+l1)] -= dkx*zt1;
            cu[1+4*(j+k1+l1)] += dky*zt1;
            cu[2+4*(j+k1+l1)] += dkz*zt1;
         }
      }
/* mode numbers kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         dky = dny*(float) k;
         at1 = 1.0/(dky*dky + dkz2);
         zt1 = at1*(dky*cu[1+4*(kj+lj)] + dkz*cu[2+4*(kj+lj)]);
         cu[1+4*(kj+lj)] -= dky*zt1;
         cu[2+4*(kj+lj)] -= dkz*zt1;
         cu[4*(k1+lj)] = zero;
         cu[1+4*(k1+lj)] = zero;
         cu[2+4*(k1+lj)] = zero;
         zt1 = at1*(dky*cu[1+4*(kj+l1)] - dkz*cu[2+4*(kj+l1)]);
         cu[1+4*(kj+l1)] -= dky*zt1;
         cu[2+4*(kj+l1)] += dkz*zt1;
         cu[4*(k1+l1)] = zero;
         cu[1+4*(k1+l1)] = zero;
         cu[2+4*(k1+l1)] = zero;
      }
/* mode numbers ky = 0, ny/2 */
      k1 = nxvh*nyh;
/* add ky to gradient operator */
      v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(3084),v_zero);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*       at1 = 1.0/(dkx*dkx + dkz2); */
         v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkz2);
         v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*       zt1 = at1*(dkx*cu[4*(j+lj)] + dkz*cu[2+4*(j+lj)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
         v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
         }
/*       cu[4*(j+lj)] -= dkx*zt1;   */
/*       cu[2+4*(j+lj)] -= dkz*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+lj)],v_zt2);
/*       cu[4*(j+k1+lj)] = zero;    */
/*       cu[1+4*(j+k1+lj)] = zero; */
/*       cu[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+k1+lj)],v_zero);
/*       zt1 = at1*(dkx*cu[4*(j+l1)] - dkz*cu[2+4*(j+l1)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
         v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(12336),
                v_zero,v_dk);
         v_zt1 = _mm512_mul_ps(v_as,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
         }
/*       cu[4*(j+l1)] -= dkx*zt1;   */
/*       cu[2+4*(j+l1)] += dkz*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+l1)],v_zt2);
/*       cu[4*(j+k1+l1)] = zero;   */
/*       cu[1+4*(j+k1+l1)] = zero; */
/*       cu[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dkz2);
         zt1 = at1*(dkx*cu[4*(j+lj)] + dkz*cu[2+4*(j+lj)]);
         cu[4*(j+lj)] -= dkx*zt1;
         cu[2+4*(j+lj)] -= dkz*zt1;
         cu[4*(j+k1+lj)] = zero;
         cu[1+4*(j+k1+lj)] = zero;
         cu[2+4*(j+k1+lj)] = zero;
         zt1 = at1*(dkx*cu[4*(j+l1)] - dkz*cu[2+4*(j+l1)]);
         cu[4*(j+l1)] -= dkx*zt1;
         cu[2+4*(j+l1)] += dkz*zt1;
         cu[4*(j+k1+l1)] = zero;
         cu[1+4*(j+k1+l1)] = zero;
         cu[2+4*(j+k1+l1)] = zero;
      }
/* mode numbers kx = 0, nx/2 */
      cu[2+4*lj] = zero;
      cu[4*(k1+lj)] = zero;
      cu[1+4*(k1+lj)] = zero;
      cu[2+4*(k1+lj)] = zero;
      cu[4*l1] = zero;
      cu[1+4*l1] = zero;
      cu[2+4*l1] = zero;
      cu[4*(k1+l1)] = zero;
      cu[1+4*(k1+l1)] = zero;
      cu[2+4*(k1+l1)] = zero;
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      dky2 = dky*dky;
      v_dkyz2 = _mm512_mul_ps(v_dky,v_dky);
/* add ky to gradient operator */
      v_dk  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/*       at1 = 1.0/(dkx*dkx + dky2); */
         v_at1 = _mm512_fmadd_ps(v_dkx,v_dkx,v_dkyz2);
         v_at1 = _mm512_div_ps(v_one,v_at1);
/* add kx to gradient operator */
         v_dk  = _mm512_mask_mov_ps(v_dk,_mm512_int2mask(771),v_dkx);
/*       zt1 = at1*(dkx*cu[4*(j+kj)] + dky*cu[1+4*(j+kj)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt1 = _mm512_mul_ps(v_dk,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
         }
/*       cu[4*(j+kj)] -= dkx*zt1;   */
/*       cu[1+4*(j+kj)] -= dky*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_dk,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+kj)],v_zt2);
/*       zt1 = at1*(dkx*cu[4*(j+k1)]- dky*cu[1+4*(j+k1)]); */
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_as = _mm512_mask_sub_ps(v_dk,_mm512_int2mask(3084),v_zero,
                v_dk);
         v_zt1 = _mm512_mul_ps(v_as,v_zt2);
         v_at = (__m512)_mm512_shuffle_epi32((__m512i)v_zt1,78);
         v_zt1 = _mm512_add_ps(v_at,v_zt1);
         v_at = _mm512_permute4f128_ps(v_zt1,177);
         v_zt1 = _mm512_mul_ps(v_at1,_mm512_add_ps(v_at,v_zt1));
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);;
         }
/*       cu[4*(j+k1)] -= dkx*zt1;   */
/*       cu[1+4*(j+k1)] += dky*zt1; */
         v_zt2 = _mm512_sub_ps(v_zt2,_mm512_mul_ps(v_as,v_zt1));
         _mm512_store_ps((float *)&cu[4*(j+k1)],v_zt2);
/*       cu[4*(j+kj+l1)] = zero;   */
/*       cu[1+4*(j+kj+l1)] = zero; */
/*       cu[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+kj+l1)],v_zero);
/*       cu[4*(j+k1+l1)] = zero;   */
/*       cu[1+4*(j+k1+l1)] = zero; */
/*       cu[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;
         at1 = 1.0/(dkx*dkx + dky2);
         zt1 = at1*(dkx*cu[4*(j+kj)] + dky*cu[1+4*(j+kj)]);
         cu[4*(j+kj)] -= dkx*zt1;
         cu[1+4*(j+kj)] -= dky*zt1;
         zt1 = at1*(dkx*cu[4*(j+k1)]- dky*cu[1+4*(j+k1)]);
         cu[4*(j+k1)] -= dkx*zt1;
         cu[1+4*(j+k1)] += dky*zt1;
         cu[4*(j+kj+l1)] = zero;
         cu[1+4*(j+kj+l1)] = zero;
         cu[2+4*(j+kj+l1)] = zero;
         cu[4*(j+k1+l1)] = zero;
         cu[1+4*(j+k1+l1)] = zero;
         cu[2+4*(j+k1+l1)] = zero;
      }

   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      cu[1+4*kj] = zero;
      cu[4*k1] = zero;
      cu[1+4*k1] = zero;
      cu[2+4*k1] = zero;
      cu[4*(kj+l1)] = zero;
      cu[1+4*(kj+l1)] = zero;
      cu[2+4*(kj+l1)] = zero;
      cu[4*(k1+l1)] = zero;
      cu[1+4*(k1+l1)] = zero;
      cu[2+4*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
      v_zt2 = _mm512_load_ps((float *)&cu[4*j]);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(255),v_zero);
      }
/*    cu[4*j] = zero;      */
      v_zt2 = _mm512_mask_mov_ps(v_zt2,_mm512_int2mask(771),v_zero);
      _mm512_store_ps((float *)&cu[4*j],v_zt2);
/*    cu[4*(j+k1)] = zero; */
/*    cu[1+4*(j+k1)] = zero; */
/*    cu[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+k1)],v_zero);
/*    cu[4*(j+l1)] = zero;   */
/*    cu[1+4*(j+l1)] = zero; */
/*    cu[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+l1)],v_zero);
/*    cu[4*(j+k1+l1)] = zero;   */
/*    cu[1+4*(j+k1+l1)] = zero; */
/*    cu[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&cu[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
      cu[4*j] = zero;
      cu[4*(j+k1)] = zero;
      cu[1+4*(j+k1)] = zero;
      cu[2+4*(j+k1)] = zero;
      cu[4*(j+l1)] = zero;
      cu[1+4*(j+l1)] = zero;
      cu[2+4*(j+l1)] = zero;
      cu[4*(j+k1+l1)] = zero;
      cu[1+4*(j+k1+l1)] = zero;
      cu[2+4*(j+k1+l1)] = zero;
   }
   cu[0] = zero;
   cu[1] = zero;
   cu[2] = zero;
   cu[4*k1] = zero;
   cu[1+4*k1] = zero;
   cu[2+4*k1] = zero;
   cu[4*l1] = zero;
   cu[1+4*l1] = zero;
   cu[2+4*l1] = zero;
   cu[4*(k1+l1)] = zero;
   cu[1+4*(k1+l1)] = zero;
   cu[2+4*(k1+l1)] = zero;
   return;
}

/*--------------------------------------------------------------------*/
void ckncibpois33(float complex cu[], float complex bxyz[],
                  float complex ffc[], float ci, float *wm, int nx,
                  int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                  int nyhd, int nzhd) {
/* this subroutine solves 3d poisson's equation in fourier space for
   magnetic field with periodic boundary conditions.
   input: cu,ffc,ci,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
   output: bxyz, wm
   approximate flop count is:
   193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the magnetic field is calculated using the equations:
   bx[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (ky*cuz[kz][ky][kx]-kz*cuy[kz][ky][kx]),
   by[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kz*cux[kz][ky][kx]-kx*cuz[kz][ky][kx]),
   bz[kz][ky][kx] = ci*ci*sqrt(-1)*g[kz][ky][kx]*
                  (kx*cuy[kz][ky][kx]-ky*cux[kz][ky][kx]),
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
   j,k,l = fourier mode numbers,
   g[kz][ky][kx] = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
   bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
   bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
   bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
   bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
   cu[l][k][j][i] = complex current density for fourier mode (j,k,l)
   bxyz[l][k][j][i] = i component of complex magnetic field
   all for fourier mode (j,k,l)
   aimag(ffc(j,k,l)) = finite-size particle shape factor s
   for fourier mode (j,k,l)
   real(ffc(j,k,l)) = potential green's function g
   for fourier mode (j,k,l)
   ci = reciprocal of velocity of light
   magnetic field energy is also calculated, using
   wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
      |cu[kz][ky][kx]*s[kz][ky][kx]|**2)
   this expression is valid only if the current is divergence-free
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   nxhd = dimension of form factor array, must be >= nxh
   nyhd = second dimension of form factor array, must be >= nyh
   nzhd = third dimension of form factor array, must be >= nzh
   requires KNC, cu, bxyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu, bxyz need to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh;
   float dnx, dny, dnz, dky, dkz, ci2, at1, at2, at3, at4;
   float complex zero, zt1, zt2, zt3;
   double wp, d0;
   __m512i v_j, v_it, v_n, v_m;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz, v_ci2;
   __m512 v_dk1, v_dk2, v_at1, v_at2, v_at3, v_at4, v_zero;
   __m512 v_zt1, v_zt2, v_zt3, v_zt4;
   __m512d v_wp, v_d;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   zero = 0.0 + 0.0*_Complex_I;
   ci2 = ci*ci;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_n =  _mm512_set_epi32(15,14,11,10,9,8,13,12,7,6,3,2,1,0,5,4);
   v_m =  _mm512_set_epi32(15,14,9,8,13,12,11,10,7,6,1,0,5,4,3,2);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_ci2 = _mm512_set1_ps(ci2);
/* calculate magnetic field and sum field energy */
   wp = 0.0;
   v_wp = _mm512_set1_pd(0.0);
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
   for (l = 1; l < nzh; l++) {
      dkz = dnz*(float) l;
      v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
      ll = nxyhd*l;
      lj = nxvyh*l;
      l1 = nxvyh*nz - lj;
/* add kz to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dkz);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkz);
      for (k = 1; k < nyh; k++) {
         dky = dny*(float) k;
         v_it = _mm512_set1_epi32(k);
         v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dky = _mm512_mul_ps(v_dny,v_dky);
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* add ky to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),
                  v_dky);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),
                  v_dky);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = ci2*crealf(ffc[j+kk+ll]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(43690),(__m512i)v_at4,177);
            v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*          at2 = at1*dnx*(float) j; */
/*          at3 = dky*at1;           */
/*          at4 = dkz*at1;           */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                     v_dkx);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                     v_dkx);
/* normalize curl operators */
            v_at2 = _mm512_mul_ps(v_at1,v_dk1);
            v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*          at1 = at1*cimagf(ffc[j+kk+ll]); */
            v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(21845),(__m512i)v_at4,177);
            v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*          zt1 = -cimagf(cu[2+4*(j+kj+lj)])              */
/*               + crealf(cu[2+4*(j+kj+lj)])*_Complex_I;/ */
/*          zt2 = -cimagf(cu[1+4*(j+kj+lj)])              */
/*               + crealf(cu[1+4*(j+kj+lj)])*_Complex_I;  */
/*          zt3 = -cimagf(cu[4*(j+kj+lj)])                */
/*               + crealf(cu[4*(j+kj+lj)])*_Complex_I;    */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                    v_zero,v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          bxyz[4*(j+kj+lj)] = at3*zt1 - at4*zt2;   */
/*          bxyz[1+4*(j+kj+lj)] = at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+kj+lj)] = at2*zt2 - at3*zt3; */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+kj+lj)],v_zt1);
/*          wp += at1*(cu[4*(j+kj+lj)]*conjf(cu[4*(j+kj+lj)]) */
/*             + cu[1+4*(j+kj+lj)]*conjf(cu[1+4*(j+kj+lj)])   */
/*             + cu[2+4*(j+kj+lj)]*conjf(cu[2+4*(j+kj+lj)])); */
            v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3));
/*          zt1 = -cimagf(cu[2+4*(j+k1+lj)])             */
/*               + crealf(cu[2+4*(j+k1+lj)])*_Complex_I; */
/*          zt2 = -cimagf(cu[1+4*(j+k1+lj)])             */
/*               + crealf(cu[1+4*(j+k1+lj)])*_Complex_I; */
/*          zt3 = -cimagf(cu[4*(j+k1+lj)])               */
/*               + crealf(cu[4*(j+k1+lj)])*_Complex_I;   */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                    v_zero,v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(12336),
                    v_zero,v_at2);
            v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(771),
                    v_zero,v_at3);
/*          bxyz[4*(j+k1+lj)] = -at3*zt1 - at4*zt2;  */
/*          bxyz[1+4*(j+k1+lj)] = at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+k1+lj)] = at2*zt2 + at3*zt3; */
            v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zt1);
/*          wp += at1*(cu[4*(j+k1+lj)]*conjf(cu[4*(j+k1+lj)]) */
/*             + cu[1+4*(j+k1+lj)]*conjf(cu[1+4*(j+k1+lj)])   */
/*             + cu[2+4*(j+k1+lj)]*conjf(cu[2+4*(j+k1+lj)])); */
            v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/*          zt1 = -cimagf(cu[2+4*(j+kj+l1)])             */
/*               + crealf(cu[2+4*(j+kj+l1)])*_Complex_I; */
/*          zt2 = -cimagf(cu[1+4*(j+kj+l1)])             */
/*               + crealf(cu[1+4*(j+kj+l1)])*_Complex_I; */
/*          zt3 = -cimagf(cu[4*(j+kj+l1)])               */
/*               + crealf(cu[4*(j+kj+l1)])*_Complex_I;   */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                    v_zero,v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(771),
                    v_zero,v_at2);
            v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3084),
                    v_zero,v_at3);
/*          bxyz[4*(j+kj+l1)] = at3*zt1 + at4*zt2;    */
/*          bxyz[1+4*(j+kj+l1)] = -at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+kj+l1)] = at2*zt2 - at3*zt3;  */
            v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zt1);
/*          wp += at1*(cu[4*(j+kj+l1)]*conjf(cu[4*(j+kj+l1)]) */
/*             + cu[1+4*(j+kj+l1)]*conjf(cu[1+4*(j+kj+l1)])   */
/*             + cu[2+4*(j+kj+l1)]*conjf(cu[2+4*(j+kj+l1)])); */
            v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/*          zt1 = -cimagf(cu[2+4*(j+k1+l1)])             */
/*               + crealf(cu[2+4*(j+k1+l1)])*_Complex_I; */
/*          zt2 = -cimagf(cu[1+4*(j+k1+l1)])             */
/*               + crealf(cu[1+4*(j+k1+l1)])*_Complex_I; */
/*          zt3 = -cimagf(cu[4*(j+k1+l1)])               */
/*               + crealf(cu[4*(j+k1+l1)])*_Complex_I;   */
            v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),
                    v_zero,v_zt3);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(13107),
                    v_zero,v_at2);
            v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3855),
                    v_zero,v_at3);
/*          bxyz[4*(j+k1+l1)] = -at3*zt1 + at4*zt2;   */
/*          bxyz[1+4*(j+k1+l1)] = -at4*zt3 - at2*zt1; */
/*          bxyz[2+4*(j+k1+l1)] = at2*zt2 + at3*zt3;  */
            v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                       v_zero);
               v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                       v_zero);
            }
            _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zt1);
/*          wp += at1*(cu[4*(j+k1+l1)]*conjf(cu[4*(j+k1+l1)]) */
/*             + cu[1+4*(j+k1+l1)]*conjf(cu[1+4*(j+k1+l1)])   */
/*             + cu[2+4*(j+k1+l1)]*conjf(cu[2+4*(j+k1+l1)])); */
            v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                    _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = ci2*crealf(ffc[j+kk+ll]);
            at2 = at1*dnx*(float) j;
            at3 = dky*at1;
            at4 = dkz*at1;
            at1 = at1*cimagf(ffc[j+kk+ll]);
            zt1 = -cimagf(cu[2+4*(j+kj+lj)])
                 + crealf(cu[2+4*(j+kj+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+kj+lj)])
                 + crealf(cu[1+4*(j+kj+lj)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+kj+lj)])
                 + crealf(cu[4*(j+kj+lj)])*_Complex_I;
            bxyz[4*(j+kj+lj)] = at3*zt1 - at4*zt2;
            bxyz[1+4*(j+kj+lj)] = at4*zt3 - at2*zt1;
            bxyz[2+4*(j+kj+lj)] = at2*zt2 - at3*zt3;
            zt1 = -cimagf(cu[2+4*(j+k1+lj)])
                 + crealf(cu[2+4*(j+k1+lj)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+k1+lj)])
                 + crealf(cu[1+4*(j+k1+lj)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+k1+lj)])
                 + crealf(cu[4*(j+k1+lj)])*_Complex_I;
            bxyz[4*(j+k1+lj)] = -at3*zt1 - at4*zt2;
            bxyz[1+4*(j+k1+lj)] = at4*zt3 - at2*zt1;
            bxyz[2+4*(j+k1+lj)] = at2*zt2 + at3*zt3;
            zt1 = -cimagf(cu[2+4*(j+kj+l1)])
                 + crealf(cu[2+4*(j+kj+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+kj+l1)])
                 + crealf(cu[1+4*(j+kj+l1)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+kj+l1)])
                 + crealf(cu[4*(j+kj+l1)])*_Complex_I;
            bxyz[4*(j+kj+l1)] = at3*zt1 + at4*zt2;
            bxyz[1+4*(j+kj+l1)] = -at4*zt3 - at2*zt1;
            bxyz[2+4*(j+kj+l1)] = at2*zt2 - at3*zt3;
            zt1 = -cimagf(cu[2+4*(j+k1+l1)])
                 + crealf(cu[2+4*(j+k1+l1)])*_Complex_I;
            zt2 = -cimagf(cu[1+4*(j+k1+l1)])
                 + crealf(cu[1+4*(j+k1+l1)])*_Complex_I;
            zt3 = -cimagf(cu[4*(j+k1+l1)])
                 + crealf(cu[4*(j+k1+l1)])*_Complex_I;
            bxyz[4*(j+k1+l1)] = -at3*zt1 + at4*zt2;
            bxyz[1+4*(j+k1+l1)] = -at4*zt3 - at2*zt1;
            bxyz[2+4*(j+k1+l1)] = at2*zt2 + at3*zt3;
            at1 = at1*(cu[4*(j+kj+lj)]*conjf(cu[4*(j+kj+lj)])
                + cu[1+4*(j+kj+lj)]*conjf(cu[1+4*(j+kj+lj)])
                + cu[2+4*(j+kj+lj)]*conjf(cu[2+4*(j+kj+lj)])
                + cu[4*(j+k1+lj)]*conjf(cu[4*(j+k1+lj)])
                + cu[1+4*(j+k1+lj)]*conjf(cu[1+4*(j+k1+lj)])
                + cu[2+4*(j+k1+lj)]*conjf(cu[2+4*(j+k1+lj)])
                + cu[4*(j+kj+l1)]*conjf(cu[4*(j+kj+l1)])
                + cu[1+4*(j+kj+l1)]*conjf(cu[1+4*(j+kj+l1)])
                + cu[2+4*(j+kj+l1)]*conjf(cu[2+4*(j+kj+l1)])
                + cu[4*(j+k1+l1)]*conjf(cu[4*(j+k1+l1)])
                + cu[1+4*(j+k1+l1)]*conjf(cu[1+4*(j+k1+l1)])
                + cu[2+4*(j+k1+l1)]*conjf(cu[2+4*(j+k1+l1)]));
            wp += (double) at1;
         }
      }
/* mode numbers kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         at1 = ci2*crealf(ffc[kk+ll]);
         at3 = at1*dny*(float) k;
         at4 = dkz*at1;
         at1 = at1*cimagf(ffc[kk+ll]);
         zt1 = -cimagf(cu[2+4*(kj+lj)])
              + crealf(cu[2+4*(kj+lj)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(kj+lj)])
              + crealf(cu[1+4*(kj+lj)])*_Complex_I;
         zt3 = -cimagf(cu[4*(kj+lj)])
              + crealf(cu[4*(kj+lj)])*_Complex_I;
         bxyz[4*(kj+lj)] = at3*zt1 - at4*zt2;
         bxyz[1+4*(kj+lj)] = at4*zt3;
         bxyz[2+4*(kj+lj)] = -at3*zt3;
         bxyz[4*(k1+lj)] = zero;
         bxyz[1+4*(k1+lj)] = zero;
         bxyz[2+4*(k1+lj)] = zero;
         zt1 = -cimagf(cu[2+4*(kj+l1)])
              + crealf(cu[2+4*(kj+l1)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(kj+l1)])
              + crealf(cu[1+4*(kj+l1)])*_Complex_I;
         zt3 = -cimagf(cu[4*(kj+l1)])
              + crealf(cu[4*(kj+l1)])*_Complex_I;
         bxyz[4*(kj+l1)] = at3*zt1 + at4*zt2;
         bxyz[1+4*(kj+l1)] = -at4*zt3;
         bxyz[2+4*(kj+l1)] = -at3*zt3;
         bxyz[4*(k1+l1)] = zero;
         bxyz[1+4*(k1+l1)]  = zero;
         bxyz[2+4*(k1+l1)] = zero;
         at1 = at1*(cu[4*(kj+lj)]*conjf(cu[4*(kj+lj)])
             + cu[1+4*(kj+lj)]*conjf(cu[1+4*(kj+lj)])
             + cu[2+4*(kj+lj)]*conjf(cu[2+4*(kj+lj)])
             + cu[4*(kj+l1)]*conjf(cu[4*(kj+l1)])
             + cu[1+4*(kj+l1)]*conjf(cu[1+4*(kj+l1)])
             + cu[2+4*(kj+l1)]*conjf(cu[2+4*(kj+l1)]));
         wp += (double) at1;
      }
/* mode numbers ky = 0, ny/2 */
      k1 = nxvh*nyh;
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),v_zero);
      v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),v_zero);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = ci2*crealf(ffc[j+ll]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j+ll]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
                 (float *)&ffc[j+ll+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                 _mm512_int2mask(43690),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*       at2 = at1*dnx*(float) j;     */
/*       at4 = dkz*at1;               */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/* normalize curl operators */
         v_at2 = _mm512_mul_ps(v_at1,v_dk1);
         v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*       at1 = at1*cimagf(ffc[j+ll]); */
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                 _mm512_int2mask(21845),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*       zt1 = -cimagf(cu[2+4*(j+lj)])             */
/*            + crealf(cu[2+4*(j+lj)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+lj)])             */
/*            + crealf(cu[1+4*(j+lj)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+lj)])               */
/*            + crealf(cu[4*(j+lj)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       bxyz[4*(j+lj)] = -at4*zt2;            */
/*       bxyz[1+4*(j+lj)] = at4*zt3 - at2*zt1; */
/*       bxyz[2+4*(j+lj)] = at2*zt2;           */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+lj)],v_zt1);
/*       wp += at1*(cu[4*(j+lj)]*conjf(cu[4*(j+lj)]) */
/*          + cu[1+4*(j+lj)]*conjf(cu[1+4*(j+lj)])   */
/*          + cu[2+4*(j+lj)]*conjf(cu[2+4*(j+lj)])   */
         v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3));
/*       bxyz[4*(j+k1+lj)] = zero;   */
/*       bxyz[1+4*(j+k1+lj)] = zero; */
/*       bxyz[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zero);
/*       zt1 = -cimagf(cu[2+4*(j+l1)])             */
/*            + crealf(cu[2+4*(j+l1)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+l1)])             */
/*            + crealf(cu[1+4*(j+l1)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+l1)])               */
/*            + crealf(cu[4*(j+l1)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(771),v_zero,
                 v_at2);
         v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(3084),v_zero,
                 v_at3);
/*       bxyz[4*(j+l1)] = at4*zt2;              */
/*       bxyz[1+4*(j+l1)] = -at4*zt3 - at2*zt1; */
/*       bxyz[2+4*(j+l1)] = at2*zt2;            */
         v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zt1);
/*       wp += at1*(cu[4*(j+l1)]*conjf(cu[4*(j+l1)]) */ 
/*          + cu[1+4*(j+l1)]*conjf(cu[1+4*(j+l1)])   */
/*          + cu[2+4*(j+l1)]*conjf(cu[2+4*(j+l1)])); */
         v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+k1+l1)] = zero;   */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+ll]);
         at2 = at1*dnx*(float) j;  
         at4 = dkz*at1;
         at1 = at1*cimagf(ffc[j+ll]);
         zt1 = -cimagf(cu[2+4*(j+lj)])
              + crealf(cu[2+4*(j+lj)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+lj)])
              + crealf(cu[1+4*(j+lj)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+lj)])
              + crealf(cu[4*(j+lj)])*_Complex_I;
         bxyz[4*(j+lj)] = -at4*zt2;
         bxyz[1+4*(j+lj)] = at4*zt3 - at2*zt1;
         bxyz[2+4*(j+lj)] = at2*zt2;
         bxyz[4*(j+k1+lj)] = zero;
         bxyz[1+4*(j+k1+lj)] = zero;
         bxyz[2+4*(j+k1+lj)] = zero;
         zt1 = -cimagf(cu[2+4*(j+l1)])
              + crealf(cu[2+4*(j+l1)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+l1)])
              + crealf(cu[1+4*(j+l1)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+l1)])
              + crealf(cu[4*(j+l1)])*_Complex_I;
         bxyz[4*(j+l1)] = at4*zt2;
         bxyz[1+4*(j+l1)] = -at4*zt3 - at2*zt1;
         bxyz[2+4*(j+l1)] = at2*zt2;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(cu[4*(j+lj)]*conjf(cu[4*(j+lj)])
             + cu[1+4*(j+lj)]*conjf(cu[1+4*(j+lj)])
             + cu[2+4*(j+lj)]*conjf(cu[2+4*(j+lj)])
             + cu[4*(j+l1)]*conjf(cu[4*(j+l1)])
             + cu[1+4*(j+l1)]*conjf(cu[1+4*(j+l1)])
             + cu[2+4*(j+l1)]*conjf(cu[2+4*(j+l1)]));
         wp += (double) at1;
      }
/* mode numbers kx = 0, nx/2 */
      at1 = ci2*crealf(ffc[ll]);
      at4 = dkz*at1;
      at1 = at1*cimagf(ffc[ll]);
      zt2 = -cimagf(cu[1+4*(lj)]) + crealf(cu[1+4*(lj)])*_Complex_I;
      zt3 = -cimagf(cu[4*(lj)]) + crealf(cu[4*(lj)])*_Complex_I;
      bxyz[4*lj] = -at4*zt2;
      bxyz[1+4*lj] = at4*zt3;
      bxyz[2+4*lj] = zero;
      bxyz[4*(k1+lj)] = zero;
      bxyz[1+4*(k1+lj)] = zero;
      bxyz[2+4*(k1+lj)] = zero;
      bxyz[4*l1] = zero;
      bxyz[1+4*l1] = zero;
      bxyz[2+4*l1] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(cu[4*lj]*conjf(cu[4*lj])
          + cu[1+4*lj]*conjf(cu[1+4*lj])
          + cu[2+4*lj]*conjf(cu[2+4*lj]));
      wp += (double) at1;
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dky);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = ci2*crealf(ffc[j+kk]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j+kk]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
                 (float *)&ffc[j+kk+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                 _mm512_int2mask(43690),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*       at2 = at1*dnx*(float) j; */
/*       at3 = dky*at1;           */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/* normalize curl operators */
         v_at2 = _mm512_mul_ps(v_at1,v_dk1);
         v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*       at1 = at1*cimagf(ffc[j+kk]); */
         v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
                    _mm512_int2mask(21845),(__m512i)v_at4,177);
         v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*       zt1 = -cimagf(cu[2+4*(j+kj)])             */
/*            + crealf(cu[2+4*(j+kj)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+kj)])             */
/*            + crealf(cu[1+4*(j+kj)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+kj)])               */
/*            + crealf(cu[4*(j+kj)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       bxyz[4*(j+kj)] = at3*zt1;             */
/*       bxyz[1+4*(j+kj)] = -at2*zt1;          */
/*       bxyz[2+4*(j+kj)] = at2*zt2 - at3*zt3; */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+kj)],v_zt1);
/*       wp += at1*(cu[4*(j+kj)]*conjf(cu[4*(j+kj)]) */
/*          + cu[1+4*(j+kj)]*conjf(cu[1+4*(j+kj)])   */
/*          + cu[2+4*(j+kj)]*conjf(cu[2+4*(j+kj)])); */
         v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3));
/*       zt1 = -cimagf(cu[2+4*(j+k1)])             */
/*            + crealf(cu[2+4*(j+k1)])*_Complex_I; */
/*       zt2 = -cimagf(cu[1+4*(j+k1)])             */
/*            + crealf(cu[1+4*(j+k1)])*_Complex_I; */
/*       zt3 = -cimagf(cu[4*(j+k1)])               */
/*            + crealf(cu[4*(j+k1)])*_Complex_I;   */
         v_zt3 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
                 v_zt3);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_zt1 = _mm512_mask_sub_ps(v_at2,_mm512_int2mask(12336),v_zero,
                 v_at2);
         v_zt2 = _mm512_mask_sub_ps(v_at3,_mm512_int2mask(771),v_zero,
                 v_at3);
/*       bxyz[4*(j+k1)] = -at3*zt1;            */
/*       bxyz[1+4*(j+k1)] = -at2*zt1;          */
/*       bxyz[2+4*(j+k1)] = at2*zt2 + at3*zt3; */
         v_zt1 = _mm512_mul_ps(v_zt1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_zt2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),
                    v_zero);
            v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),
                    v_zero);
         }
         _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zt1);
/*       wp += at1*(cu[4*(j+k1)]*conjf(cu[4*(j+k1)]) */
/*          + cu[1+4*(j+k1)]*conjf(cu[1+4*(j+k1)])   */
/*          + cu[2+4*(j+k1)]*conjf(cu[2+4*(j+k1)])); */
         v_zt4 = _mm512_fmadd_ps(v_at1,_mm512_mask_mul_ps(v_zero,
                 _mm512_int2mask(16191),v_zt3,v_zt3),v_zt4);
/* convert to double precision before accumulating */
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+kj+l1)] = zero;   */
/*       bxyz[1+4*(j+kj+l1)] = zero; */
/*       bxyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zero);
/*       bxyz[4*(j+k1+l1)] = zero;   */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = ci2*crealf(ffc[j+kk]);
         at2 = at1*dnx*(float) j;
         at3 = dky*at1;
         at1 = at1*cimagf(ffc[j+kk]);
         zt1 = -cimagf(cu[2+4*(j+kj)])
              + crealf(cu[2+4*(j+kj)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+kj)])
              + crealf(cu[1+4*(j+kj)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+kj)])
              + crealf(cu[4*(j+kj)])*_Complex_I;
         bxyz[4*(j+kj)] = at3*zt1;
         bxyz[1+4*(j+kj)] = -at2*zt1;
         bxyz[2+4*(j+kj)] = at2*zt2 - at3*zt3;
         zt1 = -cimagf(cu[2+4*(j+k1)])
              + crealf(cu[2+4*(j+k1)])*_Complex_I;
         zt2 = -cimagf(cu[1+4*(j+k1)])
              + crealf(cu[1+4*(j+k1)])*_Complex_I;
         zt3 = -cimagf(cu[4*(j+k1)])
              + crealf(cu[4*(j+k1)])*_Complex_I;
         bxyz[4*(j+k1)] = -at3*zt1;
         bxyz[1+4*(j+k1)] = -at2*zt1;
         bxyz[2+4*(j+k1)] = at2*zt2 + at3*zt3;
         bxyz[4*(j+kj+l1)] = zero;
         bxyz[1+4*(j+kj+l1)] = zero;
         bxyz[2+4*(j+kj+l1)] = zero;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         at1 = at1*(cu[4*(j+kj)]*conjf(cu[4*(j+kj)])
             + cu[1+4*(j+kj)]*conjf(cu[1+4*(j+kj)])
             + cu[2+4*(j+kj)]*conjf(cu[2+4*(j+kj)])
             + cu[4*(j+k1)]*conjf(cu[4*(j+k1)])
             + cu[1+4*(j+k1)]*conjf(cu[1+4*(j+k1)])
             + cu[2+4*(j+k1)]*conjf(cu[2+4*(j+k1)]));
         wp += (double) at1;
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      at1 = ci2*crealf(ffc[kk]);
      at3 = at1*dny*(float) k;
      at1 = at1*cimagf(ffc[kk]);
      zt1 = -cimagf(cu[2+4*(kj)]) + crealf(cu[2+4*(kj)])*_Complex_I;
      zt3 = -cimagf(cu[4*(kj)]) + crealf(cu[4*(kj)])*_Complex_I;
      bxyz[4*kj] = at3*zt1;
      bxyz[1+4*kj] = zero;
      bxyz[2+4*kj] = -at3*zt3;
      bxyz[4*k1] = zero;
      bxyz[1+4*k1] = zero;
      bxyz[2+4*k1] = zero;
      bxyz[4*(kj+l1)] = zero;
      bxyz[1+4*(kj+l1)] = zero;
      bxyz[2+4*(kj+l1)] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      at1 = at1*(cu[4*kj]*conjf(cu[4*kj])
          + cu[1+4*kj]*conjf(cu[1+4*kj])
          + cu[2+4*kj]*conjf(cu[2+4*kj]));
      wp += (double) at1;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
/*    at1 = ci2*crealf(ffc[j]); */
      v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
              (float *)&ffc[j]);
      v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,_mm512_int2mask(15),
              (float *)&ffc[j+8]);
      v_at1 = _mm512_permute4f128_ps(v_at1,0);
      v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
              _mm512_int2mask(13260),(__m512i)v_at1,78);
      v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
               _mm512_int2mask(43690),(__m512i)v_at4,177);
      v_at1 = _mm512_mul_ps(v_ci2,v_at1);
/*    at2 = at1*dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkx);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkx);
/* normalize curl operators */
      v_at2 = _mm512_mul_ps(v_at1,v_dk1);
      v_at3 = _mm512_mul_ps(v_at1,v_dk2);
/*    at1 = at1*cimagf(ffc[j]); */
      v_at4 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at4,
              _mm512_int2mask(21845),(__m512i)v_at4,177);
      v_at1 = _mm512_mul_ps(v_at1,v_at4);
/*    zt1 = -cimagf(cu[2+4*j]) + crealf(cu[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(cu[1+4*j]) + crealf(cu[1+4*j])*_Complex_I; */
      v_zt3 = _mm512_load_ps((float *)&cu[4*j]);
      v_zt3 = _mm512_mask_sub_ps(v_zt3,_mm512_int2mask(43690),v_zero,
              v_zt3);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    bxyz[4*j] = zero;       */
/*    bxyz[1+4*j] = -at2*zt1; */
/*    bxyz[2+4*j] = at2*zt2;  */
      v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
      v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_sub_ps(v_zt1,v_zt2);
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt1 = _mm512_mask_mov_ps(v_zt1,_mm512_int2mask(255),v_zero);
         v_zt3 = _mm512_mask_mov_ps(v_zt3,_mm512_int2mask(255),v_zero);
      }
      _mm512_store_ps((float *)&bxyz[4*j],v_zt1);
/*    wp += at1*(cu[4*j]*conjf(cu[4*j]) */
/*       + cu[1+4*j]*conjf(cu[1+4*j])   */
/*       + cu[2+4*j]*conjf(cu[2+4*j])); */
      v_zt4 = _mm512_mul_ps(v_at1,_mm512_mask_mul_ps(v_zero,
              _mm512_int2mask(16191),v_zt3,v_zt3));
/* convert to double precision before accumulating */
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt4));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt4,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
/*    bxyz[4*(j+k1)] = zero;   */
/*    bxyz[1+4*(j+k1)] = zero; */
/*    bxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zero);
/*    bxyz[4*(j+l1)] = zero;   */
/*    bxyz[1+4*(j+l1)] = zero; */
/*    bxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zero);
/*    bxyz[4*(j+k1+l1)] = zero;    */
/*    bxyz[1+4*(j+k1+l1)] = zero; */
/*    bxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
      at1 = ci2*crealf(ffc[j]);
      at2 = at1*dnx*(float) j;
      at1 = at1*cimagf(ffc[j]);
      zt1 = -cimagf(cu[2+4*j]) + crealf(cu[2+4*j])*_Complex_I;
      zt2 = -cimagf(cu[1+4*j]) + crealf(cu[1+4*j])*_Complex_I;
      bxyz[4*j] = zero;
      bxyz[1+4*j] = -at2*zt1;
      bxyz[2+4*j] = at2*zt2;
      bxyz[4*(j+k1)] = zero;
      bxyz[1+4*(j+k1)] = zero;
      bxyz[2+4*(j+k1)] = zero;
      bxyz[4*(j+l1)] = zero;
      bxyz[1+4*(j+l1)] = zero;
      bxyz[2+4*(j+l1)] = zero;
      bxyz[4*(j+k1+l1)] = zero;
      bxyz[1+4*(j+k1+l1)] = zero;
      bxyz[2+4*(j+k1+l1)] = zero;
      at1 = at1*(cu[4*j]*conjf(cu[4*j])
          + cu[1+4*j]*conjf(cu[1+4*j])
          + cu[2+4*j]*conjf(cu[2+4*j]));
      wp += (double) at1;
   }
   bxyz[0] = zero;
   bxyz[1] = zero;
   bxyz[2] = zero;
   bxyz[4*k1] = zero;
   bxyz[1+4*k1] = zero;
   bxyz[2+4*k1] = zero;
   bxyz[4*l1] = zero;
   bxyz[1+4*l1] = zero;
   bxyz[2+4*l1] = zero;
   bxyz[4*(k1+l1)] = zero;
   bxyz[1+4*(k1+l1)] = zero;
   bxyz[2+4*(k1+l1)] = zero;
   d0 = _mm512_reduce_add_pd(v_wp);
   *wm = (wp + d0)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmaxwel3(float complex exyz[], float complex bxyz[],
                 float complex cu[], float complex ffc[], float ci,
                 float dt, float *wf, float *wm, int nx, int ny, int nz,
                 int nxvh, int nyv, int nzv, int nxhd, int nyhd,
                 int nzhd) {
/* this subroutine solves 3d maxwell's equation in fourier space for
   transverse electric and magnetic fields with periodic boundary
   conditions.
   input: all, output: wf, wm, exyz, bxyz
   approximate flop count is:
   680*nxc*nyc*nzc + 149*(nxc*nyc + nxc*nzc + nyc*nzc)
   plus nxc*nyc*nzc divides
   where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
   the magnetic field is first updated half a step using the equations:
   bx[kz][ky][kx] = bx[kz][ky][kx] - .5*dt*sqrt(-1)*
                    (ky*ez[kz][ky][kx]-kz*ey[kz][ky][kx])
   by[kz][ky][kx] = by[kz][ky][kx] - .5*dt*sqrt(-1)*
                   (kz*ex[kz][ky][kx]-kx*ez[kz][ky][kx])
   bz[kz][ky][kx] = bz[kz][ky][kx] - .5*dt*sqrt(-1)*
                   (kx*ey[kz][ky][kx]-ky*ex[kz][ky][kx])
   the electric field is then updated a whole step using the equations:
   ex[kz][ky][kx] = ex[kz][ky][kx] + c2*dt*sqrt(-1)
                    *(ky*bz[kz][ky][kx]-kz*by[kz][ky][kx])
                    - affp*dt*cux[kz][ky][kx]*s[kz][ky][kx]
   ey[kz][ky][kx] = ey[kz][ky][kx] + c2*dt*sqrt(-1)*
                    *(kz*bx[kz][ky][kx]-kx*bz[kz][ky][kx])
                    - affp*dt*cuy[kz][ky][kx]*s[kz][ky][kx]
   ez[kz][ky][kx] = ez[kz][ky][kx] + c2*dt*sqrt(-1)
                    *(kx*by[kz][ky][kx]-ky*bx[kz][ky][kx])
                    - affp*dt*cuz[kz][ky][kx]*s[kz][ky][kx]
   the magnetic field is finally updated the remaining half step with
   the new electric field and the previous magnetic field equations.
   where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
   and s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
   j,k,l = fourier mode numbers, except for
   ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
   ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
   ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
   ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
   and similarly for bx, by, bz.
   cu[l][k][j][i] = complex current density
   exyz[l][k][j][i] = complex transverse electric field
   bxyz[l][k][j][i] = complex magnetic field
   for component i, all for fourier mode (j1,k,l)
   real(ffc[0][0][0]) = affp = normalization constant = nx*ny*nz/np,
   where np=number of particles
   aimag(ffc[l][k][j]) = finite-size particle shape factor s,
   s[kz][ky][kx] = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
   for fourier mode (j,k,l)
   ci = reciprocal of velocity of light
   dt = time interval between successive calculations
   transverse electric field energy is also calculated, using
   wf = nx*ny*nz**sum((1/affp)*|exyz[kz][ky][kx]|**2)
   magnetic field energy is also calculated, using
   wm = nx*ny*nz**sum((c2/affp)*|bxyz[kz][ky][kx]|**2)
   nx/ny/nz = system length in x/y/z direction
   nxvh = second dimension of field arrays, must be >= nxh
   nyv = third dimension of field arrays, must be >= ny
   nzv = fourth dimension of field arrays, must be >= nz
   nxhd = second dimension of form factor array, must be >= nxh
   nyhd = third dimension of form factor array, must be >= nyh
   nzhd = fourth dimension of form factor array, must be >= nzh
   requires KNC, cu, exyz, bxyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   cu, exyz, bxyz needs to have 4 components
local data                                                 */
   int nxh, nyh, nzh, nxhs, itn, j, k, l, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh;
   float dnx, dny, dnz, dth, c2, cdt, affp, anorm, dkx, dky, dkz;
   float adt, afdt;
   float at1;
   float complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9;
   double wp, ws, d0;
   __m512i v_j, v_it, v_n, v_m;
   __m512 v_dnx, v_dny, v_dnz, v_dkx, v_dky, v_dkz;
   __m512 v_zero, v_cdt, v_adt, v_afdt, v_dth, v_anorm;
   __m512 v_dk1, v_dk2, v_at1, v_at2, v_at3, v_at4;
   __m512 v_zt1, v_zt2, v_zt3, v_zt4, v_zt5, v_zt6, v_zt7;
   __m512d v_wp, v_ws, v_d;
   if (ci <= 0.0)
      return;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   dnx = 6.28318530717959/(float) nx;
   dny = 6.28318530717959/(float) ny;
   dnz = 6.28318530717959/(float) nz;
   dth = 0.5*dt;
   c2 = 1.0/(ci*ci);
   cdt = c2*dt;
   affp = creal(ffc[0]);
   adt = affp*dt;
   zero = 0.0 + 0.0*_Complex_I;
   anorm = 1.0/affp;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   v_n =  _mm512_set_epi32(15,14,11,10,9,8,13,12,7,6,3,2,1,0,5,4);
   v_m =  _mm512_set_epi32(15,14,9,8,13,12,11,10,7,6,1,0,5,4,3,2);
   v_dnx = _mm512_set1_ps(dnx);
   v_dny = _mm512_set1_ps(dny);
   v_dnz = _mm512_set1_ps(dnz);
   v_zero = _mm512_setzero_ps();
   v_cdt = _mm512_set1_ps(cdt);
   v_adt = _mm512_set1_ps(adt);
   v_dth = _mm512_set1_ps(dth);
   v_anorm = _mm512_set1_ps(anorm);
/* update electromagnetic field and sum field energies */
   ws = 0.0;
   wp = 0.0;
   v_wp = _mm512_set1_pd(0.0);
   v_ws = _mm512_set1_pd(0.0);
/* calculate the electromagnetic fields */
/* mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2 */
   for (l = 1; l < nzh; l++) {
      dkz = dnz*(float) l;
      v_dkz = _mm512_cvtfxpnt_round_adjustepi32_ps(_mm512_set1_epi32(l),
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkz = _mm512_mul_ps(v_dnz,v_dkz);
      ll = nxyhd*l;
      lj = nxvyh*l;
      l1 = nxvyh*nz - lj;
/* add kz to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dkz);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkz);
      for (k = 1; k < nyh; k++) {
         dky = dny*(float) k;
         v_it = _mm512_set1_epi32(k);
         v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dky = _mm512_mul_ps(v_dny,v_dky);
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* add ky to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),
                  v_dky);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),
                  v_dky);
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          dkx = dnx*(float) j; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                    _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
            v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
            v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                     v_dkx);
            v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                     v_dkx);
/*          afdt = adt*cimagf(ffc[j+kk+ll]); */
            v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                     _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
            v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                     _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
            v_afdt = _mm512_permute4f128_ps(v_afdt,0);
            v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                    _mm512_int2mask(13260),(__m512i)v_afdt,78);
            v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                     _mm512_int2mask(21845),(__m512i)v_afdt,177);
            v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, ky > 0, kz > 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+kj+lj)])             */
/*               + crealf(exyz[2+4*(j+kj+lj)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+kj+lj)])             */
/*               + crealf(exyz[1+4*(j+kj+lj)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+kj+lj)])               */
/*               + crealf(exyz[4*(j+kj+lj)])*_Complex_I;   */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt4 = bxyz[4*(j+kj+lj)] - dth*(dky*zt1 - dkz*zt2);   */
/*          zt5 = bxyz[1+4*(j+kj+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+kj+lj)] - dth*(dkx*zt2 - dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj+lj)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                    v_zero,v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+kj+lj)] + cdt*(dky*zt1 - dkz*zt2)   */
/*              - afdt*cu[4*(j+kj+lj)];                         */
/*          zt8 = exyz[1+4*(j+kj+lj)] + cdt*(dkz*zt3 - dkx*zt1) */
/*              - afdt*cu[1+4*(j+kj+lj)];                       */
/*          zt9 = exyz[2+4*(j+kj+lj)] + cdt*(dkx*zt2 - dky*zt3) */
/*              - afdt*cu[2+4*(j+kj+lj)];                       */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                     v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+lj)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                    v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+kj+lj)] = zt7;   */
/*          exyz[1+4*(j+kj+lj)] = zt8; */
/*          exyz[2+4*(j+kj+lj)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+kj+lj)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+kj+lj)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*             + zt9*conjf(zt9));                        */
            v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*          zt4 -= dth*(dky*zt1 - dkz*zt2); */
/*          zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2 - dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
            v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+kj+lj)] = zt4;   */
/*          bxyz[1+4*(j+kj+lj)] = zt5; */
/*          bxyz[2+4*(j+kj+lj)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+kj+lj)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+kj+lj)],v_zt5);
            }
/*          wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*             + zt6*conjf(zt6));                        */
            v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* update magnetic field half time step, ky < 0, kz > 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+k1+lj)])             */
/*               + crealf(exyz[2+4*(j+k1+lj)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+k1+lj)])             */
/*               + crealf(exyz[1+4*(j+k1+lj)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+k1+lj)])               */
/*               + crealf(exyz[4*(j+k1+lj)])*_Complex_I;   */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(12336),
                    v_zero,v_dk1);
            v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(771),
                    v_zero,v_dk2);
/*          zt4 = bxyz[4*(j+k1+lj)] + dth*(dky*zt1 + dkz*zt2);   */
/*          zt5 = bxyz[1+4*(j+k1+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+k1+lj)] - dth*(dkx*zt2 + dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1+lj)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                    v_zero,v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+k1+lj)] - cdt*(dky*zt1 + dkz*zt2)   */
/*              - afdt*cu[4*(j+k1+lj)];                         */
/*          zt8 = exyz[1+4*(j+k1+lj)] + cdt*(dkz*zt3 - dkx*zt1) */
/*              - afdt*cu[1+4*(j+k1+lj)];                       */
/*          zt9 = exyz[2+4*(j+k1+lj)] + cdt*(dkx*zt2 + dky*zt3) */
/*              - afdt*cu[2+4*(j+k1+lj)];                       */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                     v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+lj)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                    v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+k1+lj)] = zt7;   */
/*          exyz[1+4*(j+k1+lj)] = zt8; */
/*          exyz[2+4*(j+k1+lj)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+k1+lj)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+k1+lj)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*             + zt9*conjf(zt9));                        */
            v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                    v_zt6);
/*          zt4 += dth*(dky*zt1 + dkz*zt2); */
/*          zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2 + dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+k1+lj)] = zt4;   */
/*          bxyz[1+4*(j+k1+lj)] = zt5; */
/*          bxyz[2+4*(j+k1+lj)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+k1+lj)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zt5);
            }
/*          wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                     + zt6*conjf(zt6));                */
            v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                    v_zt7);
/* update magnetic field half time step, ky > 0, kz < 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+kj+l1)])             */
/*               + crealf(exyz[2+4*(j+kj+l1)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+kj+l1)])             */
/*               + crealf(exyz[1+4*(j+kj+l1)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+kj+l1)])               */
/*               + crealf(exyz[4*(j+kj+l1)])*_Complex_I;   */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(771),
                    v_zero,v_dk1);
            v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3084),
                    v_zero,v_dk2);
/*          zt4 = bxyz[4*(j+kj+l1)] - dth*(dky*zt1 + dkz*zt2);   */
/*          zt5 = bxyz[1+4*(j+kj+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+kj+l1)] - dth*(dkx*zt2 - dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj+l1)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                    v_zero,v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+kj+l1)] + cdt*(dky*zt1 + dkz*zt2)   */
/*              - afdt*cu[4*(j+kj+l1)];                         */
/*          zt8 = exyz[1+4*(j+kj+l1)] - cdt*(dkz*zt3 + dkx*zt1) */
/*              - afdt*cu[1+4*(j+kj+l1)];                       */
/*          zt9 = exyz[2+4*(j+kj+l1)] + cdt*(dkx*zt2 - dky*zt3) */
/*              - afdt*cu[2+4*(j+kj+l1)];                       */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                     v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj+l1)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                    v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+kj+l1)] = zt7;   */
/*          exyz[1+4*(j+kj+l1)] = zt8; */
/*          exyz[2+4*(j+kj+l1)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+kj+l1)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+kj+l1)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                     + zt9*conjf(zt9));                */
            v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                    v_zt6);
/*          zt4 -= dth*(dky*zt1 + dkz*zt2); */
/*          zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2 - dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+kj+l1)] = zt4;   */
/*          bxyz[1+4*(j+kj+l1)] = zt5; */
/*          bxyz[2+4*(j+kj+l1)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+kj+l1)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zt5);
            }
/*          wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                     + zt6*conjf(zt6));                */
            v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                    v_zt7);
/* update magnetic field half time step, ky < 0, kz < 0 */
/*          zt1 = -cimagf(exyz[2+4*(j+k1+l1)])             */
/*               + crealf(exyz[2+4*(j+k1+l1)])*_Complex_I; */
/*          zt2 = -cimagf(exyz[1+4*(j+k1+l1)])             */
/*               + crealf(exyz[1+4*(j+k1+l1)])*_Complex_I; */
/*          zt3 = -cimagf(exyz[4*(j+k1+l1)])              */
/*               + crealf(exyz[4*(j+k1+l1)])*_Complex_I;  */
            v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
            v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(13107),
                    v_zero,v_dk1);
            v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3855),
                    v_zero,v_dk2);
/*          zt4 = bxyz[4*(j+k1+l1)] + dth*(dky*zt1 - dkz*zt2);   */
/*          zt5 = bxyz[1+4*(j+k1+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 = bxyz[2+4*(j+k1+l1)] - dth*(dkx*zt2 + dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1+l1)]);
            v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*          zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*          zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*          zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),
                    v_zero,v_zt5);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          zt7 = exyz[4*(j+k1+l1)] - cdt*(dky*zt1 - dkz*zt2)   */
/*              - afdt*cu[4*(j+k1+l1)];                         */
/*          zt8 = exyz[1+4*(j+k1+l1)] - cdt*(dkz*zt3 + dkx*zt1) */
/*              - afdt*cu[1+4*(j+k1+l1)];                       */
/*          zt9 = exyz[2+4*(j+k1+l1)] + cdt*(dkx*zt2 + dky*zt3) */
/*              - afdt*cu[2+4*(j+k1+l1)];                       */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                     v_zt4);
            v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1+l1)]);
            v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                    v_afdt,v_zt2);
            v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*          zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*          zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*          zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
            v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),
                    v_zero,v_zt4);
            v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*          exyz[4*(j+k1+l1)] = zt7;   */
/*          exyz[1+4*(j+k1+l1)] = zt8; */
/*          exyz[2+4*(j+k1+l1)] = zt9; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&exyz[4*(j+k1+l1)],
               _mm512_int2mask(65280),v_zt4);
            }
            else {
               _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zt4);
            }
/*          ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) */
/*                     + zt9*conjf(zt9));                */
            v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                    v_zt6);
/*          zt4 += dth*(dky*zt1 - dkz*zt2); */
/*          zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*          zt6 -= dth*(dkx*zt2 + dky*zt3); */
            v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
            v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
            v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
            v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
            v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
            v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*          bxyz[4*(j+k1+l1)] = zt4;   */
/*          bxyz[1+4*(j+k1+l1)] = zt5; */
/*          bxyz[2+4*(j+k1+l1)] = zt6; */
/* zero out kx = 0 mode */
            if (j==0) {
               v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                       v_zero);
               _mm512_mask_store_ps((float *)&bxyz[4*(j+k1+l1)],
               _mm512_int2mask(65280),v_zt5);
            }
            else {
               _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zt5);
            }
/*          wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) */
/*                     + zt6*conjf(zt6));                */
            v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                    v_zt7);
/* convert to double precision before accumulating */
            v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
            v_ws = _mm512_add_pd(v_ws,v_d);
            v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
            v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
            v_wp = _mm512_add_pd(v_wp,v_d);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            dkx = dnx*(float) j;
            afdt = adt*cimagf(ffc[j+kk+ll]);
/* update magnetic field half time step, ky > 0, kz > 0 */
            zt1 = -cimagf(exyz[2+4*(j+kj+lj)])
                 + crealf(exyz[2+4*(j+kj+lj)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+kj+lj)])
                 + crealf(exyz[1+4*(j+kj+lj)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+kj+lj)])
                 + crealf(exyz[4*(j+kj+lj)])*_Complex_I;
            zt4 = bxyz[4*(j+kj+lj)] - dth*(dky*zt1 - dkz*zt2);
            zt5 = bxyz[1+4*(j+kj+lj)] - dth*(dkz*zt3 - dkx*zt1);
            zt6 = bxyz[2+4*(j+kj+lj)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+kj+lj)] + cdt*(dky*zt1 - dkz*zt2)
                - afdt*cu[4*(j+kj+lj)];
            zt8 = exyz[1+4*(j+kj+lj)] + cdt*(dkz*zt3 - dkx*zt1)
                - afdt*cu[1+4*(j+kj+lj)];
            zt9 = exyz[2+4*(j+kj+lj)] + cdt*(dkx*zt2 - dky*zt3)
                - afdt*cu[2+4*(j+kj+lj)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+kj+lj)] = zt7;
            exyz[1+4*(j+kj+lj)] = zt8;
            exyz[2+4*(j+kj+lj)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                       + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 -= dth*(dky*zt1 - dkz*zt2);
            zt5 -= dth*(dkz*zt3 - dkx*zt1);
            zt6 -= dth*(dkx*zt2 - dky*zt3);
            bxyz[4*(j+kj+lj)] = zt4;
            bxyz[1+4*(j+kj+lj)] = zt5;
            bxyz[2+4*(j+kj+lj)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                       + zt6*conjf(zt6));
            wp += (double) at1;
/* update magnetic field half time step, ky < 0, kz > 0 */
            zt1 = -cimagf(exyz[2+4*(j+k1+lj)])
                 + crealf(exyz[2+4*(j+k1+lj)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+k1+lj)])
                 + crealf(exyz[1+4*(j+k1+lj)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+k1+lj)])
                 + crealf(exyz[4*(j+k1+lj)])*_Complex_I;
            zt4 = bxyz[4*(j+k1+lj)] + dth*(dky*zt1 + dkz*zt2);
            zt5 = bxyz[1+4*(j+k1+lj)] - dth*(dkz*zt3 - dkx*zt1);
            zt6 = bxyz[2+4*(j+k1+lj)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+k1+lj)] - cdt*(dky*zt1 + dkz*zt2)
                - afdt*cu[4*(j+k1+lj)];
            zt8 = exyz[1+4*(j+k1+lj)] + cdt*(dkz*zt3 - dkx*zt1)
                - afdt*cu[1+4*(j+k1+lj)];
            zt9 = exyz[2+4*(j+k1+lj)] + cdt*(dkx*zt2 + dky*zt3)
                - afdt*cu[2+4*(j+k1+lj)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+k1+lj)] = zt7;
            exyz[1+4*(j+k1+lj)] = zt8;
            exyz[2+4*(j+k1+lj)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                       + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 += dth*(dky*zt1 + dkz*zt2);
            zt5 -= dth*(dkz*zt3 - dkx*zt1);
            zt6 -= dth*(dkx*zt2 + dky*zt3);
            bxyz[4*(j+k1+lj)] = zt4;
            bxyz[1+4*(j+k1+lj)] = zt5;
            bxyz[2+4*(j+k1+lj)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                       + zt6*conjf(zt6));
            wp += (double) at1;
/* update magnetic field half time step, ky > 0, kz < 0 */
            zt1 = -cimagf(exyz[2+4*(j+kj+l1)])
                 + crealf(exyz[2+4*(j+kj+l1)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+kj+l1)])
                 + crealf(exyz[1+4*(j+kj+l1)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+kj+l1)])
                 + crealf(exyz[4*(j+kj+l1)])*_Complex_I;
            zt4 = bxyz[4*(j+kj+l1)] - dth*(dky*zt1 + dkz*zt2);
            zt5 = bxyz[1+4*(j+kj+l1)] + dth*(dkz*zt3 + dkx*zt1);
            zt6 = bxyz[2+4*(j+kj+l1)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+kj+l1)] + cdt*(dky*zt1 + dkz*zt2)
                - afdt*cu[4*(j+kj+l1)];
            zt8 = exyz[1+4*(j+kj+l1)] - cdt*(dkz*zt3 + dkx*zt1)
                - afdt*cu[1+4*(j+kj+l1)];
            zt9 = exyz[2+4*(j+kj+l1)] + cdt*(dkx*zt2 - dky*zt3)
                - afdt*cu[2+4*(j+kj+l1)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+kj+l1)] = zt7;
            exyz[1+4*(j+kj+l1)] = zt8;
            exyz[2+4*(j+kj+l1)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                       + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 -= dth*(dky*zt1 + dkz*zt2);
            zt5 += dth*(dkz*zt3 + dkx*zt1);
            zt6 -= dth*(dkx*zt2 - dky*zt3);
            bxyz[4*(j+kj+l1)] = zt4;
            bxyz[1+4*(j+kj+l1)] = zt5;
            bxyz[2+4*(j+kj+l1)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                       + zt6*conjf(zt6));
            wp += (double) at1;
/* update magnetic field half time step, ky < 0, kz < 0 */
            zt1 = -cimagf(exyz[2+4*(j+k1+l1)])
                 + crealf(exyz[2+4*(j+k1+l1)])*_Complex_I;
            zt2 = -cimagf(exyz[1+4*(j+k1+l1)])
                 + crealf(exyz[1+4*(j+k1+l1)])*_Complex_I;
            zt3 = -cimagf(exyz[4*(j+k1+l1)])
                 + crealf(exyz[4*(j+k1+l1)])*_Complex_I;
            zt4 = bxyz[4*(j+k1+l1)] + dth*(dky*zt1 - dkz*zt2);
            zt5 = bxyz[1+4*(j+k1+l1)] + dth*(dkz*zt3 + dkx*zt1);
            zt6 = bxyz[2+4*(j+k1+l1)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
            zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
            zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
            zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
            zt7 = exyz[4*(j+k1+l1)] - cdt*(dky*zt1 - dkz*zt2)
                - afdt*cu[4*(j+k1+l1)];
            zt8 = exyz[1+4*(j+k1+l1)] - cdt*(dkz*zt3 + dkx*zt1)
                - afdt*cu[1+4*(j+k1+l1)];
            zt9 = exyz[2+4*(j+k1+l1)] + cdt*(dkx*zt2 + dky*zt3)
                - afdt*cu[2+4*(j+k1+l1)];
/* update magnetic field half time step and store electric field */
            zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
            zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
            zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
            exyz[4*(j+k1+l1)] = zt7;
            exyz[1+4*(j+k1+l1)] = zt8;
            exyz[2+4*(j+k1+l1)] = zt9;
            at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8)
                       + zt9*conjf(zt9));
            ws += (double) at1;
            zt4 += dth*(dky*zt1 - dkz*zt2);
            zt5 += dth*(dkz*zt3 + dkx*zt1);
            zt6 -= dth*(dkx*zt2 + dky*zt3);
            bxyz[4*(j+k1+l1)] = zt4;
            bxyz[1+4*(j+k1+l1)] = zt5;
            bxyz[2+4*(j+k1+l1)] = zt6;
            at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5)
                       + zt6*conjf(zt6));
            wp += (double) at1;
         }
      }
/* mode numbers kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         dky = dny*(float) k;
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
         afdt = adt*cimagf(ffc[kk+ll]);
/* update magnetic field half time step, kz > 0 */
         zt1 = -cimagf(exyz[2+4*(kj+lj)])
              + crealf(exyz[2+4*(kj+lj)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(kj+lj)])
              + crealf(exyz[1+4*(kj+lj)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(kj+lj)])
              + crealf(exyz[4*(kj+lj)])*_Complex_I;
         zt4 = bxyz[4*(kj+lj)] - dth*(dky*zt1 - dkz*zt2);
         zt5 = bxyz[1+4*(kj+lj)] - dth*(dkz*zt3);
         zt6 = bxyz[2+4*(kj+lj)] + dth*(dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(kj+lj)] + cdt*(dky*zt1 - dkz*zt2)
             - afdt*cu[4*(kj+lj)];
         zt8 = exyz[1+4*(kj+lj)] + cdt*(dkz*zt3) - afdt*cu[1+4*(kj+lj)];
         zt9 = exyz[2+4*(kj+lj)] - cdt*(dky*zt3) - afdt*cu[2+4*(kj+lj)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(kj+lj)] = zt7;
         exyz[1+4*(kj+lj)] = zt8;
         exyz[2+4*(kj+lj)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dky*zt1 - dkz*zt2);
         zt5 -= dth*(dkz*zt3);
         zt6 += dth*(dky*zt3);
         bxyz[4*(kj+lj)] = zt4;
         bxyz[1+4*(kj+lj)] = zt5;
         bxyz[2+4*(kj+lj)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(k1+lj)] = zero;
         bxyz[1+4*(k1+lj)] = zero;
         bxyz[2+4*(k1+lj)] = zero;
         exyz[4*(k1+lj)] = zero;
         exyz[1+4*(k1+lj)] = zero;
         exyz[2+4*(k1+lj)] = zero;
/* update magnetic field half time step, kz < 0 */
         zt1 = -cimagf(exyz[2+4*(kj+l1)])
              + crealf(exyz[2+4*(kj+l1)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(kj+l1)])
              + crealf(exyz[1+4*(kj+l1)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(kj+l1)])
              + crealf(exyz[4*(kj+l1)])*_Complex_I;
         zt4 = bxyz[4*(kj+l1)] - dth*(dky*zt1 + dkz*zt2);
         zt5 = bxyz[1+4*(kj+l1)] + dth*(dkz*zt3);
         zt6 = bxyz[2+4*(kj+l1)] + dth*(dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(kj+l1)] + cdt*(dky*zt1 + dkz*zt2)
             - afdt*cu[4*(kj+l1)];
         zt8 = exyz[1+4*(kj+l1)] - cdt*(dkz*zt3) - afdt*cu[1+4*(kj+l1)];
         zt9 = exyz[2+4*(kj+l1)] - cdt*(dky*zt3) - afdt*cu[2+4*(kj+l1)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(kj+l1)] = zt7;
         exyz[1+4*(kj+l1)] = zt8;
         exyz[2+4*(kj+l1)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dky*zt1 + dkz*zt2);
         zt5 += dth*(dkz*zt3);
         zt6 += dth*(dky*zt3);
         bxyz[4*(kj+l1)] = zt4;
         bxyz[1+4*(kj+l1)] = zt5;
         bxyz[2+4*(kj+l1)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(k1+l1)] = zero;
         bxyz[1+4*(k1+l1)] = zero;
         bxyz[2+4*(k1+l1)] = zero;
         exyz[4*(k1+l1)] = zero;
         exyz[1+4*(k1+l1)] = zero;
         exyz[2+4*(k1+l1)] = zero;
      }
/* mode numbers ky = 0, ny/2 */
      k1 = nxvh*nyh;
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(12336),v_zero);
      v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(771),v_zero);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/*       afdt = adt*cimagf(ffc[j+ll]); */
         v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                  _mm512_int2mask(15),(float *)&ffc[j+ll]);
         v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                  _mm512_int2mask(15),(float *)&ffc[j+ll+8]);
         v_afdt = _mm512_permute4f128_ps(v_afdt,0);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(13260),(__m512i)v_afdt,78);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(21845),(__m512i)v_afdt,177);
         v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, kz > 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+lj)])             */
/*            + crealf(exyz[2+4*(j+lj)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+lj)])             */
/*            + crealf(exyz[1+4*(j+lj)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+lj)])               */
/*            + crealf(exyz[4*(j+lj)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt4 = bxyz[4*(j+lj)] + dth*(dkz*zt2);             */
/*       zt5 = bxyz[1+4*(j+lj)] - dth*(dkz*zt3 - dkx*zt1); */
/*       zt6 = bxyz[2+4*(j+lj)] - dth*(dkx*zt2);           */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+lj)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+lj)] - cdt*(dkz*zt2) - afdt*cu[4*(j+lj)];     */
/*       zt8 = exyz[1+4*(j+lj)] + cdt*(dkz*zt3 - dkx*zt1)              */
/*           - afdt*cu[1+4*(j+lj)];                                    */
/*       zt9 = exyz[2+4*(j+lj)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+lj)]; */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                 v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+lj)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                 v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+lj)] = zt7;   */
/*       exyz[1+4*(j+lj)] = zt8; */
/*       exyz[2+4*(j+lj)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                     v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+lj)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+lj)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*       zt4 += dth*(dkz*zt2);           */
/*       zt5 -= dth*(dkz*zt3 - dkx*zt1); */
/*       zt6 -= dth*(dkx*zt2);           */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+lj)] = zt4;   */
/*       bxyz[1+4*(j+lj)] = zt5; */
/*       bxyz[2+4*(j+lj)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+lj)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+lj)],v_zt5);
         }
/*       wp +=  anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/*       bxyz[4*(j+k1+lj)] = zero;   */
/*       bxyz[1+4*(j+k1+lj)] = zero; */
/*       bxyz[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+lj)],v_zero);
/*       exyz[4*(j+k1+lj)] = zero;   */
/*       exyz[1+4*(j+k1+lj)] = zero; */
/*       exyz[2+4*(j+k1+lj)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+k1+lj)],v_zero);
/* update magnetic field half time step, kz > 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+l1)])             */
/*            + crealf(exyz[2+4*(j+l1)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+l1)])             */
/*            + crealf(exyz[1+4*(j+l1)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+l1)])               */
/*            + crealf(exyz[4*(j+l1)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                  v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(771),v_zero,
                  v_dk1);
         v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(3084),v_zero,
                 v_dk2);
/*       zt4 = bxyz[4*(j+l1)] - dth*(dkz*zt2);             */
/*       zt5 = bxyz[1+4*(j+l1)] + dth*(dkz*zt3 + dkx*zt1); */
/*       zt6 = bxyz[2+4*(j+l1)] - dth*(dkx*zt2);           */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+l1)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+l1)] + cdt*(dkz*zt2) - afdt*cu[4*(j+l1)];     */
/*       zt8 = exyz[1+4*(j+l1)] - cdt*(dkz*zt3 + dkx*zt1)              */
/*           - afdt*cu[1+4*(j+l1)];                                    */
/*       zt9 = exyz[2+4*(j+l1)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+l1)]; */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                 v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+l1)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                  v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+l1)] = zt7;   */
/*       exyz[1+4*(j+l1)] = zt8; */
/*       exyz[2+4*(j+l1)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
             v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                     v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+l1)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+l1)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                 v_zt6);
/*       zt4 -= dth*(dkz*zt2);           */
/*       zt5 += dth*(dkz*zt3 + dkx*zt1); */
/*       zt6 -= dth*(dkx*zt2);           */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+l1)] = zt4;   */
/*       bxyz[1+4*(j+l1)] = zt5; */
/*       bxyz[2+4*(j+l1)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+l1)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zt5);
         }
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                  v_zt7);
/* convert to double precision before accumulating */
         v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
         v_ws = _mm512_add_pd(v_ws,v_d);
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+k1+l1)] = zero;   */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*       exyz[4*(j+k1+l1)] = zero;   */
/*       exyz[1+4*(j+k1+l1)] = zero; */
/*       exyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;  
         afdt = adt*cimagf(ffc[j+ll]);
/* update magnetic field half time step, kz > 0 */
         zt1 = -cimagf(exyz[2+4*(j+lj)])
              + crealf(exyz[2+4*(j+lj)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+lj)])
              + crealf(exyz[1+4*(j+lj)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+lj)])
              + crealf(exyz[4*(j+lj)])*_Complex_I;
         zt4 = bxyz[4*(j+lj)] + dth*(dkz*zt2);
         zt5 = bxyz[1+4*(j+lj)] - dth*(dkz*zt3 - dkx*zt1);
         zt6 = bxyz[2+4*(j+lj)] - dth*(dkx*zt2);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+lj)] - cdt*(dkz*zt2) - afdt*cu[4*(j+lj)];
         zt8 = exyz[1+4*(j+lj)] + cdt*(dkz*zt3 - dkx*zt1)
             - afdt*cu[1+4*(j+lj)];
         zt9 = exyz[2+4*(j+lj)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+lj)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+lj)] = zt7;
         exyz[1+4*(j+lj)] = zt8;
         exyz[2+4*(j+lj)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 += dth*(dkz*zt2);
         zt5 -= dth*(dkz*zt3 - dkx*zt1);
         zt6 -= dth*(dkx*zt2);
         bxyz[4*(j+lj)] = zt4;
         bxyz[1+4*(j+lj)] = zt5;
         bxyz[2+4*(j+lj)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(j+k1+lj)] = zero;
         bxyz[1+4*(j+k1+lj)] = zero;
         bxyz[2+4*(j+k1+lj)] = zero;
         exyz[4*(j+k1+lj)] = zero;
         exyz[1+4*(j+k1+lj)] = zero;
         exyz[2+4*(j+k1+lj)] = zero;
/* update magnetic field half time step, kz > 0 */
         zt1 = -cimagf(exyz[2+4*(j+l1)])
              + crealf(exyz[2+4*(j+l1)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+l1)])
              + crealf(exyz[1+4*(j+l1)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+l1)])
              + crealf(exyz[4*(j+l1)])*_Complex_I;
         zt4 = bxyz[4*(j+l1)] - dth*(dkz*zt2);
         zt5 = bxyz[1+4*(j+l1)] + dth*(dkz*zt3 + dkx*zt1);
         zt6 = bxyz[2+4*(j+l1)] - dth*(dkx*zt2);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+l1)] + cdt*(dkz*zt2) - afdt*cu[4*(j+l1)];
         zt8 = exyz[1+4*(j+l1)] - cdt*(dkz*zt3 + dkx*zt1)
             - afdt*cu[1+4*(j+l1)];
         zt9 = exyz[2+4*(j+l1)] + cdt*(dkx*zt2) - afdt*cu[2+4*(j+l1)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+l1)] = zt7;
         exyz[1+4*(j+l1)] = zt8;
         exyz[2+4*(j+l1)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dkz*zt2);
         zt5 += dth*(dkz*zt3 + dkx*zt1);
         zt6 -= dth*(dkx*zt2);
         bxyz[4*(j+l1)] = zt4;
         bxyz[1+4*(j+l1)] = zt5;
         bxyz[2+4*(j+l1)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         exyz[4*(j+k1+l1)] = zero;
         exyz[1+4*(j+k1+l1)] = zero;
         exyz[2+4*(j+k1+l1)] = zero;
      }
/* mode numbers kx = 0, nx/2 */
      afdt = adt*cimagf(ffc[ll]);
/* update magnetic field half time step */
      zt2 = -cimagf(exyz[1+4*(lj)]) + crealf(exyz[1+4*(lj)])*_Complex_I;
      zt3 = -cimagf(exyz[4*(lj)]) + crealf(exyz[4*(lj)])*_Complex_I;
      zt4 = bxyz[4*lj] + dth*(dkz*zt2);
      zt5 = bxyz[1+4*lj] - dth*(dkz*zt3);
/* update electric field whole time step */
      zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exyz[4*lj] - cdt*(dkz*zt2) - afdt*cu[4*lj];
      zt8 = exyz[1+4*lj] + cdt*(dkz*zt3) - afdt*cu[1+4*lj];
/* update magnetic field half time step and store electric field */
      zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exyz[4*lj] = zt7;
      exyz[1+4*lj] = zt8;
      exyz[2+4*lj] = zero;
      at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8));
      ws += (double) at1;
      zt4 += dth*(dkz*zt2);
      zt5 -= dth*(dkz*zt3);
      bxyz[4*lj] = zt4;
      bxyz[1+4*lj] = zt5;
      bxyz[2+4*lj] = zero;
      at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5));
      wp += (double) at1;
      bxyz[4*(k1+lj)] = zero;
      bxyz[1+4*(k1+lj)] = zero;
      bxyz[2+4*(k1+lj)] = zero;
      exyz[4*(k1+lj)] = zero;
      exyz[1+4*(k1+lj)] = zero;
      exyz[2+4*(k1+lj)] = zero;
      bxyz[4*l1] = zero;
      bxyz[1+4*l1] = zero;
      bxyz[2+4*l1] = zero;
      exyz[4*l1] = zero;
      exyz[1+4*l1] = zero;
      exyz[2+4*l1] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      exyz[4*(k1+l1)] = zero;
      exyz[1+4*(k1+l1)] = zero;
      exyz[2+4*(k1+l1)]= zero;
   }
/* mode numbers kz = 0, nz/2 */
   l1 = nxvyh*nzh;
   for (k = 1; k < nyh; k++) {
/*    dky = dny*(float) k; */
      v_it = _mm512_set1_epi32(k);
      v_dky = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dky = _mm512_mul_ps(v_dny,v_dky);
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
/* add ky to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dky);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(771),v_dky);
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       dkx = dnx*(float) j; */
         v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
         v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
                 _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
         v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
         v_dk1  = _mm512_mask_mov_ps(v_dk1,_mm512_int2mask(3084),
                  v_dkx);
         v_dk2  = _mm512_mask_mov_ps(v_dk2,_mm512_int2mask(12336),
                  v_dkx);
/*       afdt = adt*cimagf(ffc[j+kk]); */
         v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                  _mm512_int2mask(15),(float *)&ffc[j+kk]);
         v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
                  _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
         v_afdt = _mm512_permute4f128_ps(v_afdt,0);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(13260),(__m512i)v_afdt,78);
         v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
                  _mm512_int2mask(21845),(__m512i)v_afdt,177);
         v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step, ky > 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+kj)])             */
/*            + crealf(exyz[2+4*(j+kj)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+kj)])             */
/*            + crealf(exyz[1+4*(j+kj)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+kj)])               */
/*            + crealf(exyz[4*(j+kj)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt4 = bxyz[4*(j+kj)] - dth*(dky*zt1);             */
/*       zt5 = bxyz[1+4*(j+kj)] + dth*(dkx*zt1);           */
/*       zt6 = bxyz[2+4*(j+kj)] - dth*(dkx*zt2 - dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+kj)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+kj)] + cdt*(dky*zt1) - afdt*cu[4*(j+kj)];     */
/*       zt8 = exyz[1+4*(j+kj)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+kj)]; */
/*       zt9 = exyz[2+4*(j+kj)] + cdt*(dkx*zt2 - dky*zt3)              */
/*           - afdt*cu[2+4*(j+kj)];                                    */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),
                 v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+kj)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                 v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+kj)] = zt7;   */
/*       exyz[1+4*(j+kj)] = zt8; */
/*       exyz[2+4*(j+kj)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+kj)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+kj)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*       zt4 -= dth*(dky*zt1);           */
/*       zt5 += dth*(dkx*zt1);           */
/*       zt6 -= dth*(dkx*zt2 - dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
         v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+kj)] = zt4;   */
/*       bxyz[1+4*(j+kj)] = zt5; */
/*       bxyz[2+4*(j+kj)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+kj)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+kj)],v_zt5);
         }
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* update magnetic field half time step, ky < 0 */
/*       zt1 = -cimagf(exyz[2+4*(j+k1)])             */
/*            + crealf(exyz[2+4*(j+k1)])*_Complex_I; */
/*       zt2 = -cimagf(exyz[1+4*(j+k1)])             */
/*            + crealf(exyz[1+4*(j+k1)])*_Complex_I; */
/*       zt3 = -cimagf(exyz[4*(j+k1)])               */
/*            + crealf(exyz[4*(j+k1)])*_Complex_I;   */
         v_zt4 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
         v_at2 = _mm512_mask_sub_ps(v_dk1,_mm512_int2mask(12336),v_zero,
                  v_dk1);
         v_at3 = _mm512_mask_sub_ps(v_dk2,_mm512_int2mask(771),v_zero,
                 v_dk2);
/*       zt4 = bxyz[4*(j+k1)] + dth*(dky*zt1);             */
/*       zt5 = bxyz[1+4*(j+k1)] + dth*(dkx*zt1);           */
/*       zt6 = bxyz[2+4*(j+k1)] - dth*(dkx*zt2 + dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt2 = _mm512_load_ps((float *)&bxyz[4*(j+k1)]);
         v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*       zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*       zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
/*       zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
                 v_zt5);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       zt7 = exyz[4*(j+k1)] - cdt*(dky*zt1) - afdt*cu[4*(j+k1)];     */
/*       zt8 = exyz[1+4*(j+k1)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+k1)]; */
/*       zt9 = exyz[2+4*(j+k1)] + cdt*(dkx*zt2 + dky*zt3)              */
/*           - afdt*cu[2+4*(j+k1)];                                    */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),v_zt4);
         v_zt2 = _mm512_load_ps((float *)&cu[4*(j+k1)]);
         v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),
                 v_afdt,v_zt2);
         v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*       zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*       zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
/*       zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I; */
         v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
                 v_zt4);
         v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*       exyz[4*(j+k1)] = zt7;   */
/*       exyz[1+4*(j+k1)] = zt8; */
/*       exyz[2+4*(j+k1)] = zt9; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&exyz[4*(j+k1)],
            _mm512_int2mask(65280),v_zt4);
         }
         else {
            _mm512_store_ps((float *)&exyz[4*(j+k1)],v_zt4);
         }
/*       ws += anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9)); */
         v_zt6 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4),
                 v_zt6);
/*       zt4 += dth*(dky*zt1);           */
/*       zt5 += dth*(dkx*zt1);           */
/*       zt6 -= dth*(dkx*zt2 + dky*zt3); */
         v_zt1 = _mm512_mul_ps(v_at2,v_zt3);
         v_zt2 = _mm512_mul_ps(v_at3,v_zt3);
         v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
         v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
         v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
         v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*       bxyz[4*(j+k1)] = zt4;   */
/*       bxyz[1+4*(j+k1)] = zt5; */
/*       bxyz[2+4*(j+k1)] = zt6; */
/* zero out kx = 0 mode */
         if (j==0) {
            v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),
                    v_zero);
            _mm512_mask_store_ps((float *)&bxyz[4*(j+k1)],
            _mm512_int2mask(65280),v_zt5);
         }
         else {
            _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zt5);
         }
/*       wp += anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6)); */
         v_zt7 = _mm512_fmadd_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5),
                 v_zt7);
/* convert to double precision before accumulating */
         v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
         v_ws = _mm512_add_pd(v_ws,v_d);
         v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
         v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
         v_wp = _mm512_add_pd(v_wp,v_d);
/*       bxyz[4*(j+kj+l1)] = zero; */
/*       bxyz[1+4*(j+kj+l1)] = zero; */
/*       bxyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+kj+l1)],v_zero);
/*       exyz[4*(j+kj+l1)] = zero; */
/*       exyz[1+4*(j+kj+l1)] = zero; */
/*       exyz[2+4*(j+kj+l1)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+kj+l1)],v_zero);
/*       bxyz[4*(j+k1+l1)] = zero; */
/*       bxyz[1+4*(j+k1+l1)] = zero; */
/*       bxyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*       exyz[4*(j+k1+l1)] = zero; */
/*       exyz[1+4*(j+k1+l1)] = zero; */
/*       exyz[2+4*(j+k1+l1)] = zero; */
         _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         dkx = dnx*(float) j;
         afdt = adt*cimagf(ffc[j+kk]);
/* update magnetic field half time step, ky > 0 */
         zt1 = -cimagf(exyz[2+4*(j+kj)])
              + crealf(exyz[2+4*(j+kj)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+kj)])
              + crealf(exyz[1+4*(j+kj)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+kj)])
              + crealf(exyz[4*(j+kj)])*_Complex_I;
         zt4 = bxyz[4*(j+kj)] - dth*(dky*zt1);
         zt5 = bxyz[1+4*(j+kj)] + dth*(dkx*zt1);
         zt6 = bxyz[2+4*(j+kj)] - dth*(dkx*zt2 - dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+kj)] + cdt*(dky*zt1) - afdt*cu[4*(j+kj)];
         zt8 = exyz[1+4*(j+kj)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+kj)];
         zt9 = exyz[2+4*(j+kj)] + cdt*(dkx*zt2 - dky*zt3) 
             - afdt*cu[2+4*(j+kj)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+kj)] = zt7;
         exyz[1+4*(j+kj)] = zt8;
         exyz[2+4*(j+kj)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 -= dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 - dky*zt3);
         bxyz[4*(j+kj)] = zt4;
         bxyz[1+4*(j+kj)] = zt5;
         bxyz[2+4*(j+kj)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
/* update magnetic field half time step, ky < 0 */
         zt1 = -cimagf(exyz[2+4*(j+k1)])
              + crealf(exyz[2+4*(j+k1)])*_Complex_I;
         zt2 = -cimagf(exyz[1+4*(j+k1)])
              + crealf(exyz[1+4*(j+k1)])*_Complex_I;
         zt3 = -cimagf(exyz[4*(j+k1)])
              + crealf(exyz[4*(j+k1)])*_Complex_I;
         zt4 = bxyz[4*(j+k1)] + dth*(dky*zt1);
         zt5 = bxyz[1+4*(j+k1)] + dth*(dkx*zt1);
         zt6 = bxyz[2+4*(j+k1)] - dth*(dkx*zt2 + dky*zt3);
/* update electric field whole time step */
         zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
         zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
         zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
         zt7 = exyz[4*(j+k1)] - cdt*(dky*zt1) - afdt*cu[4*(j+k1)];
         zt8 = exyz[1+4*(j+k1)] - cdt*(dkx*zt1) - afdt*cu[1+4*(j+k1)];
         zt9 = exyz[2+4*(j+k1)] + cdt*(dkx*zt2 + dky*zt3)
             - afdt*cu[2+4*(j+k1)];
/* update magnetic field half time step and store electric field */
         zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
         zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
         zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
         exyz[4*(j+k1)] = zt7;
         exyz[1+4*(j+k1)] = zt8;
         exyz[2+4*(j+k1)] = zt9;
         at1 = anorm*(zt7*conjf(zt7) + zt8*conjf(zt8) + zt9*conjf(zt9));
         ws += (double) at1;
         zt4 += dth*(dky*zt1);
         zt5 += dth*(dkx*zt1);
         zt6 -= dth*(dkx*zt2 + dky*zt3);
         bxyz[4*(j+k1)] = zt4;
         bxyz[1+4*(j+k1)] = zt5;
         bxyz[2+4*(j+k1)] = zt6;
         at1 = anorm*(zt4*conjf(zt4) + zt5*conjf(zt5) + zt6*conjf(zt6));
         wp += (double) at1;
         bxyz[4*(j+kj+l1)] = zero;
         bxyz[1+4*(j+kj+l1)] = zero;
         bxyz[2+4*(j+kj+l1)] = zero;
         exyz[4*(j+kj+l1)] = zero;
         exyz[1+4*(j+kj+l1)] = zero;
         exyz[2+4*(j+kj+l1)] = zero;
         bxyz[4*(j+k1+l1)] = zero;
         bxyz[1+4*(j+k1+l1)] = zero;
         bxyz[2+4*(j+k1+l1)] = zero;
         exyz[4*(j+k1+l1)] = zero;
         exyz[1+4*(j+k1+l1)] = zero;
         exyz[2+4*(j+k1+l1)] = zero;
      }
   }
/* mode numbers kx = 0, nx/2 */
   for (k = 1; k < nyh; k++) {
      dky = dny*(float) k;
      kk = nxhd*k;
      kj = nxvh*k;
      k1 = nxvh*ny - kj;
      afdt = adt*cimagf(ffc[kk]);
/* update magnetic field half time step */
      zt1 = -cimagf(exyz[2+4*(kj)]) + crealf(exyz[2+4*(kj)])*_Complex_I;
      zt3 = -cimagf(exyz[4*(kj)]) + crealf(exyz[4*(kj)])*_Complex_I;
      zt4 = bxyz[4*kj] - dth*(dky*zt1);
      zt6 = bxyz[2+4*kj] + dth*(dky*zt3);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt3 = -cimagf(zt4) + crealf(zt4)*_Complex_I;
      zt7 = exyz[4*kj] + cdt*(dky*zt1) - afdt*cu[4*kj];
      zt9 = exyz[2+4*kj] - cdt*(dky*zt3) - afdt*cu[2+4*kj];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt3 = -cimagf(zt7) + crealf(zt7)*_Complex_I;
      exyz[4*kj] = zt7;
      exyz[1+4*kj] = zero;
      exyz[2+4*kj] = zt9;
      at1 = anorm*(zt7*conjf(zt7) + zt9*conjf(zt9));
      ws += (double) at1;
      zt4 -= dth*(dky*zt1);
      zt6 += dth*(dky*zt3);
      bxyz[4*kj] = zt4;
      bxyz[1+4*kj] = zero;
      bxyz[2+4*kj] = zt6;
      at1 = anorm*(zt4*conjf(zt4) + zt6*conjf(zt6));
      wp += (double) at1;
      bxyz[4*k1] = zero;
      bxyz[1+4*k1] = zero;
      bxyz[2+4*k1] = zero;
      exyz[4*k1] = zero;
      exyz[1+4*k1] = zero;
      exyz[2+4*k1] = zero;
      bxyz[4*(kj+l1)] = zero;
      bxyz[1+4*(kj+l1)] = zero;
      bxyz[2+4*(kj+l1)]= zero;
      exyz[4*(kj+l1)] = zero;
      exyz[1+4*(kj+l1)] = zero;
      exyz[2+4*(kj+l1)] = zero;
      bxyz[4*(k1+l1)] = zero;
      bxyz[1+4*(k1+l1)] = zero;
      bxyz[2+4*(k1+l1)] = zero;
      exyz[4*(k1+l1)] = zero;
      exyz[1+4*(k1+l1)] = zero;
      exyz[2+4*(k1+l1)] = zero;
   }
/* mode numbers ky = 0, ny/2 */
   k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
   for (j = 0; j < nxhs; j+=2) {
/*    dkx = dnx*(float) j; */
      v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
      v_dkx = _mm512_cvtfxpnt_round_adjustepi32_ps(v_it,
              _MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
      v_dkx = _mm512_mul_ps(v_dnx,v_dkx);
/* add kx to curl operators */
      v_dk1  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(3084),v_dkx);
      v_dk2  = _mm512_mask_mov_ps(v_zero,_mm512_int2mask(12336),v_dkx);
/*    afdt = adt*cimagf(ffc[j]); */
      v_afdt = _mm512_mask_loadunpacklo_ps(v_zero,
                _mm512_int2mask(15),(float *)&ffc[j]);
      v_afdt = _mm512_mask_loadunpackhi_ps(v_afdt,
               _mm512_int2mask(15),(float *)&ffc[j+8]);
      v_afdt = _mm512_permute4f128_ps(v_afdt,0);
      v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
               _mm512_int2mask(13260),(__m512i)v_afdt,78);
      v_afdt = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_afdt,
               _mm512_int2mask(21845),(__m512i)v_afdt,177);
      v_afdt = _mm512_mul_ps(v_adt,v_afdt);
/* update magnetic field half time step */
/*    zt1 = -cimagf(exyz[2+4*j]) + crealf(exyz[2+4*j])*_Complex_I; */
/*    zt2 = -cimagf(exyz[1+4*j]) + crealf(exyz[1+4*j])*_Complex_I; */
      v_zt4 = _mm512_load_ps((float *)&exyz[4*j]);
      v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
               v_zt4);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    zt5 = bxyz[1+4*j] + dth*(dkx*zt1); */
/*    zt6 = bxyz[2+4*j] - dth*(dkx*zt2); */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
      v_zt2 = _mm512_load_ps((float *)&bxyz[4*j]);
      v_zt5 = _mm512_sub_ps(v_zt2,v_zt1);
/* update electric field whole time step */
/*    zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I; */
/*    zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I; */
      v_zt3 = _mm512_mask_sub_ps(v_zt5,_mm512_int2mask(43690),v_zero,
              v_zt5);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    zt8 = exyz[1+4*j] - cdt*(dkx*zt1) - afdt*cu[1+4*j]; */
/*    zt9 = exyz[2+4*j] + cdt*(dkx*zt2) - afdt*cu[2+4*j]; */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_fmadd_ps(v_cdt,_mm512_sub_ps(v_zt1,v_zt2),v_zt4);
      v_zt2 = _mm512_load_ps((float *)&cu[4*j]);
      v_zt2 = _mm512_mask_mul_ps(v_zero,_mm512_int2mask(16191),v_afdt,
              v_zt2);
      v_zt4 = _mm512_sub_ps(v_zt1,v_zt2);
/* update magnetic field half time step and store electric field */
/*    zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I; */
/*    zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I; */
      v_zt3 = _mm512_mask_sub_ps(v_zt4,_mm512_int2mask(43690),v_zero,
              v_zt4);
      v_zt3 = (__m512)_mm512_shuffle_epi32((__m512i)v_zt3,177);
/*    exyz[4*j] = zero;  */
/*    exyz[1+4*j] = zt8; */
/*    exyz[2+4*j] = zt9; */
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt4 = _mm512_mask_mov_ps(v_zt4,_mm512_int2mask(255),v_zero);
         _mm512_mask_store_ps((float *)&exyz[4*j],
         _mm512_int2mask(65280),v_zt4);
      }
      else {
         _mm512_store_ps((float *)&exyz[4*j],v_zt4);
      }
/*    ws += anorm*(zt8*conjf(zt8) + zt9*conjf(zt9)); */
      v_zt6 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt4,v_zt4));
/*    zt5 += dth*(dkx*zt1); */
/*    zt6 -= dth*(dkx*zt2); */
      v_zt1 = _mm512_mul_ps(v_dk1,v_zt3);
      v_zt2 = _mm512_mul_ps(v_dk2,v_zt3);
      v_zt1 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_zt1);
      v_zt2 = (__m512)_mm512_permutevar_epi32(v_m,(__m512i)v_zt2);
      v_zt1 = _mm512_mul_ps(v_dth,_mm512_sub_ps(v_zt1,v_zt2));
      v_zt5 = _mm512_sub_ps(v_zt5,v_zt1);
/*    bxyz[4*j] = zero;  */
/*    bxyz[1+4*j] = zt5; */
/*    bxyz[2+4*j] = zt6; */
/* zero out kx = 0 mode */
      if (j==0) {
         v_zt5 = _mm512_mask_mov_ps(v_zt5,_mm512_int2mask(255),v_zero);
         _mm512_mask_store_ps((float *)&bxyz[4*j],
         _mm512_int2mask(65280),v_zt5);
      }
      else {
         _mm512_store_ps((float *)&bxyz[4*j],v_zt5);
      }
/*    wp += anorm*(zt5*conjf(zt5) + zt6*conjf(zt6)); */
      v_zt7 = _mm512_mul_ps(v_anorm,_mm512_mul_ps(v_zt5,v_zt5));
/* convert to double precision before accumulating */
      v_ws = _mm512_add_pd(v_ws,_mm512_cvtpslo_pd(v_zt6));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt6,78));
      v_ws = _mm512_add_pd(v_ws,v_d);
      v_wp = _mm512_add_pd(v_wp,_mm512_cvtpslo_pd(v_zt7));
      v_d = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(v_zt7,78));
      v_wp = _mm512_add_pd(v_wp,v_d);
/*    bxyz[4*(j+k1)] = zero;   */
/*    bxyz[1+4*(j+k1)] = zero; */
/*    bxyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1)],v_zero);
/*    exyz[4*(j+k1)] = zero;   */
/*    exyz[1+4*(j+k1)] = zero; */
/*    exyz[2+4*(j+k1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+k1)],v_zero);
/*    bxyz[4*(j+l1)] = zero;   */
/*    bxyz[1+4*(j+l1)] = zero; */
/*    bxyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+l1)],v_zero);
/*    exyz[4*(j+l1)] = zero;   */
/*    exyz[1+4*(j+l1)] = zero; */
/*    exyz[2+4*(j+l1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+l1)],v_zero);
/*    bxyz[4*(j+k1+l1)] = zero;   */
/*    bxyz[1+4*(j+k1+l1)] = zero; */
/*    bxyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&bxyz[4*(j+k1+l1)],v_zero);
/*    exyz[4*(j+k1+l1)] = zero;   */
/*    exyz[1+4*(j+k1+l1)] = zero; */
/*    exyz[2+4*(j+k1+l1)] = zero; */
      _mm512_store_ps((float *)&exyz[4*(j+k1+l1)],v_zero);
   }
/* loop over remaining elements */
   for (j = itn; j < nxh; j++) {
      dkx = dnx*(float) j;
      afdt = adt*cimagf(ffc[j]);
/* update magnetic field half time step */
      zt1 = -cimagf(exyz[2+4*j]) + crealf(exyz[2+4*j])*_Complex_I;
      zt2 = -cimagf(exyz[1+4*j]) + crealf(exyz[1+4*j])*_Complex_I;
      zt5 = bxyz[1+4*j] + dth*(dkx*zt1);
      zt6 = bxyz[2+4*j] - dth*(dkx*zt2);
/* update electric field whole time step */
      zt1 = -cimagf(zt6) + crealf(zt6)*_Complex_I;
      zt2 = -cimagf(zt5) + crealf(zt5)*_Complex_I;
      zt8 = exyz[1+4*j] - cdt*(dkx*zt1) - afdt*cu[1+4*j];
      zt9 = exyz[2+4*j] + cdt*(dkx*zt2) - afdt*cu[2+4*j];
/* update magnetic field half time step and store electric field */
      zt1 = -cimagf(zt9) + crealf(zt9)*_Complex_I;
      zt2 = -cimagf(zt8) + crealf(zt8)*_Complex_I;
      exyz[4*j] = zero;
      exyz[1+4*j] = zt8;
      exyz[2+4*j] = zt9;
      at1 = anorm*(zt8*conjf(zt8) + zt9*conjf(zt9));
      ws += (double) at1;
      zt5 += dth*(dkx*zt1);
      zt6 -= dth*(dkx*zt2);
      bxyz[4*j] = zero;
      bxyz[1+4*j] = zt5;
      bxyz[2+4*j] = zt6;
      at1 = anorm*(zt5*conjf(zt5) + zt6*conjf(zt6));
      wp += (double) at1;
      bxyz[4*(j+k1)] = zero;
      bxyz[1+4*(j+k1)] = zero;
      bxyz[2+4*(j+k1)] = zero;
      exyz[4*(j+k1)] = zero;
      exyz[1+4*(j+k1)] = zero;
      exyz[2+4*(j+k1)] = zero;
      bxyz[4*(j+l1)] = zero;
      bxyz[1+4*(j+l1)] = zero;
      bxyz[2+4*(j+l1)] = zero;
      exyz[4*(j+l1)] = zero;
      exyz[1+4*(j+l1)] = zero;
      exyz[2+4*(j+l1)] = zero;
      bxyz[4*(j+k1+l1)] = zero;
      bxyz[1+4*(j+k1+l1)] = zero;
      bxyz[2+4*(j+k1+l1)] = zero;
      exyz[4*(j+k1+l1)] = zero;
      exyz[1+4*(j+k1+l1)] = zero;
      exyz[2+4*(j+k1+l1)] = zero;
   }
   bxyz[0] = zero;
   bxyz[1] = zero;
   bxyz[2] = zero;
   exyz[0] = zero;
   exyz[1] = zero;
   exyz[2]= zero;
   bxyz[4*k1] = zero;
   bxyz[1+4*k1] = zero;
   bxyz[2+4*k1] = zero;
   exyz[4*k1] = zero;
   exyz[1+4*k1] = zero;
   exyz[2+4*k1] = zero;
   bxyz[4*l1] = zero;
   bxyz[1+4*l1] = zero;
   bxyz[2+4*l1] = zero;
   exyz[4*l1] = zero;
   exyz[1+4*l1] = zero;
   exyz[2+4*l1] = zero;
   bxyz[4*(k1+l1)] = zero;
   bxyz[1+4*(k1+l1)] = zero;
   bxyz[2+4*(k1+l1)] = zero;
   exyz[4*(k1+l1)] = zero;
   exyz[1+4*(k1+l1)] = zero;
   exyz[2+4*(k1+l1)] = zero;
   d0 = _mm512_reduce_add_pd(v_ws);
   *wf = (ws + d0)*((float) nx)*((float) ny)*((float) nz);
   d0 = _mm512_reduce_add_pd(v_wp);
   *wm = c2*(wp + d0)*((float) nx)*((float) ny)*((float) nz);
   return;
}

/*--------------------------------------------------------------------*/
void ckncemfield3(float complex fxyz[], float complex exyz[],
                  float complex ffc[], int isign, int nx, int ny,
                  int nz, int nxvh, int nyv, int nzv, int nxhd,
                  int nyhd, int nzhd) {
/* this subroutine either adds complex vector fields if isign > 0
   or copies complex vector fields if isign < 0
   includes additional smoothing
   requires KNC, fxyz, exyz, ffc need to be 64 byte aligned
   nxhd needs to be a multiple of 8
   nxvh needs to be a multiple of 2
   fxyz, exyz needs to have 4 components
local data                                                 */
   int j, k, l, nxh, nyh, nzh, nxhs, itn, k1, l1, kk, kj, ll, lj;
   int nxyhd, nxvyh; 
   float at1;
   __m512 v_at1, v_zero, v_zt1, v_zt2;
   nxh = nx/2;
   nyh = 1 > ny/2 ? 1 : ny/2;
   nzh = 1 > nz/2 ? 1 : nz/2;
   nxhs = 2*(nxh/2);
   itn = 1 > nxhs ? 1 : nxhs;
   nxyhd = nxhd*nyhd;
   nxvyh = nxvh*nyv;
   v_zero = _mm512_setzero_ps();
/* add the fields */
   if (isign > 0) {
      for (l = 1; l < nzh; l++) {
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             at1 = cimagf(ffc[j+kk+ll]); */
               v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
               v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
               v_at1 = _mm512_permute4f128_ps(v_at1,0);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                      _mm512_int2mask(13260),(__m512i)v_at1,78);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(21845),(__m512i)v_at1,177);
/*             fxyz[4*(j+kj+lj)] += exyz[4*(j+kj+lj)]*at1;     */
/*             fxyz[1+4*(j+kj+lj)] += exyz[1+4*(j+kj+lj)]*at1; */
/*             fxyz[2+4*(j+kj+lj)] += exyz[2+4*(j+kj+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+lj)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],v_zt2);
/*             fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;     */
/*             fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1; */
/*             fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*             fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;     */
/*             fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1; */
/*             fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+l1)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*             fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*             fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*             fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = cimagf(ffc[j+kk+ll]);
               fxyz[4*(j+kj+lj)] += exyz[4*(j+kj+lj)]*at1;
               fxyz[1+4*(j+kj+lj)] += exyz[1+4*(j+kj+lj)]*at1;
               fxyz[2+4*(j+kj+lj)] += exyz[2+4*(j+kj+lj)]*at1;
               fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;
               fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1;
               fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1;
               fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;
               fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1;
               fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1;
               fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
               fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
               fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
            }
         }
         k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+ll]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+ll]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+lj)] += exyz[4*(j+lj)]*at1;     */
/*          fxyz[1+4*(j+lj)] += exyz[1+4*(j+lj)]*at1; */
/*          fxyz[2+4*(j+lj)] += exyz[2+4*(j+lj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+lj)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+lj)],v_zt2);
/*          fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;     */
/*          fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1; */
/*          fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+lj)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*          fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;     */
/*          fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1; */
/*          fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+ll]);
            fxyz[4*(j+lj)] += exyz[4*(j+lj)]*at1;
            fxyz[1+4*(j+lj)] += exyz[1+4*(j+lj)]*at1;
            fxyz[2+4*(j+lj)] += exyz[2+4*(j+lj)]*at1;
            fxyz[4*(j+k1+lj)] += exyz[4*(j+k1+lj)]*at1;
            fxyz[1+4*(j+k1+lj)] += exyz[1+4*(j+k1+lj)]*at1;
            fxyz[2+4*(j+k1+lj)] += exyz[2+4*(j+k1+lj)]*at1;
            fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;
            fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1;
            fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1;
            fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      l1 = nxvyh*nzh;
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+kk]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                     _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+kj)] += exyz[4*(j+kj)]*at1;     */
/*          fxyz[1+4*(j+kj)] += exyz[1+4*(j+kj)]*at1; */
/*          fxyz[2+4*(j+kj)] += exyz[2+4*(j+kj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+kj)],v_zt2);
/*          fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;     */
/*          fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1; */
/*          fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*          fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;     */
/*          fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1; */
/*          fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxyz[4*(j+kj)] += exyz[4*(j+kj)]*at1;
            fxyz[1+4*(j+kj)] += exyz[1+4*(j+kj)]*at1;
            fxyz[2+4*(j+kj)] += exyz[2+4*(j+kj)]*at1;
            fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;
            fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1;
            fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1;
            fxyz[4*(j+kj+l1)] += exyz[4*(j+kj+l1)]*at1;
            fxyz[1+4*(j+kj+l1)] += exyz[1+4*(j+kj+l1)]*at1;
            fxyz[2+4*(j+kj+l1)] += exyz[2+4*(j+kj+l1)]*at1;
            fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                 _mm512_int2mask(15),(float *)&ffc[j+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                  _mm512_int2mask(21845),(__m512i)v_at1,177);
/*       fxyz[4*j] += exyz[4*j]*at1;     */
/*       fxyz[1+4*j] += exyz[1+4*j]*at1; */
/*       fxyz[2+4*j] += exyz[2+4*j]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*j]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*j]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*j],v_zt2);
/*       fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;     */
/*       fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1; */
/*       fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*       fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;     */
/*       fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1; */
/*       fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+l1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*       fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;     */
/*       fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1; */
/*       fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_load_ps((float *)&fxyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_fmadd_ps(v_zt1,v_at1,v_zt2);
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxyz[4*j] += exyz[4*j]*at1;
         fxyz[1+4*j] += exyz[1+4*j]*at1;
         fxyz[2+4*j] += exyz[2+4*j]*at1;
         fxyz[4*(j+k1)] += exyz[4*(j+k1)]*at1;
         fxyz[1+4*(j+k1)] += exyz[1+4*(j+k1)]*at1;
         fxyz[2+4*(j+k1)] += exyz[2+4*(j+k1)]*at1;
         fxyz[4*(j+l1)] += exyz[4*(j+l1)]*at1;
         fxyz[1+4*(j+l1)] += exyz[1+4*(j+l1)]*at1;
         fxyz[2+4*(j+l1)] += exyz[2+4*(j+l1)]*at1;
         fxyz[4*(j+k1+l1)] += exyz[4*(j+k1+l1)]*at1;
         fxyz[1+4*(j+k1+l1)] += exyz[1+4*(j+k1+l1)]*at1;
         fxyz[2+4*(j+k1+l1)] += exyz[2+4*(j+k1+l1)]*at1;
      }
   }
/* copy the fields */
   else if (isign < 0) {
      for (l = 1; l < nzh; l++) {
         ll = nxyhd*l;
         lj = nxvyh*l;
         l1 = nxvyh*nz - lj;
         for (k = 1; k < nyh; k++) {
            kk = nxhd*k;
            kj = nxvh*k;
            k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
            for (j = 0; j < nxhs; j+=2) {
/*             at1 = cimagf(ffc[j+kk+ll]); */
               v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll]);
               v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                       _mm512_int2mask(15),(float *)&ffc[j+kk+ll+8]);
               v_at1 = _mm512_permute4f128_ps(v_at1,0);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                      _mm512_int2mask(13260),(__m512i)v_at1,78);
               v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                       _mm512_int2mask(21845),(__m512i)v_at1,177);
/*             fxyz[4*(j+kj+lj)] = exyz[4*(j+kj+lj)]*at1;     */
/*             fxyz[1+4*(j+kj+lj)] = exyz[1+4*(j+kj+lj)]*at1; */
/*             fxyz[2+4*(j+kj+lj)] = exyz[2+4*(j+kj+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+lj)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+lj)],v_zt2);
/*             fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;     */
/*             fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1; */
/*             fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*             fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;     */
/*             fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1; */
/*             fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*             fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*             fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*             fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
               v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
               v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
               _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
            }
/* loop over remaining elements */
            for (j = itn; j < nxh; j++) {
               at1 = cimagf(ffc[j+kk+ll]);
               fxyz[4*(j+kj+lj)] = exyz[4*(j+kj+lj)]*at1;
               fxyz[1+4*(j+kj+lj)] = exyz[1+4*(j+kj+lj)]*at1;
               fxyz[2+4*(j+kj+lj)] = exyz[2+4*(j+kj+lj)]*at1;
               fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;
               fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1;
               fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1;
               fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;
               fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1;
               fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1;
               fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
               fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
               fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
            }
         }
         k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+ll]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+ll]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+lj)] = exyz[4*(j+lj)]*at1;     */
/*          fxyz[1+4*(j+lj)] = exyz[1+4*(j+lj)]*at1; */
/*          fxyz[2+4*(j+lj)] = exyz[2+4*(j+lj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+lj)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+lj)],v_zt2);
/*          fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;     */
/*          fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1; */
/*          fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+lj)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+lj)],v_zt2);
/*          fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;     */
/*          fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1; */
/*          fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+ll]);
            fxyz[4*(j+lj)] = exyz[4*(j+lj)]*at1;
            fxyz[1+4*(j+lj)] = exyz[1+4*(j+lj)]*at1;
            fxyz[2+4*(j+lj)] = exyz[2+4*(j+lj)]*at1;
            fxyz[4*(j+k1+lj)] = exyz[4*(j+k1+lj)]*at1;
            fxyz[1+4*(j+k1+lj)] = exyz[1+4*(j+k1+lj)]*at1;
            fxyz[2+4*(j+k1+lj)] = exyz[2+4*(j+k1+lj)]*at1;
            fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;
            fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1;
            fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1;
            fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      l1 = nxvyh*nzh;
      for (k = 1; k < nyh; k++) {
         kk = nxhd*k;
         kj = nxvh*k;
         k1 = nxvh*ny - kj;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(ffc[j+kk]); */
            v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,
                    _mm512_int2mask(15),(float *)&ffc[j+kk]);
            v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                    _mm512_int2mask(15),(float *)&ffc[j+kk+8]);
            v_at1 = _mm512_permute4f128_ps(v_at1,0);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                    _mm512_int2mask(13260),(__m512i)v_at1,78);
            v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                     _mm512_int2mask(21845),(__m512i)v_at1,177);
/*          fxyz[4*(j+kj)] = exyz[4*(j+kj)]*at1;     */
/*          fxyz[1+4*(j+kj)] = exyz[1+4*(j+kj)]*at1; */
/*          fxyz[2+4*(j+kj)] = exyz[2+4*(j+kj)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+kj)],v_zt2);
/*          fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;     */
/*          fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1; */
/*          fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*          fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;     */
/*          fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1; */
/*          fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+kj+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+kj+l1)],v_zt2);
/*          fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*          fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*          fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
            v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
            v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
            _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
         }
/* loop over remaining elements */
         for (j = itn; j < nxh; j++) {
            at1 = cimagf(ffc[j+kk]);
            fxyz[4*(j+kj)] = exyz[4*(j+kj)]*at1;
            fxyz[1+4*(j+kj)] = exyz[1+4*(j+kj)]*at1;
            fxyz[2+4*(j+kj)] = exyz[2+4*(j+kj)]*at1;
            fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;
            fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1;
            fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1;
            fxyz[4*(j+kj+l1)] = exyz[4*(j+kj+l1)]*at1;
            fxyz[1+4*(j+kj+l1)] = exyz[1+4*(j+kj+l1)]*at1;
            fxyz[2+4*(j+kj+l1)] = exyz[2+4*(j+kj+l1)]*at1;
            fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
            fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
            fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
         }
      }
      k1 = nxvh*nyh;
/* vector loop over elements in blocks of 2 */
      for (j = 0; j < nxhs; j+=2) {
/*       at1 = cimagf(ffc[j]); */
         v_at1 = _mm512_mask_loadunpacklo_ps(v_zero,_mm512_int2mask(15),
                 (float *)&ffc[j]);
         v_at1 = _mm512_mask_loadunpackhi_ps(v_at1,
                 _mm512_int2mask(15),(float *)&ffc[j+8]);
         v_at1 = _mm512_permute4f128_ps(v_at1,0);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                 _mm512_int2mask(13260),(__m512i)v_at1,78);
         v_at1 = (__m512)_mm512_mask_shuffle_epi32((__m512i)v_at1,
                  _mm512_int2mask(21845),(__m512i)v_at1,177);
/*       fxyz[4*j] = exyz[4*j]*at1;     */
/*       fxyz[1+4*j] = exyz[1+4*j]*at1; */
/*       fxyz[2+4*j] = exyz[2+4*j]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*j]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*j],v_zt2);
/*       fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;     */
/*       fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1; */
/*       fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+k1)],v_zt2);
/*       fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;     */
/*       fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1; */
/*       fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+l1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+l1)],v_zt2);
/*       fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;     */
/*       fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1; */
/*       fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1; */
         v_zt1 = _mm512_load_ps((float *)&exyz[4*(j+k1+l1)]);
         v_zt2 = _mm512_mul_ps(v_zt1,v_at1);
         _mm512_store_ps((float *)&fxyz[4*(j+k1+l1)],v_zt2);
      }
/* loop over remaining elements */
      for (j = itn; j < nxh; j++) {
         at1 = cimagf(ffc[j]);
         fxyz[4*j] = exyz[4*j]*at1;
         fxyz[1+4*j] = exyz[1+4*j]*at1;
         fxyz[2+4*j] = exyz[2+4*j]*at1;
         fxyz[4*(j+k1)] = exyz[4*(j+k1)]*at1;
         fxyz[1+4*(j+k1)] = exyz[1+4*(j+k1)]*at1;
         fxyz[2+4*(j+k1)] = exyz[2+4*(j+k1)]*at1;
         fxyz[4*(j+l1)] = exyz[4*(j+l1)]*at1;
         fxyz[1+4*(j+l1)] = exyz[1+4*(j+l1)]*at1;
         fxyz[2+4*(j+l1)] = exyz[2+4*(j+l1)]*at1;
         fxyz[4*(j+k1+l1)] = exyz[4*(j+k1+l1)]*at1;
         fxyz[1+4*(j+k1+l1)] = exyz[1+4*(j+k1+l1)]*at1;
         fxyz[2+4*(j+k1+l1)] = exyz[2+4*(j+k1+l1)]*at1;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rvxy(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nzi, int nzp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of z,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform is performed
   f[i][m][n] = (1/nx*ny*nz)*sum(f[i][k][j]*exp(-sqrt(-1)*2pi*n*j/nx)*
         exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, a forward fourier transform is performed
   f[l][k][j] = sum(f[l][m][n]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nzi = initial z index used
   nzp = number of z indices used
   nxhd = first dimension of f
   nyd,nzd = second and third dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 8
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nxhyd;
   int i, j, k, l, n, nn, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nss, nxhs, nxhhs, itn;
   float ani;
   float complex t1, t2, t3;
   __m512i v_j, v_kmr, v_m, v_n, v_it;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4, v_t5, v_ani;
   v_j = _mm512_set_epi32(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0);
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nzt = nzi + nzp - 1;
   nxhyd = nxhd*nyd;
   nxhs = 8*(nxh/8);
   nxhhs = 8*(nxhh/8);
   itn = 1 > nxhhs ? 1 : nxhhs;
   v_m = _mm512_set_epi32(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
   v_n = _mm512_set_epi32(1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L180;
/* inverse fourier transform */
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* bit-reverse array elements in x */
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd*i + nn;
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* first transform in x */
      nrx = nxyz/nxh;
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 8*(ns/8);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
/* vector loop over elements in blocks of 8 */
               for (j = 0; j < nss; j+=8) {
/*                t1 = sct[kmr*j]; */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
/*                t2 = t1*f[j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[j+k2+joff] = f[j+k1+joff] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k2+joff],v_t4);
/*                f[j+k1+joff] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = sct[kmr*j];
                  t2 = t1*f[j+k2+joff];
                  f[j+k2+joff] = f[j+k1+joff] - t2;
                  f[j+k1+joff] += t2;
               }
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      v_ani = _mm512_set1_ps(ani);
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhhs; j+=8) {
/*          t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21845),
                   v_zero,v_t3);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          t2 = conjf(f[nxh-j+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[nxh-j+joff-7]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[nxh-j+joff+1]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[j+joff] = ani*(t1 + t2); */
            v_t3 = _mm512_mul_ps(v_ani,_mm512_add_ps(v_t1,v_t2));
/*          f[nxh-j+joff] = ani*conjf(t1 - t2); */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
            v_t4 = _mm512_mul_ps(v_ani,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[j+joff],
                  _mm512_int2mask(65532),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[nxh-j+joff-7],
                  _mm512_int2mask(16383),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[nxh-j+joff+1],
                  _mm512_int2mask(16383),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[nxh-j+joff-7],v_t4);
               _mm512_packstorehi_ps((float *)&f[nxh-j+joff+1],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
            t2 = conjf(f[nxh-j+joff]);
            t1 = f[j+joff] + t2;
            t2 = (f[j+joff] - t2)*t3;
            f[j+joff] = ani*(t1 + t2);
            f[nxh-j+joff] = ani*conjf(t1 - t2);
         }
      }
      ani = 2.0*ani;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         f[nxhh+joff] = ani*conjf(f[nxhh+joff]);
         f[joff] = ani*((crealf(f[joff]) + cimagf(f[joff]))
                   + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I);
      }
/* bit-reverse array elements in y */
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd*k1 + nn;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*          f[i+k1] = f[i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[i+joff]);
            _mm512_store_ps((float *)&f[i+k1],v_t2);
/*          f[i+joff] = t1; */
            _mm512_store_ps((float *)&f[i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+k1];
            f[i+k1] = f[i+joff];
            f[i+joff] = t1;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd*(j + k1) + nn;
               j2 = nxhd*(j + k2) + nn;
               t1 = sct[kmr*j];
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+j2] = f[i+j1] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j2],v_t4);
/*                f[i+j1] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[i+j2];
                  f[i+j2] = f[i+j1] - t2;
                  f[i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
/* unscramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff + nn;
         joff += nn;
         t1 = f[k1];
         f[k1] = 0.5*(cimagf(f[joff] + t1)
                  + crealf(f[joff] - t1)*_Complex_I);
         f[joff] = 0.5*(crealf(f[joff] + t1)
                    + cimagf(f[joff] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
L180: for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* scramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd*k;
         k1 = nxhd*ny - joff + nn;
         joff += nn;
         t1 = cimagf(f[k1]) + crealf(f[k1])*_Complex_I;
         f[k1] = conjf(f[joff] - t1);
         f[joff] += t1;
      }
/* bit-reverse array elements in y */
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd*k1 + nn;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+k1]);
/*          f[i+k1] = f[i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[i+joff]);
            _mm512_store_ps((float *)&f[i+k1],v_t2);
/*          f[i+joff] = t1; */
            _mm512_store_ps((float *)&f[i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+k1];
            f[i+k1] = f[i+joff];
            f[i+joff] = t1;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd*(j + k1) + nn;
               j2 = nxhd*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+j2] = f[i+j1] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j2],v_t4);
/*                f[i+j1] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[i+j2];
                  f[i+j2] = f[i+j1] - t2;
                  f[i+j1] += t2;
               }
            }
         }
         ns = ns2;
      }
/* scramble coefficients */
      kmr = nxyz/nx;
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
/* vector loop over elements in blocks of 8 */
         for (j = 0; j < nxhhs; j+=8) {
/*          t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          t2 = conjf(f[nxh-j+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[nxh-j+joff-7]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[nxh-j+joff+1]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[j+joff] = t1 + t2; */
            v_t3 = _mm512_add_ps(v_t1,v_t2);
/*          f[nxh-j+joff] = conjf(t1 - t2); */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[j+joff],
                  _mm512_int2mask(65532),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[nxh-j+joff-7],
                  _mm512_int2mask(16383),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[nxh-j+joff+1],
                  _mm512_int2mask(16383),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[nxh-j+joff-7],v_t4);
               _mm512_packstorehi_ps((float *)&f[nxh-j+joff+1],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
            t2 = conjf(f[nxh-j+joff]);
            t1 = f[j+joff] + t2;
            t2 = (f[j+joff] - t2)*t3;
            f[j+joff] = t1 + t2;
            f[nxh-j+joff] = conjf(t1 - t2);
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nxhd*k + nn;
         f[nxhh+joff] = 2.0*conjf(f[nxhh+joff]);
         f[joff] = (crealf(f[joff]) + cimagf(f[joff]))
                   + (crealf(f[joff]) - cimagf(f[joff]))*_Complex_I;
      }
/* bit-reverse array elements in x */
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd*i + nn;
            t1 = f[j1+joff];
            f[j1+joff] = f[j+joff];
            f[j+joff] = t1;
         }
      }
/* finally transform in x */
      nrx = nxyz/nxh;
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 8*(ns/8);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd*i + nn;
/* vector loop over elements in blocks of 8 */
               for (j = 0; j < nss; j+=8) {
/*                t1 = conjf(sct[kmr*j]); */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
                  v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(43690),
                         v_zero,v_t1);
/*                t2 = t1*f[j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[j+k2+joff] = f[j+k1+joff] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k2+joff],v_t4);
/*                f[j+k1+joff] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = conjf(sct[kmr*j]);
                  t2 = t1*f[j+k2+joff];
                  f[j+k2+joff] = f[j+k1+joff] - t2;
                  f[j+k1+joff] += t2;
               }
            }
         }
         ns = ns2;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rxz(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nyi, int nyp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of a three dimensional real to
   complex fast fourier transform and its inverse, for a subset of y,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, an inverse fourier transform is performed
   f[l][k][j] = sum(f[i][k][j]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, a forward fourier transform is performed
   f[i][m][n] = sum(f[l][m][n]*exp(sqrt(-1)*2pi*l*i/nz))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = first dimension of f
   nyd,nzd = second and third dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 8
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nxhyd;
   int i, j, k, l, n, ll, j1, j2, k1, k2, l1, ns, ns2, km, kmr, i0, i1;
   int nss, nxhs;
   float complex t1, t2;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nzh = nz/2;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nyt = nyi + nyp - 1;
   nxhyd = nxhd*nyd;
   nxhs = 8*(nxh/8);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L100;
/* inverse fourier transform */
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*          f[i+i1] = f[i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[i+i0]);
            _mm512_store_ps((float *)&f[i+i1],v_t2);
/*          f[i+i0] = t1; */
            _mm512_store_ps((float *)&f[i+i0],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+i1];
            f[i+i1] = f[i+i0];
            f[i+i0] = t1;
         }
      }
   }
/* finally transform in z */
   nrz = nxyz/nz;
   ns = 1;
   for (l = 0; l < indz; l++) {
      ns2 = ns + ns;
      km = nzh/ns;
      kmr = km*nrz;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = nxhyd*(j + k1);
            j2 = nxhyd*(j + k2);
            t1 = sct[kmr*j];
            v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                   crealf(t1));
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd*n;
               i0 = i1 + j1;
               i1 += j2;
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+i1] = f[i+i0] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i1],v_t4);
/*                f[i+i0] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[i+i1];
                  f[i+i1] = f[i+i0] - t2;
                  f[i+i0] += t2;
               }
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         t1 = f[l1];
         f[l1] = 0.5*(cimagf(f[ll] + t1)
                    + crealf(f[ll] - t1)*_Complex_I);
         f[ll] = 0.5*(crealf(f[ll] + t1)
                    + cimagf(f[ll] - t1)*_Complex_I);
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         i1 = nxhd*nyh;
         i0 = i1 + ll;
         i1 += l1;
         t1 = f[i1];
         f[i1] = 0.5*(cimagf(f[i0] + t1)
                  +   crealf(f[i0] - t1)*_Complex_I);
         f[i0] = 0.5*(crealf(f[i0] + t1)
                    + cimagf(f[i0] - t1)*_Complex_I);
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L100: for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         t1 = cimagf(f[l1]) + crealf(f[l1])*_Complex_I;
         f[l1] = conjf(f[ll] - t1);
         f[ll] += t1;
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         i1 = nxhd*nyh;
         i0 = i1 + ll;
         i1 += l1;
         t1 = cimagf(f[i1]) + crealf(f[i1])*_Complex_I;
         f[i1] = conjf(f[i0] - t1);
         f[i0] += t1;
      }
   }
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 8 */
         for (i = 0; i < nxhs; i+=8) {
/*          t1 = f[i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[i+i1]);
/*          f[i+i1] = f[i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[i+i0]);
            _mm512_store_ps((float *)&f[i+i1],v_t2);
/*          f[i+i0] = t1; */
            _mm512_store_ps((float *)&f[i+i0],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[i+i1];
            f[i+i1] = f[i+i0];
            f[i+i0] = t1;
         }
      }
   }
/* first transform in z */
   nrz = nxyz/nz;
   ns = 1;
   for (l = 0; l < indz; l++) {
      ns2 = ns + ns;
      km = nzh/ns;
      kmr = km*nrz;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = nxhyd*(j + k1);
            j2 = nxhyd*(j + k2);
            t1 = conjf(sct[kmr*j]);
            v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                   crealf(t1));
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd*n;
               i0 = i1 + j1;
               i1 += j2;
/* vector loop over elements in blocks of 8 */
               for (i = 0; i < nxhs; i+=8) {
/*                t2 = t1*f[i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[i+i1] = f[i+i0] - t2; */
                  v_t3 = _mm512_load_ps((float *)&f[i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i1],v_t4);
/*                f[i+i0] += t2; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[i+i1];
                  f[i+i1] = f[i+i0] - t2;
                  f[i+i0] += t2;
               }
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rv3xy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int nzi, int nzp, int nxhd, int nyd, int nzd,
                   int nxhyzd, int nxyzhd) {
/* this subroutine performs the x-y part of 3 three dimensional complex
   to real fast fourier transforms and their inverses, for a subset of z,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms is performed
   f[i][m][n][0:2] = (1/nx*ny*nz)*sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*n*j/nx)
         *exp(-sqrt(-1)*2pi*m*k/ny))
   if isign = 1, three forward fourier transforms are performed
   f[l][k][j][0:2] = sum(f[l][m][n][0:2]*exp(sqrt(-1)*2pi*n*j/nx)*
         exp(sqrt(-1)*2pi*m*k/ny))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nzi = initial z index used
   nzp = number of z indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:2] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:2] = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:2] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:2] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:2]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:2]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:2]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:2]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 2
   f needs to have 4 components
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, nxhh, ny, nyh;
   int nz, nxyz, nxhyz, nzt, nrx, nry, nxhd4, nxhyd;
   int i, j, k, l, n, nn, jj, j1, j2, k1, k2, ns, ns2, km, kmr, joff;
   int nss, nxhs, nxhhs, itn;
   float at1, at2, ani;
   float complex t1, t2, t3, t4;
   __m512i v_j, v_kmr, v_m, v_n, v_l, v_it;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4, v_t5, v_ani, v_half;
   v_j = _mm512_set_epi32(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   nxhh = nx/4;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nzt = nzi + nzp - 1;
   nxhd4 = 4*nxhd;
   nxhyd = nxhd4*nyd;
   nxhs = 2*(nxh/2);
   nxhhs = 2*(nxhh/2);
   itn = 1 > nxhhs ? 1 : nxhhs;
   v_m = _mm512_set_epi32(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
   v_n = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   v_half = _mm512_set1_ps(0.5f);
   if (isign > 0)
      goto L230;
/* inverse fourier transform */
   v_l = _mm512_set_epi32(15,11,14,10,13,9,12,8,7,3,6,2,5,1,4,0);
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          at1 = cimagf(f[2+4*j+joff]);                          */
/*          at2 = crealf(f[2+4*j+joff]);                          */
/*          f[2+4*j+joff] = crealf(f[1+4*j+joff])                 */
/*                          + crealf(f[3+4*j+joff])*_Complex_I;   */
/*          f[1+4*j+joff] = cimagf(f[4*j+joff]) + at1*_Complex_I; */
/*          f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;   */
            v_t1 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = (__m512)_mm512_permutevar_epi32(v_l,(__m512i)v_t1);
            _mm512_store_ps((float *)&f[4*j+joff],v_t1);
          }
/* loop over remaining elements */
         for (j = nxhs; j < nxh; j++) {
            at1 = cimagf(f[2+4*j+joff]);
            at2 = crealf(f[2+4*j+joff]);
            f[2+4*j+joff] = crealf(f[1+4*j+joff])
                            + crealf(f[3+4*j+joff])*_Complex_I;
            f[1+4*j+joff] = cimagf(f[4*j+joff]) + at1*_Complex_I;
            f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
          }
      }
/* bit-reverse array elements in x */
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd4*i + nn;
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff]);
            v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff]);
            v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
            _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
               _mm512_int2mask(255),v_t2);
            _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
               _mm512_int2mask(255),v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
               _mm512_int2mask(255),v_t1);
            _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
               _mm512_int2mask(255),v_t1);
         }
      }
/* first transform in x */
      nrx = nxyz/nxh;
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 2*(ns/2);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nss; j+=2) {
/*                t1 = sct[kmr*j]; */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
/*                t2 = t1*f[4*j+k2+joff];   */
/*                t3 = t1*f[1+4*j+k2+joff]; */
/*                t4 = t1*f[2+4*j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*                f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
/*                f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k2+joff],v_t4);
/*                f[4*j+k1+joff] += t2;   */
/*                f[1+4*j+k1+joff] += t3; */
/*                f[2+4*j+k1+joff] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = sct[kmr*j];
                  t2 = t1*f[4*j+k2+joff];
                  t3 = t1*f[1+4*j+k2+joff];
                  t4 = t1*f[2+4*j+k2+joff];
                  f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
                  f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
                  f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
                  f[4*j+k1+joff] += t2;
                  f[1+4*j+k1+joff] += t3;
                  f[2+4*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble coefficients and normalize */
      kmr = nxyz/nx;
      ani = 0.5/(((float) nx)*((float) ny)*((float) nz));
      v_ani = _mm512_set1_ps(ani);
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhhs; j+=2) {
/*          t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21845),
                   v_zero,v_t3);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          for (jj = 0; jj < 3; jj++) {      */
/*          t2 = conjf(f[jj+4*(nxh-j)+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff+8]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[jj+4*j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[jj+4*j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);

/*          f[jj+4*j+joff] = ani*(t1 + t2); */
            v_t3 = _mm512_mul_ps(v_ani,_mm512_add_ps(v_t1,v_t2));
/*          f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2); */
/*          }                                          */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
            v_t4 = _mm512_mul_ps(v_ani,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(65280),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],
                  _mm512_int2mask(255),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],
                  _mm512_int2mask(255),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[4*j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],v_t4);
               _mm512_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) - crealf(sct[kmr*j])*_Complex_I;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+4*(nxh-j)+joff]);
               t1 = f[jj+4*j+joff] + t2;
               t2 = (f[jj+4*j+joff] - t2)*t3;
               f[jj+4*j+joff] = ani*(t1 + t2);
               f[jj+4*(nxh-j)+joff] = ani*conjf(t1 - t2);
            }
         }
      }
/*    ani = 2.0*ani; */
      v_ani = _mm512_add_ps(v_ani,v_ani);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/*       for (jj = 0; jj < 3; jj++) {                      */
/*       f[jj+4*nxhh+joff] = ani*conjf(f[jj+4*nxhh+joff]); */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[4*nxhh+joff]);
         v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(42),v_zero,
                v_t1);
         v_t1 = _mm512_mul_ps(v_ani,v_t1);
         _mm512_mask_store_ps((float *)&f[4*nxhh+joff],
               _mm512_int2mask(63),v_t1);
/*       f[jj+joff] = ani*((crealf(f[jj+joff])            */
/*                     + cimagf(f[jj+joff]))              */
/*                     + (crealf(f[jj+joff])              */
/*                     - cimagf(f[jj+joff]))*_Complex_I); */
/*       }                                                */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
         v_t3 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t1,v_t2);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(21),v_t1,v_t2);
         v_t3 = _mm512_mul_ps(v_ani,v_t3);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t3);
      }
/* bit-reverse array elements in y */
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*          f[4*i+k1] = f[4*i+joff];     */
/*          f[1+4*i+k1] = f[1+4*i+joff]; */
/*          f[2+4*i+k1] = f[2+4*i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
            _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*          f[4*i+joff] = t1;   */
/*          f[1+4*i+joff] = t2; */
/*          f[2+4*i+joff] = t3; */
            _mm512_store_ps((float *)&f[4*i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[4*i+k1];
            t2 = f[1+4*i+k1];
            t3 = f[2+4*i+k1];
            f[4*i+k1] = f[4*i+joff];
            f[1+4*i+k1] = f[1+4*i+joff];
            f[2+4*i+k1] = f[2+4*i+joff];
            f[4*i+joff] = t1;
            f[1+4*i+joff] = t2;
            f[2+4*i+joff] = t3;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd4*(j + k1) + nn;
               j2 = nxhd4*(j + k2) + nn;
               t1 = sct[kmr*j];
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+j2];   */
/*                t3 = t1*f[1+4*i+j2]; */
/*                t4 = t1*f[2+4*i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+j2] = f[4*i+j1] - t2;     */
/*                f[1+4*i+j2] = f[1+4*i+j1] - t3; */
/*                f[2+4*i+j2] = f[2+4*i+j1] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j2],v_t4);
/*                f[4*i+j1] += t2;   */
/*                f[1+4*i+j1] += t3; */
/*                f[2+4*i+j1] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+j2];
                  t3 = t1*f[1+4*i+j2];
                  t4 = t1*f[2+4*i+j2];
                  f[4*i+j2] = f[4*i+j1] - t2;
                  f[1+4*i+j2] = f[1+4*i+j1] - t3;
                  f[2+4*i+j2] = f[2+4*i+j1] - t4;
                  f[4*i+j1] += t2;
                  f[1+4*i+j1] += t3;
                  f[2+4*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* unscramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd4*k;
         k1 = nxhd4*ny - joff + nn;
         joff += nn;
/*       for (jj = 0; jj < 3; jj++) { */
/*       t1 = f[jj+k1];               */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[k1]);
/*       f[jj+k1] = 0.5*(cimagf(f[jj+joff] + t1)            */
/*                   + crealf(f[jj+joff] - t1)*_Complex_I); */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(42),v_t2,v_t1);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(21),v_t2,v_t1);
         v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
         v_t3 = _mm512_mul_ps(v_half,v_t3);
         _mm512_mask_store_ps((float *)&f[k1],_mm512_int2mask(63),v_t3);
/*       f[jj+joff] = 0.5*(crealf(f[jj+joff] + t1)            */
/*                     + cimagf(f[jj+joff] - t1)*_Complex_I); */
/*       }                                                    */
         v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t2,v_t1);
         v_t2 = _mm512_mask_add_ps(v_t2,_mm512_int2mask(21),v_t2,v_t1);
         v_t2 = _mm512_mul_ps(v_half,v_t2);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),v_t2);
      }
   }
   return;
/* forward fourier transform */
L230: v_l = _mm512_set_epi32(15,13,11,9,14,12,10,8,7,5,3,1,6,4,2,0);
   for (n = nzi-1; n < nzt; n++) {
      nn = nxhyd*n;
/* scramble modes kx = 0, nx/2 */
      for (k = 1; k < nyh; k++) {
         joff = nxhd4*k;
         k1 = nxhd4*ny - joff + nn;
         joff += nn;
/*       for (jj = 0; jj < 3; jj++) {                            */
/*       t1 = cimagf(f[jj+k1]) + crealf(f[jj+k1])*_Complex_I; */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[k1]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,177);
/*       f[jj+k1] = conjf(f[jj+joff] - t1); */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(63),v_t2,v_t1);
         v_t3 = _mm512_mask_sub_ps(v_t3,_mm512_int2mask(42),
                v_zero,v_t3);
         _mm512_mask_store_ps((float *)&f[k1],_mm512_int2mask(63),v_t3);
/*       f[jj+joff] += t1; */
/*       }                 */
         v_t2 = _mm512_mask_add_ps(v_t2,_mm512_int2mask(63),v_t2,v_t1);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t2);
      }
/* bit-reverse array elements in y */
      nry = nxhyz/ny;
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
         k1 = (mixup[k] - 1)/nry;
         if (k >= k1)
            continue;
         k1 = nxhd4*k1 + nn;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+k1];   */
/*          t2 = f[1+4*i+k1]; */
/*          t3 = f[2+4*i+k1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+k1]);
/*          f[4*i+k1] = f[4*i+joff];     */
/*          f[1+4*i+k1] = f[1+4*i+joff]; */
/*          f[2+4*i+k1] = f[2+4*i+joff]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+joff]);
            _mm512_store_ps((float *)&f[4*i+k1],v_t2);
/*          f[4*i+joff] = t1;   */
/*          f[1+4*i+joff] = t2; */
/*          f[2+4*i+joff] = t3; */
            _mm512_store_ps((float *)&f[4*i+joff],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[4*i+k1];
            t2 = f[1+4*i+k1];
            t3 = f[2+4*i+k1];
            f[4*i+k1] = f[4*i+joff];
            f[1+4*i+k1] = f[1+4*i+joff];
            f[2+4*i+k1] = f[2+4*i+joff];
            f[4*i+joff] = t1;
            f[1+4*i+joff] = t2;
            f[2+4*i+joff] = t3;
         }
      }
/* then transform in y */
      nry = nxyz/ny;
      ns = 1;
      for (l = 0; l < indy; l++) {
         ns2 = ns + ns;
         km = nyh/ns;
         kmr = km*nry;
         for (k = 0; k < km; k++) {
            k1 = ns2*k;
            k2 = k1 + ns;
            for (j = 0; j < ns; j++) {
               j1 = nxhd4*(j + k1) + nn;
               j2 = nxhd4*(j + k2) + nn;
               t1 = conjf(sct[kmr*j]);
               v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                      crealf(t1));
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+j2];   */
/*                t3 = t1*f[1+4*i+j2]; */
/*                t4 = t1*f[2+4*i+j2]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+j2]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+j2] = f[4*i+j1] - t2;     */
/*                f[1+4*i+j2] = f[1+4*i+j1] - t3; */
/*                f[2+4*i+j2] = f[2+4*i+j1] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+j1]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j2],v_t4);
/*                f[4*i+j1] += t2;   */
/*                f[1+4*i+j1] += t3; */
/*                f[2+4*i+j1] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+j1],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+j2];
                  t3 = t1*f[1+4*i+j2];
                  t4 = t1*f[2+4*i+j2];
                  f[4*i+j2] = f[4*i+j1] - t2;
                  f[1+4*i+j2] = f[1+4*i+j1] - t3;
                  f[2+4*i+j2] = f[2+4*i+j1] - t4;
                  f[4*i+j1] += t2;
                  f[1+4*i+j1] += t3;
                  f[2+4*i+j1] += t4;
               }
            }
         }
         ns = ns2;
      }
/* scramble coefficients */
      kmr = nxyz/nx;
      v_kmr = _mm512_set1_epi32(2*kmr);
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhhs; j+=2) {
/*          t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I; */
            v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
            v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
            v_t3 = _mm512_i32gather_ps(v_it,(float *)sct,4);
            v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,177);
/*          for (jj = 0; jj < 3; jj++) {      */
/*          t2 = conjf(f[jj+4*(nxh-j)+joff]); */
            v_t2 = _mm512_loadunpacklo_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff]);
            v_t2 = _mm512_loadunpackhi_ps(v_t2,
                   (float *)&f[4*(nxh-j-1)+joff+8]);
/* reverse data */
            v_t2 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t2);
            v_t2 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(43690),
                   v_zero,v_t2);
/*          t1 = f[jj+4*j+joff] + t2; */
            v_t4 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = _mm512_add_ps(v_t4,v_t2);
/*          t2 = (f[jj+4*j+joff] - t2)*t3; */
            v_t2 = _mm512_sub_ps(v_t4,v_t2);
            v_t5 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,160);
            v_t5 = _mm512_mul_ps(v_t2,v_t5);
            v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
            v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t3,245);
            v_t4 = _mm512_mul_ps(v_t2,v_t4);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                   v_zero,v_t4);
            v_t2 = _mm512_add_ps(v_t5,v_t4);
/*          f[jj+4*j+joff] = t1 + t2; */
            v_t3 = _mm512_add_ps(v_t1,v_t2);
/*          f[jj+4*(nxh-j)+joff] = conjf(t1 - t2); */
/*          }                                      */
            v_t4 = _mm512_sub_ps(v_t1,v_t2);
            v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(43690),
                   v_zero,v_t4);
/* reverse data */
            v_t4 = (__m512)_mm512_permutevar_epi32(v_n,(__m512i)v_t4);
            if (j==0) {
               _mm512_mask_store_ps((float *)&f[4*j+joff],
                  _mm512_int2mask(65280),v_t3);
               _mm512_mask_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],
                  _mm512_int2mask(255),v_t4);
               _mm512_mask_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],
                  _mm512_int2mask(255),v_t4);
            }
            else {
               _mm512_store_ps((float *)&f[4*j+joff],v_t3);
               _mm512_packstorelo_ps((float *)&f[4*(nxh-j-1)+joff],v_t4);
               _mm512_packstorehi_ps((float *)&f[4*(nxh-j-1)+joff+8],v_t4);
            }
         }
/* loop over remaining elements */
         for (j = itn; j < nxhh; j++) {
            t3 = cimagf(sct[kmr*j]) + crealf(sct[kmr*j])*_Complex_I;
            for (jj = 0; jj < 3; jj++) {
               t2 = conjf(f[jj+4*(nxh-j)+joff]);
               t1 = f[jj+4*j+joff] + t2;
               t2 = (f[jj+4*j+joff] - t2)*t3;
               f[jj+4*j+joff] = t1 + t2;
               f[jj+4*(nxh-j)+joff] = conjf(t1 - t2);
            }
         }
      }
      for (k = 0; k < ny; k++) {
         joff = nxhd4*k + nn;
/*       for (jj = 0; jj < 3; jj++) { */
/*       f[jj+4*nxhh+joff] = 2.0*conjf(f[jj+4*nxhh+joff]); */
         v_t1 = _mm512_mask_load_ps(v_t1,_mm512_int2mask(63),
                (float *)&f[4*nxhh+joff]);
         v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(42),v_zero,
                v_t1);
         v_t1 = _mm512_add_ps(v_t1,v_t1);
         _mm512_mask_store_ps((float *)&f[4*nxhh+joff],
               _mm512_int2mask(63),v_t1);
/*       f[jj+joff] = (crealf(f[jj+joff]) + cimagf(f[jj+joff])) */
/*                  + (crealf(f[jj+joff])                       */
/*                  - cimagf(f[jj+joff]))*_Complex_I;           */
/*       }                                                      */
         v_t2 = _mm512_mask_load_ps(v_t2,_mm512_int2mask(63),
                (float *)&f[joff]);
         v_t1 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
         v_t3 = _mm512_mask_sub_ps(v_t2,_mm512_int2mask(42),v_t1,v_t2);
         v_t3 = _mm512_mask_add_ps(v_t3,_mm512_int2mask(21),v_t1,v_t2);
         _mm512_mask_store_ps((float *)&f[joff],_mm512_int2mask(63),
         v_t3);
      }
/* bit-reverse array elements in x */
      nrx = nxhyz/nxh;
      for (j = 0; j < nxh; j++) {
         j1 = (mixup[j] - 1)/nrx;
         if (j >= j1)
            continue;
         for (i = 0; i < ny; i++) {
            joff = nxhd4*i + nn;
/*          t1 = f[4*j1+joff];   */
/*          t2 = f[1+4*j1+joff]; */
/*          t3 = f[2+4*j1+joff]; */
            v_t1 = _mm512_mask_loadunpacklo_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff]);
            v_t1 = _mm512_mask_loadunpackhi_ps(v_t1,
                   _mm512_int2mask(255),(float *)&f[4*j1+joff+8]);
/*          f[4*j1+joff] = f[4*j+joff];     */
/*          f[1+4*j1+joff] = f[1+4*j+joff]; */
/*          f[2+4*j1+joff] = f[2+4*j+joff]; */
            v_t2 = _mm512_mask_loadunpacklo_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff]);
            v_t2 = _mm512_mask_loadunpackhi_ps(v_t2,
                   _mm512_int2mask(255),(float *)&f[4*j+joff+8]);
            _mm512_mask_packstorelo_ps((float *)&f[4*j1+joff],
               _mm512_int2mask(255),v_t2);
            _mm512_mask_packstorehi_ps((float *)&f[4*j1+joff+8],
               _mm512_int2mask(255),v_t2);
/*          f[4*j+joff] = t1;   */
/*          f[1+4*j+joff] = t2; */
/*          f[2+4*j+joff] = t3; */
            _mm512_mask_packstorelo_ps((float *)&f[4*j+joff],
               _mm512_int2mask(255),v_t1);
            _mm512_mask_packstorehi_ps((float *)&f[4*j+joff+8],
               _mm512_int2mask(255),v_t1);
         }
      }
/* finally transform in x */
      nrx = nxyz/nxh;
      ns = 1;
      for (l = 0; l < indx1; l++) {
         ns2 = ns + ns;
         km = nxhh/ns;
         kmr = km*nrx;
         nss = 2*(ns/2);
         v_kmr = _mm512_set1_epi32(2*kmr);
         for (k = 0; k < km; k++) {
            k1 = 4*ns2*k;
            k2 = k1 + 4*ns;
            for (i = 0; i < ny; i++) {
               joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
               for (j = 0; j < nss; j+=2) {
/*                t1 = conjf(sct[kmr*j]); */
                  v_it = _mm512_add_epi32(_mm512_set1_epi32(j),v_j);
                  v_it = _mm512_fmadd_epi32(v_kmr,v_it,v_m);
                  v_t1 = _mm512_i32gather_ps(v_it,(float *)sct,4);
                  v_t1 = _mm512_mask_sub_ps(v_t1,_mm512_int2mask(43690),
                         v_zero,v_t1);
/*                t2 = t1*f[4*j+k2+joff];   */
/*                t3 = t1*f[1+4*j+k2+joff]; */
/*                t4 = t1*f[2+4*j+k2+joff]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*j+k2+joff]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*j+k2+joff] = f[4*j+k1+joff] - t2;     */
/*                f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3; */
/*                f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*j+k1+joff]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k2+joff],v_t4);
/*                f[4*j+k1+joff] += t2;   */
/*                f[1+4*j+k1+joff] += t3; */
/*                f[2+4*j+k1+joff] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*j+k1+joff],v_t4);
               }
/* loop over remaining elements */
               for (j = nss; j < ns; j++) {
                  t1 = conjf(sct[kmr*j]);
                  t2 = t1*f[4*j+k2+joff];
                  t3 = t1*f[1+4*j+k2+joff];
                  t4 = t1*f[2+4*j+k2+joff];
                  f[4*j+k2+joff] = f[4*j+k1+joff] - t2;
                  f[1+4*j+k2+joff] = f[1+4*j+k1+joff] - t3;
                  f[2+4*j+k2+joff] = f[2+4*j+k1+joff] - t4;
                  f[4*j+k1+joff] += t2;
                  f[1+4*j+k1+joff] += t3;
                  f[2+4*j+k1+joff] += t4;
               }
            }
         }
         ns = ns2;
      }
/* swap complex components */
      for (i = 0; i < ny; i++) {
         joff = nxhd4*i + nn;
/* vector loop over elements in blocks of 2 */
         for (j = 0; j < nxhs; j+=2) {
/*          f[3+4*j+joff] = cimagf(f[2+4*j+joff]) */
/*                          + cimagf(f[3+4*j+joff])*_Complex_I; */
/*          at1 = crealf(f[2+4*j+joff]); */
/*          f[2+4*j+joff] = cimagf(f[4*j+joff]) */
/*                          + cimagf(f[1+4*j+joff])*_Complex_I; */
/*          at2 = crealf(f[1+4*j+joff]); */
/*          f[1+4*j+joff] = at1 + 0.0*_Complex_I; */
/*          f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I; */
            v_t1 = _mm512_load_ps((float *)&f[4*j+joff]);
            v_t1 = (__m512)_mm512_permutevar_epi32(v_l,(__m512i)v_t1);
            _mm512_store_ps((float *)&f[4*j+joff],v_t1);
         }
/* loop over remaining elements */
         for (j = nxhs; j < nxh; j++) {
            f[3+4*j+joff] = cimagf(f[2+4*j+joff])
                            + cimagf(f[3+4*j+joff])*_Complex_I;
            at1 = crealf(f[2+4*j+joff]);
            f[2+4*j+joff] = cimagf(f[4*j+joff])
                            + cimagf(f[1+4*j+joff])*_Complex_I;
            at2 = crealf(f[1+4*j+joff]);
            f[1+4*j+joff] = at1 + 0.0*_Complex_I;
            f[4*j+joff] = crealf(f[4*j+joff]) + at2*_Complex_I;
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncfft3rv3z(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nyi, int nyp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd) {
/* this subroutine performs the z part of 3 three dimensional complex to
   real fast fourier transforms and their inverses, for a subset of y,
   using complex arithmetic
   for isign = (-1,1), input: all, output: f
   for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
   for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
   where N = (nx/2)*ny*nz
   indx/indy/indz = exponent which determines length in x/y/z direction,
   where nx=2**indx, ny=2**indy, nz=2**indz
   if isign = -1, three inverse fourier transforms is performed
   f[l][k][j][0:2] = sum(f[i][k][j][0:2]*exp(-sqrt(-1)*2pi*l*i/nz))
   if isign = 1, three forward fourier transforms are performed
   f[i][m][n][0:2] = sum(f[l][m][n][0:2]*exp(sqrt(-1)*2pi*l*i/nz))
   mixup = array of bit reversed addresses
   sct = sine/cosine table
   nyi = initial y index used
   nyp = number of y indices used
   nxhd = second dimension of f
   nyd,nzd = third and fourth dimensions of f
   nxhyzd = maximum of (nx/2,ny,nz)
   nxyzhd = maximum of (nx,ny,nz)/2
   fourier coefficients are stored as follows:
   f[l][k][j][0:2] = real, imaginary part of mode j,k,l
   where 0 <= j < nx/2, 0 <= k < ny, 0 <= l < nz, except for
   f[l][k][0][0:2], = real, imaginary part of mode nx/2,k,l,
   where ny/2+1 <= k < ny and 0 <= l < nz, and
   f[l][0][0][0:2] = real, imaginary part of mode nx/2,0,l,
   f[l][ny/2][0][0:2] = real, imaginary part mode nx/2,ny/2,l,
   where nz/2+1 <= l < nz, and
   imag(f[0][0][0][0:2]) = real part of mode nx/2,0,0
   imag(f[0][ny/2][0][0:2]) = real part of mode nx/2,ny/2,0
   imag(f[nz/2][0][0][0:2]) = real part of mode nx/2,0,nz/2
   imag(f[nz/2][ny/2][0][0:2]) = real part of mode nx/2,ny/2,nz/2
   using jpl storage convention, as described in:
   E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
   Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
   Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
   December 1993.
   requires KNC, f needs to be 64 byte aligned
   nxhd need to be a multiple of 2
   f needs to have 4 components
   written by viktor k. decyk, ucla
local data                                                            */
   int indx1, ndx1yz, nx, nxh, ny, nyh;
   int nz, nzh, nxyz, nxhyz, nyt, nrz, nxhd4, nxhyd;
   int i, j, k, l, n, ll, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr;
   int i0, i1;
   int nxhs;
   float complex t1, t2, t3, t4;
   __m512 v_zero, v_t1, v_t2, v_t3, v_t4;
   if (isign==0)
      return;
   indx1 = indx - 1;
   ndx1yz = indx1 > indy ? indx1 : indy;
   ndx1yz = ndx1yz > indz ? ndx1yz : indz;
   nx = 1L<<indx;
   nxh = nx/2;
   ny = 1L<<indy;
   nyh = ny/2;
   nz = 1L<<indz;
   nzh = nz/2;
   nxyz = nx > ny ? nx : ny;
   nxyz = nxyz > nz ? nxyz : nz;
   nxhyz = 1L<<ndx1yz;
   nyt = nyi + nyp - 1;
   nxhd4 = 4*nxhd;
   nxhyd = nxhd4*nyd;
   nxhs = 2*(nxh/2);
   v_zero = _mm512_setzero_ps();
   v_t1 = _mm512_setzero_ps();
   v_t2 = _mm512_setzero_ps();
   v_t3 = _mm512_setzero_ps();
   v_t4 = _mm512_setzero_ps();
   if (isign > 0)
      goto L120;
/* inverse fourier transform */
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd4*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+i1];   */
/*          t2 = f[1+4*i+i1]; */
/*          t3 = f[2+4*i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*          f[4*i+i1] = f[4*i+i0];     */
/*          f[1+4*i+i1] = f[1+4*i+i0]; */
/*          f[2+4*i+i1] = f[2+4*i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
            _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*          f[4*i+i0] = t1;   */
/*          f[1+4*i+i0] = t2; */
/*          f[2+4*i+i0] = t3; */
            _mm512_store_ps((float *)&f[4*i+i0],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[4*i+i1];
            t2 = f[1+4*i+i1];
            t3 = f[2+4*i+i1];
            f[4*i+i1] = f[4*i+i0];
            f[1+4*i+i1] = f[1+4*i+i0];
            f[2+4*i+i1] = f[2+4*i+i0];
            f[4*i+i0] = t1;
            f[1+4*i+i0] = t2;
            f[2+4*i+i0] = t3;
         }
      }
   }
/* finally transform in z */
   nrz = nxyz/nz;
   ns = 1;
   for (l = 0; l < indz; l++) {
      ns2 = ns + ns;
      km = nzh/ns;
      kmr = km*nrz;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = nxhyd*(j + k1);
            j2 = nxhyd*(j + k2);
            t1 = sct[kmr*j];
            v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                   crealf(t1));
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd4*n;
               i0 = i1 + j1;
               i1 += j2;
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+i1];   */
/*                t3 = t1*f[1+4*i+i1]; */
/*                t4 = t1*f[2+4*i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+i1] = f[4*i+i0] - t2;     */
/*                f[1+4*i+i1] = f[1+4*i+i0] - t3; */
/*                f[2+4*i+i1] = f[2+4*i+i0] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i1],v_t4);
/*                f[4*i+i0] += t2;   */
/*                f[1+4*i+i0] += t3; */
/*                f[2+4*i+i0] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+i1];
                  t3 = t1*f[1+4*i+i1];
                  t4 = t1*f[2+4*i+i1];
                  f[4*i+i1] = f[4*i+i0] - t2;
                  f[1+4*i+i1] = f[1+4*i+i0] - t3;
                  f[2+4*i+i1] = f[2+4*i+i0] - t4;
                  f[4*i+i0] += t2;
                  f[1+4*i+i0] += t3;
                  f[2+4*i+i0] += t4;
               }
            }
         }
      }
      ns = ns2;
   }
/* unscramble modes kx = 0, nx/2 */
   for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < 3; jj++) {
            t1 = f[jj+l1];
            f[jj+l1] = 0.5*(cimagf(f[jj+ll] + t1)
                          + crealf(f[jj+ll] - t1)*_Complex_I);
            f[jj+ll] = 0.5*(crealf(f[jj+ll] + t1)
                          + cimagf(f[jj+ll] - t1)*_Complex_I);
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < 3; jj++) {
            i1 = nxhd4*nyh;
            i0 = i1 + ll;
            i1 += l1;
            t1 = f[jj+i1];
            f[jj+i1] = 0.5*(cimagf(f[jj+i0] + t1)
                        +   crealf(f[jj+i0] - t1)*_Complex_I);
            f[jj+i0] = 0.5*(crealf(f[jj+i0] + t1)
                          + cimagf(f[jj+i0] - t1)*_Complex_I);
         }
      }
   }
   return;
/* forward fourier transform */
/* scramble modes kx = 0, nx/2 */
L120: for (n = 1; n < nzh; n++) {
      ll = nxhyd*n;
      l1 = nxhyd*nz - ll;
      if (nyi==1) {
         for (jj = 0; jj < 3; jj++) {
            t1 = cimagf(f[jj+l1]) + crealf(f[jj+l1])*_Complex_I;
            f[jj+l1] = conjf(f[jj+ll] - t1);
            f[jj+ll] += t1;
         }
      }
      if ((nyi <= (nyh+1)) && (nyt >= (nyh+1))) {
         for (jj = 0; jj < 3; jj++) {
            i1 = nxhd4*nyh;
            i0 = i1 + ll;
            i1 += l1;
            t1 = cimagf(f[jj+i1]) + crealf(f[jj+i1])*_Complex_I;
            f[jj+i1] = conjf(f[jj+i0] - t1);
            f[jj+i0] += t1;
         }
      }
   }
/* bit-reverse array elements in z */
   nrz = nxhyz/nz;
   for (l = 0; l < nz; l++) {
      ll = nxhyd*l;
      l1 = (mixup[l] - 1)/nrz;
      if (l >= l1)
         continue;
      l1 = nxhyd*l1;
      for (n = nyi-1; n < nyt; n++) {
         i1 = nxhd4*n;
         i0 = i1 + ll;
         i1 += l1;
/* vector loop over elements in blocks of 2 */
         for (i = 0; i < nxhs; i+=2) {
/*          t1 = f[4*i+i1];   */
/*          t2 = f[1+4*i+i1]; */
/*          t3 = f[2+4*i+i1]; */
            v_t1 = _mm512_load_ps((float *)&f[4*i+i1]);
/*          f[4*i+i1] = f[4*i+i0];     */
/*          f[1+4*i+i1] = f[1+4*i+i0]; */
/*          f[2+4*i+i1] = f[2+4*i+i0]; */
            v_t2 = _mm512_load_ps((float *)&f[4*i+i0]);
            _mm512_store_ps((float *)&f[4*i+i1],v_t2);
/*          f[4*i+i0] = t1;   */
/*          f[1+4*i+i0] = t2; */
/*          f[2+4*i+i0] = t3; */
            _mm512_store_ps((float *)&f[4*i+i0],v_t1);
         }
/* loop over remaining elements */
         for (i = nxhs; i < nxh; i++) {
            t1 = f[4*i+i1];
            t2 = f[1+4*i+i1];
            t3 = f[2+4*i+i1];
            f[4*i+i1] = f[4*i+i0];
            f[1+4*i+i1] = f[1+4*i+i0];
            f[2+4*i+i1] = f[2+4*i+i0];
            f[4*i+i0] = t1;
            f[1+4*i+i0] = t2;
            f[2+4*i+i0] = t3;
         }
      }
   }
/* first transform in z */
nrz = nxyz/nz;
   ns = 1;
   for (l = 0; l < indz; l++) {
      ns2 = ns + ns;
      km = nzh/ns;
      kmr = km*nrz;
      for (k = 0; k < km; k++) {
         k1 = ns2*k;
         k2 = k1 + ns;
         for (j = 0; j < ns; j++) {
            j1 = nxhyd*(j + k1);
            j2 = nxhyd*(j + k2);
            t1 = conjf(sct[kmr*j]);
            v_t1 = _mm512_set4_ps(cimagf(t1),crealf(t1),cimagf(t1),
                   crealf(t1));
            for (n = nyi-1; n < nyt; n++) {
               i1 = nxhd4*n;
               i0 = i1 + j1;
               i1 += j2;
/* vector loop over elements in blocks of 2 */
               for (i = 0; i < nxhs; i+=2) {
/*                t2 = t1*f[4*i+i1];   */
/*                t3 = t1*f[1+4*i+i1]; */
/*                t4 = t1*f[2+4*i+i1]; */
                  v_t2 = _mm512_load_ps((float *)&f[4*i+i1]);
                  v_t3 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,160);
                  v_t3 = _mm512_mul_ps(v_t2,v_t3);
                  v_t2 = (__m512)_mm512_shuffle_epi32((__m512i)v_t2,177);
                  v_t4 = (__m512)_mm512_shuffle_epi32((__m512i)v_t1,245);
                  v_t4 = _mm512_mul_ps(v_t2,v_t4);
                  v_t4 = _mm512_mask_sub_ps(v_t4,_mm512_int2mask(21845),
                         v_zero,v_t4);
                  v_t2 = _mm512_add_ps(v_t3,v_t4);
/*                f[4*i+i1] = f[4*i+i0] - t2;     */
/*                f[1+4*i+i1] = f[1+4*i+i0] - t3; */
/*                f[2+4*i+i1] = f[2+4*i+i0] - t4; */
                  v_t3 = _mm512_load_ps((float *)&f[4*i+i0]);
                  v_t4 = _mm512_sub_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i1],v_t4);
/*                f[4*i+i0] += t2;   */
/*                f[1+4*i+i0] += t3; */
/*                f[2+4*i+i0] += t4; */
                  v_t4 = _mm512_add_ps(v_t3,v_t2);
                  _mm512_store_ps((float *)&f[4*i+i0],v_t4);
               }
/* loop over remaining elements */
               for (i = nxhs; i < nxh; i++) {
                  t2 = t1*f[4*i+i1];
                  t3 = t1*f[1+4*i+i1];
                  t4 = t1*f[2+4*i+i1];
                  f[4*i+i1] = f[4*i+i0] - t2;
                  f[1+4*i+i1] = f[1+4*i+i0] - t3;
                  f[2+4*i+i1] = f[2+4*i+i0] - t4;
                  f[4*i+i0] += t2;
                  f[1+4*i+i0] += t3;
                  f[2+4*i+i0] += t4;
               }
            }
         }
      }
      ns = ns2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rvx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for real to complex fft, with packed data */
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      ckncfft3rvxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rxz(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rvxy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rv3(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
/* wrapper function for 3 2d real to complex ffts, with packed data */
/* local data */
   int ny, nz;
   static int nyi = 1, nzi = 1;
/* calculate range of indices */
   ny = 1L<<indy;
   nz = 1L<<indz;
/* inverse fourier transform */
   if (isign < 0) {
/* perform xy fft */
      ckncfft3rv3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
/* perform z fft */
      ckncfft3rv3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                   nxhyzd,nxyzhd);
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform z fft */
      ckncfft3rv3z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,nzd,
                  nxhyzd,nxyzhd);
/* perform xy fft */
      ckncfft3rv3xy(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
                    nzd,nxhyzd,nxyzhd);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void ckncxiscan2_(int *isdata, int *nths) {
   ckncxiscan2(isdata,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgbpush3lt_(float *part, float *fxyz, float *bxyz, float *qbm,
                    float *dt, float *dtc, float *ek, int *idimp,
                    int *nop, int *npe, int *nx, int *ny, int *nz,
                    int *nxv, int *nyv, int *nzv, int *ipbc) {
   ckncgbpush3lt(part,fxyz,bxyz,*qbm,*dt,*dtc,ek,*idimp,*nop,*npe,*nx,
                 *ny,*nz,*nxv,*nyv,*nzv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrbpush3lt_(float *part, float *fxyz, float *bxyz, float *qbm,
                     float *dt, float *dtc, float *ci, float *ek,
                     int *idimp, int *nop, int *npe, int *nx, int *ny,
                     int *nz, int *nxv, int *nyv, int *nzv, int *ipbc) {
   ckncgrbpush3lt(part,fxyz,bxyz,*qbm,*dt,*dtc,*ci,ek,*idimp,*nop,*npe,
                  *nx,*ny,*nz,*nxv,*nyv,*nzv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgpost3lt_(float *part, float *q, float *qm, int *nop, int *npe,
                   int *idimp, int *nxv, int *nyv, int *nzv) {
   ckncgpost3lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cknc2gpost3lt_(float *part, float *q, float *qm, int *nop,
                    int *npe, int *idimp, int *nxv, int *nyv,
                    int *nzv) {
   cknc2gpost3lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgjpost3lt_(float *part, float *cu, float *qm, float *dt,
                    int *nop, int *npe, int *idimp, int *nx, int *ny,
                    int *nz, int *nxv, int *nyv, int *nzv, int *ipbc) {
   ckncgjpost3lt(part,cu,*qm,*dt,*nop,*npe,*idimp,*nx,*ny,*nz,*nxv,*nyv,
                 *nzv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncgrjpost3lt_(float *part, float *cu, float *qm, float *dt,
                     float *ci, int *nop, int *npe, int *idimp, int *nx,
                     int *ny, int *nz, int *nxv, int *nyv, int *nzv,
                     int *ipbc) {
   ckncgrjpost3lt(part,cu,*qm,*dt,*ci,*nop,*npe,*idimp,*nx,*ny,*nz,*nxv,
                  *nyv,*nzv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void ckncdsortp3yzlt_(float *parta, float *partb, int *npic, int *idimp,
                      int *nop, int *npe, int *ny1, int *nyz1) {
   ckncdsortp3yzlt(parta,partb,npic,*idimp,*nop,*npe,*ny1,*nyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cknccguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   cknccguard3l(fxyz,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncacguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                    int *nye, int *nze) {
   ckncacguard3l(cu,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncaguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe,
                   int *nye, int *nze) {
   ckncaguard3l(q,*nx,*ny,*nz,*nxe,*nye,*nze);
   return;
}

/*--------------------------------------------------------------------*/
void ckncpois33_(float complex *q, float complex *fxyz, int *isign,
                 float complex *ffc, float *ax, float *ay, float *az,
                 float *affp, float *we, int *nx, int *ny, int *nz,
                 int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
                 int *nzhd) {
   ckncpois33(q,fxyz,*isign,ffc,*ax,*ay,*az,*affp,we,*nx,*ny,*nz,*nxvh,
              *nyv,*nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cknccuperp3_(float complex *cu, int *nx, int *ny, int *nz,
                  int *nxvh, int *nyv, int *nzv) {
   cknccuperp3(cu,*nx,*ny,*nz,*nxvh,*nyv,*nzv);
   return;
}

/*--------------------------------------------------------------------*/
void ckncibpois33_(float complex *cu, float complex *bxyz,
                   float complex *ffc, float *ci, float *wm, int *nx,
                   int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
                   int *nxhd, int *nyhd, int *nzhd) {
   ckncibpois33(cu,bxyz,ffc,*ci,wm,*nx,*ny,*nz,*nxvh,*nyv,*nzv,*nxhd,
                *nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncmaxwel3_(float complex *exyz, float complex *bxyz,
                  float complex *cu, float complex *ffc, float *ci,
                  float *dt, float *wf, float *wm, int *nx, int *ny,
                  int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                  int *nyhd, int *nzhd) {
   ckncmaxwel3(exyz,bxyz,cu,ffc,*ci,*dt,wf,wm,*nx,*ny,*nz,*nxvh,*nyv,
               *nzv,*nxhd,*nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rvx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rvx(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncemfield3_(float complex *fxyz, float complex *exyz,
                   float complex *ffc, int *isign, int *nx, int *ny,
                   int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                   int *nyhd, int *nzhd) {
   ckncemfield3(fxyz,exyz,ffc,*isign,*nx,*ny,*nz,*nxvh,*nyv,*nzv,*nxhd,
                *nyhd,*nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ckncwfft3rv3_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                   int *nxyzhd) {
   ckncwfft3rv3(f,*isign,mixup,sct,*indx,*indy,*indz,*nxhd,*nyd,*nzd,
                *nxhyzd,*nxyzhd);
   return;
}
