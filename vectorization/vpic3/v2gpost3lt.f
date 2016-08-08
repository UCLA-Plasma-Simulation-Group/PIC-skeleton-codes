c-----------------------------------------------------------------------
      subroutine V2GPOST3LT(part,q,qm,nop,npe,idimp,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c vectorizable version using guard cells
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(n,1) = position x of particle n
c part(n,2) = position y of particle n
c part(n,3) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c npe = first dimension of particle array
c idimp = size of phase space = 6
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
      implicit none
      integer nop, npe, idimp, nxv, nyv, nzv
      real qm
      real part, q
      dimension part(npe,idimp), q(nxv*nyv*nzv)
c local data
      integer npblk, lvect
      parameter(npblk=32,lvect=8)
      integer i, j, k, l, ipp, joff, nps, nn, mm, ll, nxyv
      real x, y, z, w, dx1, dxp, dyp, dzp, amx, amy, amz
c scratch arrays
      integer n, m
      real s
      dimension n(npblk), m(lvect), s(npblk,lvect)
      nxyv = nxv*nyv
      m(1) = 0
      m(2) = 1
      m(3) = nxv
      m(4) = nxv + 1
      m(5) = nxyv
      m(6) = nxyv + 1
      m(7) = nxyv + nxv
      m(8) = nxyv + nxv + 1
c loop over particles
      ipp = nop/npblk
c outer loop over number of full blocks
      do 40 k = 1, ipp
      joff = npblk*(k - 1)
c inner loop over particles in block
      do 10 j = 1, npblk
c find interpolation weights
      x = part(j+joff,1)
      y = part(j+joff,2)
      z = part(j+joff,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      n(j) = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
      s(j,1) = amx*amz
      s(j,2) = amy*amz
      s(j,3) = dyp*amz
      s(j,4) = dx1*amz
      s(j,5) = amx*dzp
      s(j,6) = amy*dzp
      s(j,7) = dyp*dzp
      s(j,8) = dx1*dzp
   10 continue
c deposit charge
      do 30 j = 1, npblk
!dir$ ivdep
      do 20 i = 1, lvect
      q(n(j)+m(i)) = q(n(j)+m(i)) + s(j,i)
   20 continue
   30 continue
   40 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 50 j = nps, nop
c find interpolation weights
      x = part(j,1)
      y = part(j,2)
      z = part(j,3)
      nn = x
      mm = y
      ll = z
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      dzp = z - real(ll)
      nn = nn + nxv*mm + nxyv*ll + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      dx1 = dxp*dyp
      dyp = amx*dyp
      amx = amx*amy
      amz = 1.0 - dzp
      amy = dxp*amy
c deposit charge
      x = q(nn) + amx*amz
      y = q(nn+1) + amy*amz
      z = q(nn+nxv) + dyp*amz
      w = q(nn+1+nxv) + dx1*amz
      q(nn) = x
      q(nn+1) = y
      q(nn+nxv) = z
      q(nn+1+nxv) = w
      mm = nn + nxyv
      x = q(mm) + amx*dzp
      y = q(mm+1) + amy*dzp
      z = q(mm+nxv) + dyp*dzp
      w = q(mm+1+nxv) + dx1*dzp
      q(mm) = x
      q(mm+1) = y
      q(mm+nxv) = z
      q(mm+1+nxv) = w
   50 continue
      return
      end
