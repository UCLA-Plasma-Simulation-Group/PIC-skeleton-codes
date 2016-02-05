import numpy as np
import math
import time
from numpy import vectorize
pi2 = 2*np.pi
float_type = np.float32
double_type = np.float64



#form factor for poisson equation
def form_factor(ffc,nx,ny,nxh,nyh,affp,ax,ay,kvector):
	dnx = pi2/float(nx)
	dny = pi2/float(ny)
	for i in xrange(nx):
		for j in xrange(ny):
			if(i > nxh):
				dkx = dnx * (nx-i)
			else:
				dkx = dnx* i

			if(j> nyh):
				dky=dny * (ny-j)
			else:
				dky = dny * j

			at = dky*dky+dkx*dkx
			gauss_fact = math.exp(-0.5*(dky*dky*ay*ay+dkx*dkx*ax*ax))
			#eliminate zero mode and nyquist mode (nxh,nyh)
			if(j == nyh):
				kvector.itemset((1,i,j),0)
			else:
				if(j>nyh):
					kvector.itemset((1,i,j),1.0j*dky)
				else:
					kvector.itemset((1,i,j),-1.0j*dky)

			if(i == nxh):
				kvector.itemset((0,i,j),0)
			else:
				if(i>nxh):
					kvector.itemset((0,i,j),1.0j*dkx)
				else:
					kvector.itemset((0,i,j),-1.0j*dkx)
			
			if(at==0):
				ffc.itemset((i,j),affp+1.0j)

			else:
				ffc.itemset((i,j),affp*gauss_fact/at+1.0j*gauss_fact)

				


#initialize particles on grid
def init_particles(part,vtx,vty,vx0,vy0,npx,npy,idimp,nptot,nx,ny,ipbc):
	edgelx =0.0
	edgely = 0.0
	at1 = np.float(nx)/np.float(npx)
	at2 = np.float(ny)/np.float(npy)
	# boundary conditions
	if(ipbc==2):
		edgelx = 1.0;
		edgely = 1.0
		at1 = np.float(nx-2)/npx
		at2 = np.float(ny-2)/npy
	elif(ipbc ==3):
		edgelx =1.0
		at1 = np.float(nx-2)/npx

	att = edgelx + at1 * (np.arange(npx)+0.5)
	for y in xrange(npy):
		at3 = edgely + at2 * (np.float(y)+0.5)
		part[0,npx*y:(npx*(y+1))] = np.copy(att)
		part[1,npx*y:(npx*(y+1))] = at3
		
	part[2,:]=vtx * np.random.normal(0,1,nptot)
	part[3,:]=vty * np.random.normal(0,1,nptot)
		
	#correct for drift
	vdrift_x = np.sum(part[2,:],dtype=double_type)/nptot -vx0
	vdrift_y = np.sum(part[3,:],dtype=double_type)/nptot - vy0

	np.subtract(part[2,:],vdrift_x,part[2,:])
	np.subtract(part[3,:],vdrift_y,part[3,:])
		


#deposit particles
def deposit2(part,qe,qme,nptot,idimp,nxe,nye):
	floor = np.int_(part)
	floor_x = floor[0,:]
	floor_y = floor[1,:]

	#calculate distances from nearest grid points
	dx = np.subtract(part[0,:],floor[0,:])
	dy = np.subtract(part[1,:],floor[1,:])
	adx = np.subtract(np.float(1.0),dx)
	ady = np.subtract(np.float(1.0),dy)

	#area weighting
	a1 = np.multiply(adx,ady)
	a2 = np.multiply(adx,dy)
	a3 = np.multiply(dx,ady)
	a4 = np.multiply(dx,dy)

	#python sucks at individual element manipulation, time to vectorize!
	ngrid = nxe*nye
	a1f = np.bincount(floor_y*nxe+floor_x,a1,ngrid)
	a2f = np.bincount((floor_y+1)*nxe+floor_x,a2,ngrid)
	a3f = np.bincount((floor_y)*nxe+floor_x+1,a3,ngrid)
	a4f = np.bincount((floor_y+1)*nxe+floor_x+1,a4,ngrid)
	qf = np.add(a1f,a2f)
	np.add(qf,a3f,qf)
	np.add(qf,a4f,qf)
	yy = np.arange(nye)
	xx = np.arange(nxe)
	for x in xrange(nxe):
		qe[x,yy] = qf[yy*nxe+x]
	np.multiply(qe,qme,qe)

 	
#guard cells for charge interpolation
def guard_cells(qe,nx,ny,nxe,nye):
	qe[0,:ny]+= qe[nx,:ny]
	qe[nx,:ny]=0.0
	qe[:nx,0] += qe[:nx,ny]
	qe[:nx,ny] = 0
	qe[0,0] += qe[nx,ny]
	qe[nx,ny]=0
	

#fourier analyze rho
def fourier_rho(qe_fft,qe,nx,ny):
	qe_fft[:nx,:ny] = np.fft.fft2(qe[:nx,:ny]/float(nx*ny))


#force in fourier space
def fourier_force(qe_fft,fxye,ffc,nx,ny,nxh,nyh,we,kvector):
	wp = 0.0
	#unwrap vectors
	kx = kvector[0,:,:]
	ky = kvector[1,:,:]
	at1 = np.multiply(ffc.real,ffc.imag,dtype = double_type)
	temp = np.multiply(qe_fft,at1)
	#print temp
	np.multiply(temp,kx,fxye[0,:,:])
	np.multiply(temp,ky,fxye[1,:,:])
	rho_squared = np.multiply(qe_fft,qe_fft.conj())
	#zero freq 
	rho_squared[0,0]=0
	#nyquist freq
	rho_squared[nxh,:] = 0
	rho_squared[:,nyh] = 0
	wp = np.sum(np.multiply(at1,rho_squared.real),dtype= double_type)
	#kinetic energy
	we[0] = .5*wp *float(nx*ny)

	






#force to real space
def real_force(fxyre,fxye,nx,ny,nxe,nye):
	fxyre[0,:nx,:ny] = float(nx*ny)*np.fft.ifft2(fxye[0,:nx,:ny]).real
	fxyre[1,:nx,:ny] = float(nx*ny)*np.fft.ifft2(fxye[1,:nx,:ny]).real


#guard cells for force
def guard_force(fxye,nx,ny):
	for k in xrange(ny):
		fxye.itemset((0,nx,k),fxye.item(0,0,k)) 
		fxye.itemset((1,nx,k),fxye.item(1,0,k)) 

	for j in xrange(nx):
		fxye.itemset((0,j,ny),fxye.item(0,j,0)) 
		fxye.itemset((1,j,ny),fxye.item(1,j,0)) 

	fxye.itemset((0,nx,ny),fxye.item(0,0,0)) 
	fxye.itemset((1,nx,ny),fxye.item(1,0,0)) 


#push these particles 
def push_particles(part,fxye,qbme,nx,ny,dt,wke,ipbc,nptot,ndim):
	qtm = qbme * dt
	edgelx = 0.0
	edgely = 0.0
	edgerx = float(nx)
	edgery = float(ny)
	if(ipbc == 2):
		edgelx = 1.0
		edgely = 1.0
		edgerx = float(nx-1)
		edgery = float(ny-1)
	elif(ipbc == 3):
		edgelx = 1.0
		edgerx = float(nx-1)
	
	floor = part.astype(int)
	floor_x = floor[0,:]
	floor_y = floor[1,:]

	#calculate distances from lower grid points
	dx = np.subtract(part[0,:],floor[0,:])
	dy = np.subtract(part[1,:],floor[1,:])
	adx = np.subtract(np.float(1.0),dx)
	ady = np.subtract(np.float(1.0),dy)

	#area weighting
	a1 = np.multiply(adx,ady)
	a2 = np.multiply(adx,dy)
	a3 = np.multiply(dx,ady)
	a4 = np.multiply(dx,dy)

	da = np.empty((ndim,nptot))
	temp = np.empty((ndim,nptot))
	
	#weight the forces
	np.multiply(fxye[:,floor_x,floor_y],a1,da)
	np.multiply(fxye[:,floor_x,floor_y+1],a2,temp)
	np.add(da,temp,da)

	np.multiply(fxye[:,floor_x+1,floor_y],a3,temp)
	np.add(da,temp,da)
	np.multiply(fxye[:,floor_x+1,floor_y+1],a4,temp)
	np.add(da,temp,da)
	#calculate velocity
	vel = np.copy(part[2:4,:])
	v_new = np.add(vel, qtm*da)
	np.copyto(part[2:4,:],v_new)
	v_sum = np.add(v_new,vel)
	sum1 = np.sum(np.multiply(v_sum,v_sum))
	
	#calculate position
	pos_temp = np.add(part[0:2,:],v_new * dt)
	
	#boundary conditions
	x_temp = pos_temp[0,:]
	y_temp = pos_temp[1,:]
	
	#periodic
	if(ipbc == 1):
		x_temp = np.where(x_temp<edgelx,x_temp+edgerx,x_temp) 
		x_temp = np.where(x_temp>=edgerx,x_temp-edgerx,x_temp) 
		y_temp = np.where(y_temp<edgely,y_temp+edgery,y_temp) 
		y_temp = np.where(y_temp>=edgery,y_temp-edgery,y_temp) 
	#reflecting
	elif(ipbc == 2):
		indices = np.where(np.logical_or(x_temp<edgelx, x_temp>= edgerx))
		x_temp[indices] = part[0,indices]
		part[2,indices] = -part[2,indices] 
		indices = np.where(np.logical_or(y_temp<edgely, y_temp>= edgery))
		y_temp[indices] = part[1,indices]
		part[3,indices] = -part[3,indices] 
	else: 
		print 'ipbc must be periodic or reflecting'
	np.copyto(part[0,:],x_temp)
	np.copyto(part[1,:],y_temp)
	wke[0] = sum1*0.125
	








#auxiliary methods for debugging purposes

# #deposit particles
# def deposit(part,qe,qme,nptot,idimp,nxe,nye):
# 	t = time.time()
# 	f = 0
# 	f2 = 0
# 	print 'hi'
	
# 	for j in xrange(nptot):
# 		if(j% 100000 == 0):
# 			print j
# 			print f
# 			print f2 
# 		t2 = time.time()
# 		#xxyy = part[:,j]
# 		xx = part[0,j]
# 		yy = part[1,j]
# 		x_floor = int(xx)
# 		y_floor = int(yy)
# 		dx = xx-x_floor
# 		dy = yy-y_floor
# 		f+= time.time()-t2
# 		adx = 1.0-dx
# 		ady = 1.0-dy
# 		#mm = nxe*y_floor
# 		#mp = mm +nxe
# 		#np = x_floor+1
# 		t2 = time.time()
# 		qe.itemset((x_floor,y_floor),qe.item(x_floor,y_floor)+qme*adx*ady)
# 		qe.itemset((x_floor,y_floor+1),qe.item(x_floor,y_floor+1)+qme*adx*dy)
# 		qe.itemset((x_floor+1,y_floor),qe.item(x_floor+1,y_floor)+qme*dx*ady)
# 		qe.itemset((x_floor+1,y_floor+1),qe.item(x_floor+1,y_floor+1)+qme*dx*dy)
		
#  		f2+= time.time()-t2
#  		#f += time.time()-t2
#  	print time.time()-t


	#weight the forces
	# np.multiply(force_re[0,flo_x,flo_y],a1,da[0,:])
	# np.multiply(force_re[1,flo_x,flo_y],a1,da[1,:])
	# np.multiply(force_re[0,flo_x,flo_y+1],a2,temp[0,:])
	# np.multiply(force_re[1,flo_x,flo_y+1],a2,temp[1,:])
	# np.add(da,temp,da)
	# np.multiply(force_re[0,flo_x+1,flo_y],a3,temp[0,:])
	# np.multiply(force_re[1,flo_x+1,flo_y],a3,temp[1,:])
	# np.add(da,temp,da)
	# np.multiply(force_re[0,flo_x+1,flo_y+1],a4,temp[0,:])
	# np.multiply(force_re[1,flo_x+1,flo_y+1],a4,temp[1,:])
	# np.add(da,temp,da)