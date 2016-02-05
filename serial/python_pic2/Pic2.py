#-------------------
#Skeleton 2D Electrostatic PIC code
# written by Thamine Dalichaouch, Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import numpy as np
import math
from push import *
import time
import matplotlib.pyplot as plt
int_type = np.int32
double_type = np.float64
float_type = np.float32
complex_type = np.complex64

def main():
	#----- Same format as other Skeleton Codes written by Viktor Decyk -----

	#number of cells in each direction-- nx = 2^indx
	indx = 8; indy=8
	#indx = 3; indy=3
	#number of particles in each direction
	#npx = 2024; npy = 2024
	npx = 500;npy = 500
	#ndim = number dimensions 
	ndim =2
	#tend = time at end of plasma simulation in units of plasma freq.
	#dt = time interval
	#qme = electron charge in units of |e|
	tend =500.0; dt = 0.1; qme = -1
	#vtx/vty are thermal velocities
	#vx0/vy0 are drift velocities
	vtx = 1.0; vty = 1.0; vx0 = 0; vy0 = 0;
	# ax/ay = smoothed particle size in x/y direction
	ax = .912871; ay = .912871
	# idimp = number of particle coordinates = 4
	# ipbc = particle boundary condition: 1 = periodic, 2 = reflecting
	# sortime = number of time steps between standard electron sorting
	idimp = 4; ipbc = 2; sortime = 50
	# wke/we/wt = particle kinetic/electric field/total energy
	wke = np.zeros((1),float_type)
	we = np.zeros((1),float_type)
	wt = np.zeros((1),float_type)

	# declare and initialize timing data
	itime = np.empty((4),int_type)
	tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
	tpush = 0.0; tsort = 0.0
	dtime = np.empty((1),double_type)

	# initialize scalars for standard code
	# nptot = total number of particles in simulation
	# nx/ny = number of grid points in x/y direction
	nptot = npx*npy; nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))

	# nyquist frequencies
	nxh = int(nx/2); nyh = max(1,int(ny/2))

	#guard cell definitions; do not change guard cells
	nxe = nx + 2; nye = ny + 1; 


	# nloop = number of time steps in simulation
	# ntime = current time step
	nloop = int(tend/dt +0.001); ntime = 0

	phi = np.empty((nloop,nxh,nyh),complex_type)
	phi.fill(0.0)
	#extra params
	qbme = qme 
	
	# inverse particles per cell normalizing factor
	affp = float(nx*ny)/nptot 

	# allocate data for code
	# part, part2 = particle arrays
	part = np.empty((idimp,nptot),float_type,'F')
	if(sortime>0):
		part2 = np.empty((idimp,nptot),float_type,'F')

	# electron charge density with guard cells
	qe = np.empty((nxe,nye),float_type,'F')
	qe2 = np.empty((nxe,nye),float_type,'F')

	# wavenumbers (ndim,kx,ky) for fourier analysis
	kvector = np.empty((ndim,nxe,nye),complex_type,'F')
	kvector.fill(0.0)

	#fft electron charge density
	qe_fft = np.empty((nxe,nye),complex_type,'F')

	# smoothed electric field with guard cells in fourier space
	fxye = np.empty((ndim,nxe,nye),complex_type,'F')

	# smoothed electric field with guard cells in real space
	fxyre = np.empty((ndim,nxe,nye),float_type,'F')
	
	# form factor array for poisson solver
	ffc = np.empty((nxe,nye),complex_type,'F')
	ffc.fill(0.0)

	#form factor array do this once
	form_factor(ffc,nx,ny,nxh,nyh,affp,ax,ay,kvector)

	#initialize particle velocities and 
	init_particles(part,vtx,vty,vx0,vy0,npx,npy,idimp,nptot,nx,ny,ipbc)
	
	#main loop
	for ntime in xrange(0,nloop):
		#zero matrices
		fxye.fill(0.0)
		qe.fill(0.0)
		fxyre.fill(0.0)
		qe_fft.fill(0.0)

		#deposit charge
		t = time.time()
		deposit2(part,qe,qme,nptot,idimp,nxe,nye)
		tdpost += time.time() - t
			
		
		#add guard cells
		t = time.time()
		guard_cells(qe,nx,ny,nxe,nye)
		tguard += time.time() - t

		
		#charge to fourier space!
		t = time.time()
		fourier_rho(qe_fft,qe,nx,ny)
		tfft += time.time() - t
		

		#calculate field in fourier space
		t = time.time()
		fourier_force(qe_fft,fxye,ffc,nx,ny,nxh,nyh,we,kvector)
		tfield += time.time()-t
		np.copyto(phi[ntime,:,:],fxye[0,:nxh,:nyh])	

		#force to real space, updates fxye
		t = time.time()
		real_force(fxyre,fxye,nx,ny,nxe,nye)
		tfft += time.time() -t
		
		t = time.time()
		guard_force(fxyre,nx,ny)
		tguard += time.time()-t


		#push particles updates part,wke
		wke[0] = 0.0
		t = time.time()
		push_particles(part,fxyre,qbme,nx,ny,dt,wke,ipbc,nptot,ndim)
		tpush += time.time() -t
		if(sortime > 0):
			if(ntime%sortime == 0):
				t = time.time()
				#sort by y axis
				part = part[:,np.argsort(part[1,:])]
				tsort = tsort + time.time() -t

		if (ntime<1):
			print "Initial Field, Kinetic and Total Energies:"
			print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)
      		
      


	print "ntime = ", ntime
	print "Final Field, Kinetic and Total Energies:"
	print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)

	print ""
	print "deposit time = ", tdpost
	print "guard time = ", tguard
	print "solver time = ", tfield
	print "fft time = ", tfft
	print "push time = ", tpush
	tfield = tfield + tguard + tfft
	print "sort time = ", tsort
	print "total solver time = ", tfield
	par_time = tdpost + tpush + tsort
	print "total particle time = ", par_time
	wt = par_time + tfield
	print "total time = ", wt
	print ""

	#diagnostics for code to view a plot of the disperion relation for Ex
	# default set to false
	print_diagnostics = True

	if(print_diagnostics):
		plt.figure(1,figsize=(12,6))
		freq_start = 0.7
		freq_end = 3
		fx_st = int(freq_start*dt*nloop/(2*np.pi))
		fx_end = int(freq_end*dt*nloop/(2*np.pi))
		print fx_st,fx_end
		freq_start = (2*np.pi/(dt*nloop))*fx_st
		freq_end = (2*np.pi/(dt*nloop))*fx_end
		plt.imshow(np.abs(np.fft.fft(phi[:,:,1],axis=0)[fx_st:fx_end,:]),origin ='lower',extent=[0,np.pi,freq_start,freq_end])
		#plt.axis([0,np.pi, freq_start,freq_end])
		plt.colorbar()
		plt.savefig('disp',dpi = 400,bbox_inches = 'tight')

	#plt.show()

if __name__ == "__main__":
	main()