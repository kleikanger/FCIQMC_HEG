import sys
import os
import math
import time
from datetime import datetime, date, time

class runcode:

	'Constructor'
	def __init__(self):

		'''
		Class variables, default values.
		'''
		self.d_dt = 0.001
		self.dS = 5.0
		self.d_rs = 0.1
		self.d_mu = 0.0
		self.i_numpart = 10
		self.i_numorbs = 18
		self.s_system = '2DHEG' # '3DHEG'
		self.i_numshells_basis = -1 #if negative, k_c sets the basis
		self.i_numclosedshells = -1 #if negative, k_c sets the basis
		self.d_kc=1; #N 10 //2 shells
		self.b_writeodata=True
		self.i_numthreads=2
		self.i_initiatorlimit=3
		self.i_limit_nw=20000
		self.i_num_loops=1000000
		self.i_n_start_collecting_e=5000
		self.timestamp = self.set_time_now()
		self.i_update_s_frequency=1;
		self.s_ofilename=self.timestamp;
		self.d_pexone=0.0
		self.d_xi=0.3 #see:umrigars update formula. S(n+1) = S(n)-(d_xi/d_dt)*(N_W(n)/N_W(n+1)) 
		self.d_detweight=1.3
		self.i_maxndets=1000000
		self.i_ranseed=-1 #RNG seed, (automatically set to timestamp if -1 of <0).
		self.d_redistlimit=0.04
		self.BITSET_BOUNDSCHECK = True

	def runcode(self):
		if self.i_numshells_basis>0:
			self.loadNumShellsBasis()
		if self.i_numclosedshells>0:
			self.loadNumClosedShells()
		self.set_timestep()
		self.writeIniFiles()
		self.makerun()

	def loadNumShellsBasis(self):
		n = self.i_numshells_basis

		if self.s_system is '2DHEG':
			import shells2DHEG
			self.i_numorbs = shells2DHEG.states[n-1]
			self.d_kc = shells2DHEG.r_kc[n-1]+.001
		elif self.s_system is '3DHEG':
			import shells3DHEG
			self.i_numorbs = shells3DHEG.states[n-1]
			self.d_kc = shells3DHEG.r_kc[n-1]+.001

		else: 
			print 'err: no system set'
	
	def loadNumClosedShells(self):
		n = self.i_numclosedshells

		if self.s_system is '2DHEG':
			import shells2DHEG
			self.i_numpart = shells2DHEG.states[n-1]

		elif self.s_system is '3DHEG':
			import shells3DHEG
			self.i_numpart = shells3DHEG.states[n-1]

		else: 
			print 'err: no system set'

	#generate string with system time
	def set_time_now(self):
		return '%s'%datetime.now().strftime("%Y_%b_%d_%H%M%S")

	def compile(self):
		os.system('make --silent')

	def makerun(self):
		os.system('make run --silent')

	def set_timestep(self):
		k = 0.1 # k in (0,1). 1 represents the largest "sensible" timestep (see: Shepherd)
		pi = math.pi
		N = self.i_numpart
		M = self.i_numorbs
		if not N is M:
			if self.s_system is '3DHEG':
				l = (4./3.*N*pi*self.d_rs**3.)**(1./3.)
				self.d_dt = k * 4*pi*l / float(N*(N-1)*(M-N))
			elif self.s_system is '2DHEG':
				l = ( pi*N*self.d_rs**2. )**(1./2.)
				self.d_dt = k * 4*l / float(N*(N-1)*(M-N))

	def writeIniFiles(self):
		
		fsout = open('last_run_definitions.h', 'w')
		fsout.write('#ifndef DEFINITIONS_H\n')
		fsout.write('#define DEFINITIONS_H\n')
		fsout.write('   #define INUM_ORBITALS (int)%s\n'%self.i_numorbs)
		fsout.write('   #define BITSET_BOUNDSCHECK %i\n'%self.BITSET_BOUNDSCHECK)
		fsout.write('#endif\n')
		fsout.close()
		
		os.system('cp last_run_definitions.h inifiles/%s.h'%self.s_ofilename)

		fsout = open('imputVars.cpp', 'w')
		fsout.write('#include "imputVars.h"\n')
		fsout.write('void imputVars::init()\n')
		fsout.write('{\n')
		fsout.write(' 	d_dt=%f;\n'%self.d_dt)
		fsout.write(' 	dS=%f;\n'%self.dS)
		fsout.write(' 	d_rs=%f;\n'%self.d_rs)
		fsout.write(' 	d_mu=%f;\n'%self.d_mu)
		fsout.write(' 	i_numpart=%i;\n'%self.i_numpart)
		fsout.write(' 	s_system="%s";\n'%self.s_system)
		fsout.write(' 	d_kc=%f;\n'%self.d_kc)
		fsout.write(' 	b_writeodata=%i;\n'%self.b_writeodata)
		fsout.write(' 	i_numthreads=%i;\n'%self.i_numthreads)
		fsout.write(' 	i_initiatorlimit=%i;\n'%self.i_initiatorlimit)
		fsout.write(' 	i_limit_nw=%i;\n'%self.i_limit_nw)
		fsout.write(' 	i_num_loops=%i;\n'%self.i_num_loops)
		fsout.write(' 	i_n_start_collecting_e=%i;\n'%self.i_n_start_collecting_e)
		fsout.write(' 	i_update_s_frequency=%i;\n'%self.i_update_s_frequency)
		fsout.write(' 	s_ofpath="%s";\n'%('outputdata/'+self.s_ofilename))
		fsout.write(' 	d_pexone=%f;\n'%self.d_pexone)
		fsout.write(' 	d_xi=%f;\n'%self.d_xi)
		fsout.write(' 	d_detweight=%f;\n'%self.d_detweight)
		fsout.write(' 	i_maxndets=%i;\n'%self.i_maxndets)
		fsout.write(' 	i_ranseed=%i;\n'%self.i_ranseed)
		fsout.write(' 	d_redistlimit=%f;\n'%self.d_redistlimit)
		fsout.write('}');
		fsout.close()
		
		if self.b_writeodata:
			os.system('cp imputVars.cpp inifiles/%s.cpp'%self.s_ofilename)



	
		
