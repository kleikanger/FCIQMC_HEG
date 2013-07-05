'''
Example:
	import pyDataAnal
	#construct object with the timestamp 2012_Dec_21_182455
	a = pyDataAnal.pyDataAnal('2012_Dec_21_182455', 'd 3i d', 5)
	#load data
	a.readBF2array()
	#blockingdata start 60000, minBS 200. maxbs 1400. blockstep 20 
	#PE = projected energy
	a.blcAnalysis(60000, 200, 1400, 20, 'PE')
	#GE = generational enstomator
	a.blcAnalysis(60000, 200, 1400, 20, 'GE')
	#
	a.getEnergyAverages(60000, 0)
	#
	a.plotPopulation()
	#
	a.plotTrailingE(60000)

	a_data[:,0] = dS ( shifts )
	a_data[:,1] = num. walkers.
	a_data[:,2] = num. dets.
	a_data[:,3] = n_0 (determinants at reference)
	a_data[:,4] = e_p (projected energy)
'''

import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import sys
import csv

class pyDataAnal:

	'Constructor'
	def __init__(self, s_datetag_, s_structtype_, i_nElem_):
		'Class variables'
		self.s_datetag = s_datetag_
		self.s_structtype = s_structtype_
		self.i_nElem = i_nElem_
		'''
		Read data from the inifile
		'''
		'''
		s = []
		with open('inifiles/'+s_datetag_+'.ini') as f:
			c = csv.reader(f, delimiter='=', skipinitialspace=True)
			for line in c:
				s.append(line)

		self.i_numpart = np.int(s[1][1])
		self.d_dt = np.double(s[2][1])
		self.d_initialS = np.double(s[3][1])
		self.i_limnw = np.int(s[4][1])
		self.i_numloops = np.int(s[5][1])
		self.i_nstartcolle = np.int(s[6][1])
		self.i_update_s_frequency = np.int(s[7][1])
		self.b_energycut = bool(s[8][1])
		self.b_useveff = bool(s[9][1])
		self.d_lambd = np.double(s[10][1])
		self.s_initialstate = s[11][1]
		self.i_initiatorlimit = s[12][1]
		'''
		'''
		Read data from the definitions file
		'''
		'''
		s = []
		with open('inifiles/def'+s_datetag_+'.h') as f:
			c = csv.reader(f, delimiter=')', skipinitialspace=True)
			for line in c:
				s.append(line)
		self.i_R = int(s[2][1])
		'''
		self.d_omega = 1
		self.a_data = np.array([])
		self.a_PEdata = np.array([])

		self.la_PEdata = []
		self.la_PEdatan0 = []

		#self.a_bdataPE = []
		#self.a_bdataGE = []

		#self.i_m
		#self.i_s

	'''
	Read data from file
	'''
	def readBF2array(self, units='W'):

		structT = self.s_structtype
		filename = 'outputdata/' + self.s_datetag + 'mixed.dat'
		nElem = self.i_nElem

		with open(filename, "rb") as f:
			structLen = struct.calcsize('='+structT)
			numelem = os.fstat(f.fileno()).st_size/structLen
			rarray = np.zeros([numelem, nElem])
			#rarray = []        

			i=0
			byte = f.read(structLen)
			while byte:
				#remove this line is the length of the file is correct
				#slows down the code!
				if ( len(byte)==structLen ):
					rarray[i] = np.double(struct.unpack('='+structT, byte))
					#rarray.append(s)
					byte = f.read(structLen)
				else:
					break

				i += 1
			f.close()

			if (units == 'W'):
				rarray[:,0] *= self.d_omega
				rarray[:,4] *= self.d_omega

			self.a_data = rarray

	'''
	Read data from array with doubles. Append array to list self.la_PEdata
	'''
	def readPE2array(self, filename, units='W', ret=0):
		
		with open(filename, "rb") as f:
			
			structLen = struct.calcsize('=d')
			numelem = os.fstat(f.fileno()).st_size/structLen
			rarray = np.zeros([numelem], dtype=np.float64)

			i = 0
			byte = f.read(structLen)
			while byte:
				rarray[i] = np.float64(struct.unpack('=d', byte))
				byte = f.read(structLen)
				i += 1
			
			f.close()

		if (units == 'W'):
			rarray[:] *= self.d_omega
		if ret:
			return rarray
		else:
			self.la_PEdata.append(rarray)
	'''
	Read data from array with doubles. Append array to list self.la_PEdata
	'''
	def readN02array(self, filename, units='W', ret=0):
		
		with open(filename, "rb") as f:
			
			structLen = struct.calcsize('=i')
			numelem = os.fstat(f.fileno()).st_size/structLen
			rarray = np.zeros([numelem], dtype=np.int32) #XXX int64 some os

			i = 0
			byte = f.read(structLen)
			while byte:
				rarray[i] = np.int32(struct.unpack('=i', byte))
				byte = f.read(structLen)
				i += 1
			
			f.close()

		if (units == 'W'):
			rarray[:] *= self.d_omega
		if ret:
			return rarray
		else:
			self.la_PEdatan0.append(rarray)
	'''
	process all files starting with the imput string ifname
	on the form ifname + 'r%i'%i + .dat, where i is an integer 
	between 0 and 100
	'''
	def readAllPE(self, ifname = '', maxnumprocs = 100):
		if not ifname:
			ifname = self.s_datetag
		print('processing all datafiles %s*'%ifname)
		s=os.listdir('outputdata/') #filepath = outputfiles
		#find names starting with imp[:-1]>
		for nam in s:
			for i in range (0, maxnumprocs):
				if nam == ifname + 'er%i'%i + '.dat':
					print 'reading %s'%nam
					self.readPE2array('outputdata/'+nam)
		for nam in s:
			for i in range (0, maxnumprocs):
				if nam == ifname + 'nr%i'%i + '.dat':
					print 'reading %s'%nam
					self.readN02array('outputdata/'+nam)

	'''
	Get the blocks from start to the end of the array
	Store remove all elem before iteration start.
	Store the elem s.t. all elements from the same cycle are
	stores after each other. This is necc. to get the correct 
	blocking result.

	data = either la_PEdatan0 or la_PEdata
	'''
	def getPERarr(self, start, data):
		'''
		n = 0
		counter = np.zeros(len(self.la_PEdata))
		res = []
		for elem in self.a_PEdata:
			if start>0:
				if elem<-1e11:
					start -= 1
			else:
				if elem>-1e11:
					res.append(elem)
				else:
					n +=1
		return np.array(res)
		'''

		res = []
		num_arrays = len(data)

		sep = data[0][0]

		i = 0
		
		a_counter = np.zeros(num_arrays)
		a_start = np.ones(num_arrays)*start
		#tot nr of elem in the different data arrays
		a_totlelem = []
		for elem in data:
			a_totlelem.append(len(elem))
		a_totlelem = np.array(a_totlelem)

		while (a_counter-a_totlelem).any(): #stops when a_counter = a_totlelem
			for elem in data[i][a_counter[i]::]:
				a_counter[i] += 1
				if a_start[i]>0:
					if elem==sep: #large negative vals as separators
						a_start[i] -= 1
				else:
					if elem>sep:
						res.append(elem)
					else:
						i = (i+1)%num_arrays #cyclic
						break #only breaks the last forloop
					
		return np.array(res)

	'''
	Find the statistical error using blocking
	0ne vec (sum<i|h|j>/n0, ....)
	'''
	def blcsigma(self, data, blcsize):

		l = data.size
		s = 0
		d = 0
		#nblocks = l/np.double(blcsize)
		nblocks = l/blcsize #integer division

		for i in range(0, nblocks * blcsize, blcsize):
			tmp = np.sum(data[i:i+blcsize])
			s += tmp
			d += tmp**2

		d = d / (np.double(blcsize**2 * nblocks))
		s = s / (np.double(blcsize * nblocks))
		#print s
		#print d
		#print l
		#print data.size
		return np.sqrt( abs(s**2-d) / nblocks )

	'''
	Find the statistical error using blocking
	Two vecs (sum<i|h|j>, ....) and (n0, ....)
	'''
	def blcsigmaPER(self, e, n0, blcsize):

		l = e.size
		s = 0
		d = 0

		nblocks = l/blcsize #integer division
		
		for i in range(0, nblocks * blcsize, blcsize):
			tmp = np.sum(e[i:i+blcsize], dtype=np.float64)
			tmp /= np.sum(n0[i:i+blcsize], dtype=np.float64) 
			s += tmp
			d += tmp**2

		d = d / (np.double(nblocks))
		s = s / (np.double(nblocks))

		err = np.sqrt( abs(s**2-d) / nblocks ) 
		print 'mean:%e, err:%e, blocksize:%i'%(s,err,blcsize)

		return err

	'''
	Find the statistical error using blocking
	Callable method
	'''
	def blcAnalysis(self, start, s_estimator, minblcsize=0, maxblcsize=0, blcstep=0):

		blcresult = []

		#projected estimator exact
		if (s_estimator == 'PER'):
			e = self.getPERarr(start, self.la_PEdata)
			n0 = self.getPERarr(start, self.la_PEdatan0)

			if minblcsize is 0:
				minblcsize = 8
				maxblcsize = e.size / 10
				blcstep = ( maxblcsize - minblcsize ) /100

			#for blcsize in range(minblcsize, maxblcsize, blcstep):
			#	blcresult.append(self.blcsigmaPER(e, n0, blcsize))
			blcsize = minblcsize
			while blcsize < maxblcsize:
				blcresult.append(self.blcsigmaPER(e, n0, blcsize))
				blcsize *= 2

			#projected energy approx
		elif (s_estimator == 'PE'):
			
			e = self.a_data[start::, 4] * self.a_data[start::, 2]
			n0 = self.a_data[start::, 2]
			
			if minblcsize is 0:
				minblcsize = 8
				maxblcsize = e.size / 10
				blcstep = max(( maxblcsize - minblcsize ) /100, 1)

			#for blcsize in range(minblcsize, maxblcsize, blcstep):
			#	blcresult.append(self.blcsigmaPER(e, n0, blcsize))
			blcsize = minblcsize
			while blcsize < maxblcsize:
				blcresult.append(self.blcsigmaPER(e, n0, blcsize))
				#blcresult.append(self.blcsigma(self.a_data[start::, 4] , blcsize))
				blcsize *= 2


			#general estimator
		elif (s_estimator == 'GE'):
			data = self.a_data[start::,0]

			if minblcsize is 0:
				minblcsize = 8
				maxblcsize = data.size / 10
				blcstep = max(( maxblcsize - minblcsize ) /100, 1)

			#for blcsize in range(minblcsize, maxblcsize, blcstep):
			#	blcresult.append( self.blcsigma(data, blcsize) )
			blcsize = minblcsize
			while blcsize < maxblcsize:
				blcresult.append(self.blcsigma(data, blcsize))
				blcsize *= 2

#		elif (s_estimator == 'PEPOE'):
#			if minblcsize is 0:
#				minblcsize = 16
#				maxblcsize = data.size / 50
#				blcstep = max(( maxblcsize - minblcsize ) /100, 1)
#			
#			E = self.a_data[200000:,4] * self.a_data[200000:,3] / len(E)
#			N = self.a_data[200000:,3] / len(N)
#
#			mE = E.mean()
#			mN = N.mean()
#
#			blcE = []
#			blcN = []
#		
#			for blcsize in range(minblcsize, maxblcsize, blcstep):
#				blcE.append( self.blcsigma(E, blcsize) )
#			for blcsize in range(minblcsize, maxblcsize, blcstep):
#				blcN.append( self.blcsigma(N, blcsize) )
#
#			sE = 
#			plt.figure()
#			plt.plot(blcN)
#			sN = input('input <N> err: ')
#			err_pe = input('input <E_P> err: ')
#
#
#
#			return sE**2 / mE**2 + sN**2 / mN**2 - 2 * mE / mN * np.cov(E,N)[0,1]

		return np.array(blcresult)

	'''
	Do blocking analysis
	Plot data
	'''
	def blcPlots(self, start, minblcsize=0, maxblcsize=0, blcstep=0):
	
		res_pe = self.blcAnalysis(start, 'PE', minblcsize, maxblcsize, blcstep)
		res_ge = self.blcAnalysis(start, 'GE', minblcsize, maxblcsize, blcstep)

		#x=np.arange(len(res_pe))*blcstep
		x = []
		blcsize = 8
		maxblcsize = len(self.a_data[start:,0]) / 10
		while blcsize < maxblcsize:
			x.append(blcsize)
			blcsize *= 2
		x = np.array(x)

		nn = np.float(len(self.a_data[start:,0])) / x

		errge = res_ge / ( np.sqrt(2*(nn-1)) )
		errpe = res_pe / ( np.sqrt(2*(nn-1)) )

		plt.figure()
		plt.xscale('log')
		#plt.plot(x, res_ge, 'g-s')
		#plt.plot(x, res_pe, 'r-d')
		plt.errorbar(x,res_ge,yerr=errge,fmt='g-d'); plt.show(block=False)
		plt.errorbar(x,res_pe,yerr=errpe,fmt='r-.s'); plt.show(block=False)
		plt.legend(["$\langle S \\rangle$", "$\langle E_p \\rangle$"])#, loc=0)
		plt.ylabel('$\epsilon$', size=17)
		plt.xlabel('$N_{block}$', size=17)
		plt.show(block=False)

	'''
	Plot the trailing averages of the PE and the GE after element nr <start>
	plot the <exact> value as a line,
	'''
	def plotTrailingE(self, start, exact=999999):

		data = self.a_data

		#calculate trailing averages
		avebefore_ge = []
		n = data[start:,0].size/1000
		for a in range(start+n, data[:,0].size, n):
			avebefore_ge.append(data[start:a, 0].sum()/data[start:a, 0].size)
			#print( '%d %i', (data[a::,0].sum()/data[a::,0].size, a))
		avebefore_ge = np.array(avebefore_ge)

		avebefore_pe = []
		n = data[start:,4].size/1000
		for a in range(start+n, data[:,4].size, n):
			#samples does not have the same size
			epe = ( data[start:a, 4] * data[start:a, 3] ).sum()
			num = data[start:a, 3].sum()
			avebefore_pe.append(epe / num)

		avebefore_pe = np.array(avebefore_pe)

		x=np.array(range(0, avebefore_ge.size, 1))*n+start

		plt.figure()

		#plt.plot(x, data[start::data[start:,0].size/1000, 0])
		if (exact != 999999):
			plt.plot(x, np.ones(avebefore_ge.size)*exact, '-.')#21.12992)
		plt.plot(x, avebefore_ge, 'g-d', markevery = 100)
		plt.plot(x, avebefore_pe, 'r-s', markevery = 100)

		if (exact != 999999):
			plt.legend(("$exact=%f$"%exact,\
					"$\langle S \\rangle$",\
					"$\langle E_p \\rangle$"),\
					loc=0)
		else:
			plt.legend(("$\langle S \\rangle$",\
					"$\langle E_p \\rangle$"),\
					loc=0)
		plt.xlabel('Number of iterations')# ($1/\\tau$)', size=17)
		#plt.ylabel('$$',size=17)
		#plt.title(r'$\,N_P=2,\,m=0,\,s=0,\,\omega=1,\,R=7$', size=17)
		#plt.ylim(21.08 ,21.18)
		plt.show(block=False)
	'''
	Plot the fluctuations in E_p(t) and S(t)
	'''
	def plotInstantE(self):
		data = self.a_data

		x=range(len(data[0::,0])/1000)
		x = np.array(x)*1000
		plt.xlabel('$n=$ Number of iterations')# ($1/\\tau$)', size=17)
		plt.plot(x, data[1::1000,0],'g-d', markevery = 80)
		plt.plot(x, data[1::1000,4],'r-s', markevery = 80)
		plt.legend(['$S(n\cdot\\tau)$','$E_p(n\cdot\\tau)$'])
		plt.show(block=False)
	'''
	Plot n spawned and n killed as well?
	'''
	def plotPopulation(self):
		data = self.a_data

		x = range(len(data[0::,0])/1000)
		x = np.array(x)*1000

		plt.figure()
		plt.xlabel('iterations')
		plt.plot(x, data[1001::1000,1],'r-d', markevery = 80)
		plt.plot(x, data[1001::1000,3],'g-o', markevery = 80)
		plt.plot(x, data[1001::1000,2],'k-s', markevery = 80)
		plt.legend(['$N_W$','$N_0$','$N_{determinants}$'])
		#ax.grid()
		plt.show(block=False)
	'''
	Do all
	'''
	'''
	plotAll(self, start):
		self.getEnergyAverages(start, 0)
		self.blcPlots()
		self.plotTrailingE()
		self.plotInstantE()
		self.plotPopulation()
	'''

	'''
	Find the average energy after iteration nr <start>
	'''
	def getEnergyAverages(self, start, ret):
		data = self.a_data
		length = data[start::, 0].size

		ege = data[start::, 0].sum()/length

		epe = ( data[start::, 4] * data[start::, 3] ).sum()
		num = data[start::, 3].sum()

		epe = epe / np.double(num)

		if (ret == 0):
			#generational estimator
			print( '\nStatistical mean: generational estimator after it. nr. %i: %f'\
					%(start, ege) )
			#projected energy
			print( 'Statistical mean: projected energy after it. nr. %i: %f\n'\
					%(start, epe) )
		else:
			return ege, epe

def pb(fname, start):
	a = pyDataAnal(fname, 'd q q q  d', 5)
	a.readBF2array()
	a.blcPlots(start);
	err_ge = input('input <S> err: ')
	err_pe = input('input <E_P> err: ')
	plt.close()
	e_ge, e_pe = a.getEnergyAverages(start, 1); 

	return e_ge, e_pe, err_ge, err_pe

def listRes(fnames, start, ret=False):

	if type(fnames)==str:
		fnames = [fnames]

	e_ge = 0
	e_pe = 0
	err_ge = 0
	err_pe = 0

	l_e_ge = []
	l_e_pe = []
	l_err_ge = []
	l_err_pe = []

	for nam in fnames:
		e_ge, e_pe, err_ge, err_pe = pb(nam, start)
		l_e_ge.append(e_ge)
		l_e_pe.append(e_pe)
		l_err_ge.append(err_ge)
		l_err_pe.append(err_pe)

	print ('\n results: \n')
	for i in range(0, len(l_err_ge)):
		print(fnames[i])
		print 'GE=%f+-%f, PE=%f+-%f'%(l_e_ge[i], l_err_ge[i], l_e_pe[i], l_err_pe[i])

	if ret:
		return l_e_ge, l_err_ge, l_e_pe, l_err_pe




'''
emp1 = Employee('Frode', 2100)
emp2 = Employee('Johhny', 1300)

emp2.displayEmployee()
emp2.displaySalary()
print Employee.empcount

print "Employee.__doc__:", Employee.__doc__
print "Employee.__name__:", Employee.__name__
print "Employee.__module__:", Employee.__module__
print "Employee.__bases__:", Employee.__bases__
print "Employee.__dict__:", Employee.__dict__i
'''
