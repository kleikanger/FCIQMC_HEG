import matplotlib.pyplot as plt
import numpy as np
import pyDataAnal

num = 300
maxnw = a.i_limnw
numsims = 22
fnam = '6ptw1m0s0_'

def addsim(f):

	a  = pyDataAnal.pyDataAnal(f, 'diiid', 5); a.readBF2array()


	penum = np.zeros(num)
	pecum = np.zeros(num)

	for elem in a.a_data[:]:
		i = min(elem[1] * num // maxnw, num-1)
		pecum[i] += elem[4]
		penum[i] += 1

	return pecum/penum

res = np.zeros((numsims,num))
for i in range(1,numsims):
	f = '%s%i'%(fnam, i)
	res[i-1] = addsim(f)

err = np.zeros(num)

for i in range(num):
    err[i] = np.std(res[:-1,i])/np.sqrt(numsims)

#plt.figure()
y = np.zeros(num)
x = range(0,maxnw,maxnw/num)
for i in range (num):
    y[i] = res[:-1,i].mean()


plt.errorbar(x,y,yerr=err,fmt='k.')

plt.plot(x, y)
plt.show(block=False)


