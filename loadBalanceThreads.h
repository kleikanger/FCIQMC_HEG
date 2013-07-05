#ifndef LOADBALANCETHREADS_H
#define LOADBALANCETHREADS_H

#include "imputVars.h"
#include <vector>

using std::vector;

class loadBalanceThreads 
{
	private:
		vector<unsigned long> pul_firstdet_thread;
		vector<unsigned long> pul_ndet_thread;
		vector<long> pl_nwfirst_thread;
		vector<long> pl_nwlast_thread;
		const double d_detweight;
		const int i_numthreads;

	public:
		/*
		 * constructor
		 */
		loadBalanceThreads(imputVars &vars) :
			pul_firstdet_thread(vars.iGETnumthreads()),
			pul_ndet_thread(vars.iGETnumthreads()),
			pl_nwfirst_thread(vars.iGETnumthreads()),
			pl_nwlast_thread(vars.iGETnumthreads()),
			d_detweight(vars.dGETdetweight()),
			i_numthreads(vars.iGETnumthreads())
		{}
		
		unsigned long ulGETfirstdet_thread(int i_thread) const
			{ return pul_firstdet_thread[i_thread]; }
		unsigned long ulGETndet_thread(int i_thread) const
			{ return pul_ndet_thread[i_thread]; }
		long lGETnwfirst_thread(int i_thread) const
			{ return pl_nwfirst_thread[i_thread]; }
		long lGETnwlast_thread(int i_thread) const
			{ return pl_nwlast_thread[i_thread]; }
		/*
		 * class methods in source file
		 */	
		void init(unsigned int);
		void loadBalance(long*, long, long);
};

#endif
