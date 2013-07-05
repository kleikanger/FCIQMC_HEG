#ifndef WALKERCONTAINERCLASS_H
#define WALKERCONTAINERCLASS_H

#include "imputVars.h"
#include "sortWalkers.h"
#include "walkerDistribution.h"
#include "loadBalanceThreads.h"
#include <vector>

using std::cerr;
using std::bitset;
using std::vector;

		struct loadBalThr
		{
			unsigned long first;
			unsigned long ndet;
			long nwfirst_thrd;
			long nwlast_thrd;
		};

template<size_t N>
class walkerContainerClass
{/*startvimfold*/

	private:
		long l_numwalkers_node;
		long l_numwalkers_global;
		long l_numdets_global;
		long l_ndet_node;
		long l_nnew_node;
		int i_numthreads;
		
		long* pl_temp;

		sortWalkers<N> o_sortwalkers;
		walkerDistribution<N> o_walkerdistr;
		imputVars vars;
		loadBalanceThreads loadbal;

	public:

		//OBS public vars!	
		long** ppl_new;
		bitset<N>* pb_new;
		long *pl_population;
		bitset<N> *pb_slater;
		vector< loadBalThr > wdistthr;
		
		walkerContainerClass(int, int, imputVars&);
		~walkerContainerClass();
		void init();
		void initiateWalkers(unsigned int, bitset<N>);
		long findNumWalkers();
		void reDistributeSortMerge();
		void loadBalance();

		void SETnnew(long l_nnew_node){ this->l_nnew_node = l_nnew_node; }
		void SETndet(long l_ndet_node){ this->l_ndet_node = l_ndet_node; }
		long lGETndet() const { return l_ndet_node; }
};/*endvimfold*/
#endif
