#ifndef WALKERCONTAINERCLASS_CPP
#define WALKERCONTAINERCLASS_CPP

#include <mpi.h>
#include "walkerContainerClass.h"
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::bitset;

template<size_t N>
walkerContainerClass<N>::walkerContainerClass(int i_myrank, int i_nprocs, imputVars &vars) :
	o_sortwalkers(vars.iGETmaxndets()/5, vars.iGETmaxndets()),
	o_walkerdistr(i_myrank, i_nprocs, vars),
	loadbal(vars),
	//vlb(vars.iGETnumthreads(), vector<long>(6))
	wdistthr(vars.iGETnumthreads())
{
	this->vars = vars;
	pb_slater = NULL;
	pl_population = NULL;
	pl_temp = NULL;
	pb_new = NULL;
	ppl_new = NULL;
}
template<size_t N>
walkerContainerClass<N>::~walkerContainerClass()
{
	if (pb_slater)
	{
		delete [] pb_slater;
		delete [] pl_population;
		delete [] pl_temp;
		delete [] pb_new;
		free_matrix((void**)ppl_new);
	}
}
template<size_t N>
void walkerContainerClass<N>::init()
{
	int i_maxndets = vars.iGETmaxndets();
	i_numthreads = vars.iGETnumthreads();

	try
	{
		//containing the walkers 
		pb_slater= new bitset<N>[i_maxndets];
		pl_population = new long[i_maxndets];
		//temporary array
		pl_temp = new long[i_maxndets/5];
		pb_new = new bitset<N>[i_maxndets/5];
		ppl_new = (long**) matrix(i_maxndets/5, 2, sizeof(long));
	}
	catch (std::bad_alloc& ba)
	{
		cerr 
			<< "\nERROR: Memory allocation failed in walkerContainerClass!"
			<< "what(): " << ba.what()
			<< "\nTerminating run.\n";
		exit(0);
	}
}
/*!
 *
 * Set the initial state. One determinant if the i_initnumwalkers!=0 else
 * 0 determinants.
 *
 */
template<size_t N>
void walkerContainerClass<N>::initiateWalkers(
		unsigned int i_initnumwalkers, bitset<N> b_initial_state = 0)
{
	loadbal.init(i_initnumwalkers);
	if (i_initnumwalkers!=0)
	{
		l_ndet_node = 1;
		l_numwalkers_node = i_initnumwalkers;
		pb_slater[0] = b_initial_state;
		pl_population[0] = i_initnumwalkers;
	}
	else
	{
		l_ndet_node = 0;
		l_numwalkers_node = 0;
	}
}

template<size_t N>
long walkerContainerClass<N>::findNumWalkers()
{
	l_numwalkers_node = 0;
	for (int i=0; i<l_ndet_node; i++)
		l_numwalkers_node += abs(pl_population[i]);
	return l_numwalkers_node;
}
/*!
 *
 * Sort, merge, redistribute
 *
 */
template<size_t N>
void walkerContainerClass<N>::reDistributeSortMerge()
{
	if (l_nnew_node>1)
	{
		o_sortwalkers.sortIndices(pb_new, l_nnew_node);
		o_sortwalkers.reorderWalkers(pb_new, ppl_new, l_nnew_node);
	}
	// mpi : distrubute the new walkers
	MPI_Barrier(MPI_COMM_WORLD); //TODO NECC?
	o_walkerdistr.findNodeLimits(l_ndet_node, pb_slater); // FIXME only necc after redistribution of walkers and in the initialization !!!!
	MPI_Barrier(MPI_COMM_WORLD); //TODO NECC?
	o_walkerdistr.distributeNewWalkers(ppl_new, pb_new, l_nnew_node);
	// sort the redistributed walkers
	if (l_nnew_node>1)
	{
		o_sortwalkers.sortIndices(pb_new, l_nnew_node);
		o_sortwalkers.reorderWalkers(pb_new, ppl_new, l_nnew_node);
	}
	// merge the old and the new determiant lists and remove dets with 0 population
	o_sortwalkers.mergeNewAndOld(
			pb_slater, pl_population, l_ndet_node, 
			pb_new, ppl_new, l_nnew_node);
	// mpi : find new weights and redistribute determinants if necc.

	findNumWalkers();

	o_walkerdistr.setParam(l_ndet_node, l_numwalkers_node);//XXX NOT NECC?
	o_walkerdistr.findIdealWeight();
	
	if (o_walkerdistr.reDistribute())
	{
		long l_nrecv_first, l_nrecv_last;
		o_walkerdistr.findSendToLastNext(pl_population, pb_slater);
		o_walkerdistr.redistributeDeterminants(
				pl_population, pb_slater, pl_temp,
				pb_new, l_nrecv_first, l_nrecv_last);
		o_sortwalkers.mergeRecvAndOld(
				pl_population, pb_slater, pl_temp, pb_new,
				l_nrecv_first, l_nrecv_last, l_ndet_node);
	}
}
/*!
 *
 * Load balance on threads. Make sure that l_numwalkers_node and l_ndet_node is updated.
 *
 */
template<size_t N>
void walkerContainerClass<N>::loadBalance()
{
	loadbal.loadBalance(pl_population, l_numwalkers_node, l_ndet_node);

	for (int i=0;i<i_numthreads;i++)
		{
			//wdistthr[i] = loadbal.GETwdistthr(i);
			wdistthr[i].first = loadbal.ulGETfirstdet_thread(i);
			wdistthr[i].ndet = loadbal.ulGETndet_thread(i);
			wdistthr[i].nwfirst_thrd = loadbal.lGETnwfirst_thread(i); 
			wdistthr[i].nwlast_thrd = loadbal.lGETnwlast_thread(i);
		}

}
#endif
