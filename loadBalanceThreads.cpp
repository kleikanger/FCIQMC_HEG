#include "loadBalanceThreads.h"
#include <vector>
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::vector;
using std::abs;

/*!
 * 
 * Initiate with 1 dets and i_initnumwalkers walkers if 
 * i_initnumwalkers != 0 and 0 dets otherwise.
 *
 */
void loadBalanceThreads::init(unsigned int ui_initnumwalkers)
{/*startvimfold*/
	for (int i=0;i<i_numthreads;i++)
	{
		pul_ndet_thread[i] = 0;
		pul_firstdet_thread[i] = 0;
		pl_nwfirst_thread[i] = 0;
		pl_nwlast_thread[i] = 0;
	}
	if (ui_initnumwalkers!=0)
	{
		pul_ndet_thread[0] = 1;
		pul_firstdet_thread[0] = 0;
		pl_nwfirst_thread[0] = ui_initnumwalkers;
		pl_nwlast_thread[0] = ui_initnumwalkers;
	}
}/*endvimfold*/

/*!
 *
 * Load balancing between the threads.
 * Distribute the determinants and the walkers between the threads. 
 * The load should be the same on all threads according to the formula
 * 
 * 		weight = num_walkers + num_determinants * d_detweight.
 *
 * pul_firstdet_thread[i]	: The first determinant to be evaluated by thread i
 * pul_ndet_thread[i] 		: The number of determinants to be evaluated by thread i
 * (The number of walkers to be evaluated by a thread in not always the same as the 
 * 	occupancy of the thread pi_occ Since a determinant can be split between two threads:)
 * pl_nwfirst_thread[i] 	: The number of w's on the first det. to be eval. by thread i
 * pl_nwlast_thread[i] 		: The number of w's on the last det. to be eval. by thread i
 * (We neec to know the total nr of dets on a thread to know it the thread is an initiator)
 *
 *
 */
void loadBalanceThreads::loadBalance(
		long* pl_population, 
		long l_numwalkers_node, 
		long l_ndet_node)
{/*startvimfold*/

	unsigned long ul_nwlast_node;
		// if no dets, set all to 0 and exit method
		if (l_ndet_node==0)
		{
			for (int i=0;i<i_numthreads;i++)
			{
				pul_firstdet_thread[i] = 0;
				pl_nwlast_thread[i] = 0;
				pul_ndet_thread[i] = 0;
			}
			return;
		}

		double d_thread_weight_cum, d_leftover_c;
		long l_nwthisdet, i_countdets, i_tcd;
		int k, i_lastsplit;

		const double d_thread_weight =
			static_cast<double>(l_numwalkers_node + d_detweight*l_ndet_node)
			/ static_cast<double>(i_numthreads);
		// in the case of fewer walkers than threads 1 is set to be the the minimum value

		pul_firstdet_thread[0] = 0;
		pl_nwfirst_thread[0] = pl_population[0];

		i_tcd = 0; k = 0; l_nwthisdet = 0; i_lastsplit = 0;
		i_countdets = 1; 
		d_thread_weight_cum = d_detweight;
		d_leftover_c = d_detweight;

		while (k<i_numthreads-1)
		{

			if (i_tcd<l_ndet_node) // enough walkers on thread k
				d_thread_weight_cum += static_cast<double>(abs(pl_population[i_tcd])-l_nwthisdet);
			else // all walkers distributed
				break;
			// when the number of walkers/dets is high enough, set the number of walkers/dets on this 
			// thread and change to the next thread (++k)
			if (d_thread_weight_cum<d_thread_weight)
			{
				++i_countdets;
				++i_tcd;
				d_thread_weight_cum += d_detweight;
				i_lastsplit = 1;
			}
			else //(d_thread_weight_cum>=d_thread_weight)
			{
				//there is a possibility that l_lo<0
				long l_lo = static_cast<int>(d_thread_weight_cum - d_thread_weight - d_leftover_c);
				if (l_lo<0) l_lo = 0;

				pul_firstdet_thread[k+1] = i_tcd;
				ul_nwlast_node = abs(pl_population[i_tcd]);
				pl_nwlast_thread[k] = ul_nwlast_node - l_lo - l_nwthisdet;

				pul_ndet_thread[k] = i_countdets;
				pl_nwfirst_thread[k+1] = ul_nwlast_node - pl_nwlast_thread[k] - l_nwthisdet;

				pl_nwlast_thread[k] *= ( (pl_population[i_tcd]>0) ? 1 : -1);
				pl_nwfirst_thread[k+1] *= ( (pl_population[i_tcd]>0) ? 1 : -1);

				if (static_cast<double>(abs(pl_population[i_tcd])-l_nwthisdet-abs(pl_nwlast_thread[k]))>=d_thread_weight)
				{
					i_lastsplit = 0;
					d_thread_weight_cum = 0;
					i_countdets = 1;
					l_nwthisdet += abs(pl_nwlast_thread[k]); //must have same abs val as nw_last_thread[k+1] 
					d_leftover_c = 0;
				}
				else 
				{
					i_lastsplit = 1;
					d_thread_weight_cum = static_cast<double>(abs(pl_nwfirst_thread[k+1])) + d_detweight;
					d_leftover_c = d_detweight; 
					++i_tcd;
					if (i_tcd<l_ndet_node) //&& (k<i_numthreads-1) 
					{
						l_nwthisdet = 0;
						i_countdets = 2;
					}
					else
					{
						l_nwthisdet += abs(pl_nwlast_thread[k]); 
						i_countdets = 1;
					}
				}
				++k;
			}
		}

		pul_ndet_thread[k] = l_ndet_node-i_tcd+i_lastsplit;
		pl_nwlast_thread[k] = abs(pl_population[l_ndet_node-1]) - l_nwthisdet;
		pl_nwlast_thread[k] *= (pl_population[l_ndet_node-1]>0) ? 1. : -1.; 

		// set number of determinants in the rest of the threads to 0
		++k;
		while (k<i_numthreads)
		{
			pul_firstdet_thread[k] = 0;
			pul_ndet_thread[k] = 0;
			pl_nwfirst_thread[k] = 0; // NOT NECC
			pl_nwlast_thread[k] = 0; // NOT NECC
			++k;
		}

}/*endvimfold*/
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
