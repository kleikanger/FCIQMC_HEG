#ifndef WALKERDISTRIBUTION_H
#define WALKERDISTRIBUTION_H

#include <mpi.h>  //TODO NECC HERE?

#include <bitset>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "newmatrix.h"
#include "imputVars.h"

#include <iostream> 

//#include <omp.h> //TODO NECC HERE?

#ifndef BITSET_BOUNDSCHECK
	#define BITSET_BOUNDSCHECK true
#endif

#ifndef SET
#if BITSET_BOUNDSCHECK
	#define SET set
	#define RESET reset
	#define TEST test
#else
	#define SET _Unchecked_set
	#define RESET _Unchecked_reset
	#define TEST _Unchecked_test
#endif
#endif

using std::bitset;
using std::cerr;

template <size_t N>
class walkerDistribution
{

	private:
		int	i_numthreads_node; //const
		int	i_numthreads_global; //const
		int	i_nprocs; //const

		double d_detweight;
		double d_idealweight;
		double d_redistlimit;

		long l_numwalkers_node;
		long l_numwalkers_global;
		long l_numdets_node;
		long l_numdets_global;

		double *pd_nodeweights;

		int i_thread_id; //const
		int i_myrank; //const

		//bitset<N>* pb_firstdet;
		bitset<N>* pb_lastdet;

		long* pl_sendtolast;
		long* pl_sendtonext;

		//size_t t_bsbytes;
		bitset<N> b_empty;

		//MPI_Status mpi_status;

		long** ppl_temp;
		long** ppl_swap;
		bitset<N>* pb_temp;
		bitset<N>* pb_swap;

		MPI_Datatype mpi_bitsetn;
							
		//lessthan
		class LessThan 
		{
			private:
				unsigned int k;
				bitset<N> b_tmp;
			public:
				bool operator()	(const bitset<N>& lhs, const bitset<N>& rhs)
				{
					b_tmp = (lhs^rhs);
					k = b_tmp._Find_first();
					if (k==N) 
						return false; //true:<=, false:<
					else if (rhs.TEST(k))
						return false;
					else
						return true;
				}
		};

	public:
		/*
		 * 
		 * Constructor
		 *
		 */
		walkerDistribution(int i_myrank, int i_nprocs, imputVars &vars)
		{/*startvimfold*/
			this->i_nprocs = i_nprocs;
			this->i_myrank = i_myrank;
			i_numthreads_node = vars.iGETnumthreads();
			d_detweight = vars.dGETdetweight();
			d_redistlimit = vars.dGETredistlimit();

			long i_nndmax = (vars.iGETmaxndets())/5 ;
			ppl_temp = (long**) matrix(i_nndmax, 2, sizeof(long));
			try
			{
				pd_nodeweights = new double[i_nprocs];
				pb_lastdet = new bitset<N>[i_nprocs];
				//pb_firstdet = new bitset<N>[i_nprocs];
				pl_sendtolast = new long[i_nprocs];
				pl_sendtonext = new long[i_nprocs];
				//temporary arrays, can be large.
				pb_temp = new bitset<N>[i_nndmax];
			}
			catch (std::bad_alloc&)
			{
				cerr 
					<< "\nERROR: Memory allocation failed in walkerDistribution!"
					<< "\nTerminating run.\n";
				exit(0);
			}
			
			// for comparison in fome funcs
			b_empty.reset();
			// find the global number of threads
			MPI_Allreduce(&i_numthreads_node, &i_numthreads_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			// MPI derived type for bitset
			mpi_bitsetn = makeDerivedtype();

		};/*endvimfold*/
		/*
		 *
		 * Destructor
		 *
		 */
		~walkerDistribution()
		{/*startvimfold*/
			delete [] pd_nodeweights;
			delete [] pb_lastdet;
			//delete [] pb_firstdet;
			delete [] pl_sendtolast;
			delete [] pl_sendtonext;
			delete [] pb_temp;
			delete [] ppl_temp;
		};/*endvimfold*/
		
		/*
		 *
		 *
		 */
		void setParam(long l_nd_n, long l_nw_n)
		{/*startvimfold*/
			l_numdets_node = l_nd_n;
			l_numwalkers_node = l_nw_n;
		}/*endvimfold*/

		MPI_Datatype makeDerivedtype() {
			//int err = 0;/*startvimfold*/
			MPI_Datatype t;
			// This is always guaranteed to work.  sizeof(T) for any type T
			// must be a nonnegative integer, and thus must be a multiple of
			// sizeof(char) (which == 1 by the C++03 standard).
			//err = 
			MPI_Type_contiguous (sizeof(bitset<N>), MPI_CHAR, &t);
			/*	if (err != MPI_SUCCESS) 
			 *	{} 
			 */
			//err = 
			MPI_Type_commit (&t);
			/* if (err != MPI_SUCCESS )
			 * {}
			 * */ 
			return t;
		}/*endvimfold*/

		/*
		 *
		 * Find the weight distribution across the nodes that is stored in the global array pd_nodeweights, 
		 * and the ideal weight distribution that is stored in the global variable d_idealweight. 
		 * pd_nodeweights and d_idealweight are equal on all nodes after this.
		 *
		 * Remember to setParam before running this func.
		 *
		 */
		void findIdealWeight()
		{/*startvimfold*/
			double d_nodeweight = d_detweight*l_numdets_node + l_numwalkers_node;
			// all processes receives the messages and stores them in rank order
			MPI_Allgather(&d_nodeweight, 1, MPI_DOUBLE, pd_nodeweights, 1, MPI_DOUBLE, MPI_COMM_WORLD);

			d_idealweight = 0;
			for (int i=0;i<i_nprocs;i++) 
			{
				d_idealweight += pd_nodeweights[i];
			}
			// assuming (wrongly) that the code scales linearly with the 
			// number of threads on a node
			d_idealweight *=  i_numthreads_node / static_cast<double>(i_numthreads_global);
		}/*endvimfold*/
		/*
		 *
		 * After finding the ideal weight, check if a redistribution of the determinants are necc.
		 *
		 */
		bool reDistribute()
		{/*startvimfold*/
			//Check if the determinants should be redistributed
			for (int i=0;i<i_nprocs;i++)
			{
				if (d_idealweight>(pd_nodeweights[i]*(1.+d_redistlimit)))
					return true;
				else if (d_idealweight<(pd_nodeweights[i]*(1.-d_redistlimit)))
					return true;
				else 
					return false;
			}
			return false;
		}/*endvimfold*/

		/*
		 *
		 * Run findIdealWeight first to find the weight distribution across the nodes and d_idealweight.
		 *
		 * This method finds the number of walkers that should be sent to the next/last node
		 * with rank +1/-1.  It is a sequential algorithm, which means that all other nodes stops running
		 * till it recieves the number of walkers that should be sent to/from the last node.
		 *
		 * These numbers are stored in the class variables l_sendtolast and l_sendtonext.
		 *
		 * Remember to setParam before running this func.
		 *
		 * It would be much better to find a method which was not sequental.
		 *
		 * XXX OPTIMIZATIONS : Ideas to implement later XXX
		 *
		 * - Start sending from the first node that has a ''critical'' population
		 *   and stop sending at the firs node after the last ''critical'' with a 
		 *   weight that is ''non-critical''
		 *
		 * - Start calculations from the last and the first node simultaniously. 
		 *   Join on the middle.
		 *
		 */
		void findSendToLastNext(long* pl_occupancy, bitset<N>* pb_sd)
		{/*startvimfold*/
			MPI_Status mpi_status;
			double d_wtmp, d_nextnodediff, d_lastnodediff;
			long i, l_sendtolast, l_sendtonext;
				
			d_wtmp = 0;
			i = 0;

			if (i_myrank==0)
			{
				d_lastnodediff = 0.;
				l_sendtolast = 0;
			}
			else //if (i_myrank!=0) 
			{
				// every process except root waits here till it recieves the variable 
				// d_lastnodediff (d_nextnodediff is sent from the node with lower rank)
				MPI_Recv(&d_lastnodediff, 1, MPI_DOUBLE, i_myrank-1, 0, MPI_COMM_WORLD, &mpi_status);
				
				// find l_sendtolast 
				if (d_lastnodediff>0)
				{
					while ( (d_wtmp<d_lastnodediff) && (i<l_numdets_node) )
					{
						d_wtmp += static_cast<double>(abs(pl_occupancy[i])) + d_detweight; 
						++i;
					}
					l_sendtolast = i;
					d_wtmp = 0.;
				}
				else //if (d_lastnodediff<=0)
				{
					l_sendtolast = 0;
					d_wtmp = -d_lastnodediff; //positive
				}
			}

			// add the number of walkers to keep on this node to i
			while ( (d_wtmp<d_idealweight) && (i<l_numdets_node) )
			{
				d_wtmp += static_cast<double>(abs(pl_occupancy[i])) + d_detweight; 
				++i;
			}

			// find l_sendtonext
			if (i==l_numdets_node)
			{
				// too few walkers on node. transfer the first 
				// determiants and their population from the next node to this node
				l_sendtonext = 0;
				d_nextnodediff = std::max(d_idealweight - d_wtmp, 0.0); //TODO
			}
			else //if (d_wtmp>=d_idealweight)
			{
				// too many walkers on node
				//transfer walkers to the next node
				l_sendtonext = 0;//l_numdets_node - l_sendtolast - i;//XXX OBS - l_sendtolast ???
				d_nextnodediff = 0;// -l_sendtonext * d_detweight;
				while (i<l_numdets_node) 
				{
					d_nextnodediff -= abs(pl_occupancy[i]) + d_detweight;
					++l_sendtonext;
					++i;
				}
			}
			
			// sent to the next process
			if (i_myrank!=i_nprocs-1)
			{
				MPI_Send(&d_nextnodediff, 1, MPI_DOUBLE, i_myrank+1, 0, MPI_COMM_WORLD);
			}
			else //if (i_myrank==i_nprocs-1) 
			{
				l_sendtonext = 0;//NECC?
			}

#if 0		//XXX XXX XXX
			if (i_myrank==0)
			{
				for (i=0;i<i_nprocs;i++)
					std::cout << i << " stl " << pl_sendtolast[i] << " stn " << pl_sendtonext[i] << "\n";
			}
#endif			//XXX XXX XXX
			MPI_Allgather(&l_sendtolast, 1, MPI_LONG, pl_sendtolast, 1, MPI_LONG, MPI_COMM_WORLD);
			MPI_Allgather(&l_sendtonext, 1, MPI_LONG, pl_sendtonext, 1, MPI_LONG, MPI_COMM_WORLD); // XXX NECC?

		}/*endvimfold*/
		
		/*
		 *
		 * Redistribute the determinants to improve the work balance.
		 *
		 * NB: The vew dets are already sorted and larger or smaller than the present walkers
		 * This makes merging easy.
		 *
		 * l_sendtonext : changed to positive number if new dets are recieved at the end of determinants list, and
		 * 		to negative number if dets are sent from the end of the d.l.
		 * l_sendtolast : changed to positive number if new dets are recieved at the start of the d.l., and
		 * 		to negative number if dets are sent from the start of the d.l..
		 */
		void redistributeDeterminants(
				long* pl_occupancy, bitset<N>* pb_sd, 
				long* pl_recbuf, bitset<N>* pb_recbuf, 
				long &l_nrecv_first, long &l_nrecv_last)
		{/*startvimfold*/

			MPI_Status mpi_status;
			l_nrecv_first = l_nrecv_last = 0;

			// transmit last elem's to next proc
			if ( (pl_sendtonext[i_myrank]!=0) && (i_myrank<i_nprocs-1) )
			{
				const long l_indx = l_numdets_node - pl_sendtonext[i_myrank];
				// send pl_occupancy elem
				MPI_Send(&pl_occupancy[l_indx], pl_sendtonext[i_myrank], 
						MPI_DOUBLE, i_myrank+1, 0, MPI_COMM_WORLD);
				// send pb_sd (slater determinants) elem
				MPI_Send(&pb_sd[l_indx], pl_sendtonext[i_myrank], 
						mpi_bitsetn, i_myrank+1, 0, MPI_COMM_WORLD);
				// set sent occupations to 0
				for (int i=l_indx; i<l_indx+pl_sendtonext[i_myrank]; i++)//XXX NECC?
					pl_occupancy[i] = 0;
				l_nrecv_last = - pl_sendtonext[i_myrank];
			}
			if ( (pl_sendtonext[i_myrank-1]!=0) && (i_myrank>0) )
			{
				// recv pl_occupancy elem
				MPI_Recv(&pl_recbuf[0], pl_sendtonext[i_myrank-1], 
						MPI_DOUBLE, i_myrank-1, 0, MPI_COMM_WORLD, &mpi_status);
				// recv pb_sd (slater determinants) elem
				MPI_Recv(&pb_recbuf[0], pl_sendtonext[i_myrank-1], 
						mpi_bitsetn, i_myrank-1, 0, MPI_COMM_WORLD, &mpi_status);
				l_nrecv_first = pl_sendtonext[i_myrank-1];
			}
			// transmit first elem's to last proc
			if ( (pl_sendtolast[i_myrank]!=0) && (i_myrank>0) )
			{
				MPI_Send(&pl_occupancy[0], pl_sendtolast[i_myrank], 
						MPI_DOUBLE, i_myrank-1, 0, MPI_COMM_WORLD);
				MPI_Send(&pb_sd[0], pl_sendtolast[i_myrank], 
						mpi_bitsetn, i_myrank-1, 0, MPI_COMM_WORLD);
				for (int i=0; i<pl_sendtolast[i_myrank]; i++)//XXX NECC?
					pl_occupancy[i] = 0;
				l_nrecv_first = - pl_sendtolast[i_myrank];
			}
			if ( (pl_sendtolast[i_myrank+1]!=0) && (i_myrank<i_nprocs-1) )
			{
				const long l_indx = (l_nrecv_first>0) ? l_nrecv_first : 0;
				MPI_Recv(&pl_recbuf[l_indx], pl_sendtolast[i_myrank+1], 
						MPI_DOUBLE, i_myrank+1, 0, MPI_COMM_WORLD, &mpi_status);
				MPI_Recv(&pb_recbuf[l_indx], pl_sendtolast[i_myrank+1], 
						mpi_bitsetn, i_myrank+1, 0, MPI_COMM_WORLD, &mpi_status);
				l_nrecv_last = pl_sendtolast[i_myrank+1];
			}
		}/*endvimfold*/

		void distributeNewWalkers(long** &ppl_occ,  bitset<N>* &pb_newdets, long &l_nnew)
		{/*startvimfold*/
			bitset<N>* pb_tmp;
			bitset<N> b_tmp;
			bitset<N> b_upperbound;
			int i, k, j;
			long pl_firstto[i_nprocs];
			long pl_lastto[i_nprocs];//lastto not really neccessary XXX remove: is simply firstto + 1 
			long l_delim;
			int i_sendcount;
			int pi_sendcount_l[i_nprocs];
			int pi_sendcount_bs[i_nprocs];
			int i_displ_l[i_nprocs];
			int i_displ_bs[i_nprocs];
			bool b_pop;
			
			// find where to the first det should be sent
			for (i=0;i<i_nprocs-1;)//-1;)
			{
				// check if smallets new det is smaller than pb_firstdet[i+1]
//				if (pb_firstdet[i+1]!=b_empty)
//				{
//					b_tmp = (pb_newdets[0]^pb_firstdet[i+1]); 
//					k = b_tmp._Find_first();
//					if (k!=N) //not equal
//					{
//						if (pb_newdets[0].TEST(k)) //less than
//						{	
//							break;
//						}
//					}
//				}
				if (pb_lastdet[i]!=b_empty)
					break;
				//no walkers to this det
				pl_firstto[i] = pl_lastto[i] = 0;
				++i;
			}

			// find which walkers/dets to send to which node?
			// unpopulated nodes i where node i+1 is populated will 
			// recieve all dets between pb_lastdet and pb_firstdet 
			// of the populated nodes before and after.
			l_delim = 0;
			for (;i<i_nprocs-1;i++)
			{
				if (pb_lastdet[i]==b_empty)
				{
//					if (pb_lastdet[i+1]!=b_empty) 
//					{
//						b_upperbound = pb_firstdet[i+1];
//						pl_firstto[i] = l_delim;
//						b_pop = true;
//					}
//					else 
					{
						b_pop = false;
						pl_lastto[i] = pl_firstto[i] = l_delim;
					} 
				}
				else
				{
					b_upperbound = pb_lastdet[i];
					pl_firstto[i] = l_delim;
					b_pop = true;
				}

				// binary search. Finds the last bitset smaller than or equal to b_upperbound
				if (b_pop) 
				{	
					pb_tmp = 
						std::lower_bound(&pb_newdets[l_delim], &pb_newdets[l_nnew], 
								b_upperbound, LessThan() );
							//lambda func, overloading <
		//					[&k, &b_tmp](const bitset<N>& lhs, const bitset<N>& rhs) -> bool
		//					{
		//						b_tmp = (lhs^rhs);
		//						k = b_tmp._Find_first();
		//						if (k==N) 
		//							return false; // true:<=, false:<
		//						else if (rhs.TEST(k))
		//							return false;
		//						else
		//							return true;
		//					});
					pl_lastto[i] = ((pb_tmp - pb_newdets)>0) ? (pb_tmp - pb_newdets) :  0; //0 necc?
					//pl_lastto[i] = (pb_tmp - pb_newdets);
					if (pb_newdets[pl_lastto[i]]==b_upperbound)
						if (pl_lastto[i]<l_nnew)
							++pl_lastto[i];
					
					l_delim = pl_lastto[i];
				}
			}
			pl_firstto[i_nprocs-1] = l_delim; //XXX tro|bbel? can become larger an lastto
			pl_lastto[i_nprocs-1] = l_nnew;
//print
#if 0	
			MPI_Barrier(MPI_COMM_WORLD);
			for (k=0;k<i_nprocs;k++)
			{
				if (k==i_myrank)
				{
					for (i=0;i<i_nprocs;i++)
						std::cerr << "<" << pl_firstto[i] << ">";
					std::cerr << "\n";
					for (i=0;i<i_nprocs;i++)
						std::cerr << "<" << pl_lastto[i] << ">";
				}
				std::cout << std::endl;
				MPI_Barrier(MPI_COMM_WORLD);
			}
#endif
			//Gather the new determinants on the correct nodes	
			for (i=0;i<i_nprocs;i++)
			{
				i_sendcount = pl_lastto[i] - pl_firstto[i];
//std::cerr << "|" << i_sendcount << "|";

				// gather the num of elem sent from each proc
				MPI_Gather(
						&i_sendcount, 1, MPI_INT, 
						pi_sendcount_bs, 1, MPI_INT, 
						i, MPI_COMM_WORLD);
				// update l_nnew
				if (i_myrank==i)
				{
					l_nnew = 0;
					for (k=0;k<i_nprocs;k++)
						l_nnew += pi_sendcount_bs[k];
				}

				//std::cerr << ppl_occ[0][0] << " " << ppl_occ[0][1] << " " << ppl_occ[1][0] << "\n";
				//std::cerr << *(ppl_occ[0]) << " " << *(ppl_occ[0]+1) << " " << *(ppl_occ[0]+2) << "\n";

				// init the sendcount and dislplacement arrays for MPI_Gatherv
				j = 0;
				for (k=0;k<i_nprocs;k++)
				{
			   		i_displ_l[k] = j * 2;	
			   		i_displ_bs[k] = j; //TODO not necc use derived datatype	
					j += pi_sendcount_bs[k];
				}
				for (k=0;k<i_nprocs;k++)
					pi_sendcount_l[k] = pi_sendcount_bs[k] * 2;// 2 long's per det since pl_sendcount_l is a 2D array
				
				// gather the new determinants on the correct nodes
				j = pl_firstto[i];
				MPI_Gatherv(
					ppl_occ[j], i_sendcount * 2, MPI_LONG, // *2 since this is a 2d array
					ppl_temp[0], pi_sendcount_l, i_displ_l, MPI_LONG,
					i, MPI_COMM_WORLD);
				MPI_Gatherv(
					&pb_newdets[j], i_sendcount , mpi_bitsetn,
					pb_temp, pi_sendcount_bs, i_displ_bs, mpi_bitsetn,
					i, MPI_COMM_WORLD);
		/*
				// XXX is pb_temp changed only for the root proc??
				if (i_myrank==i)
				{
					// pointer swapping
					pb_swap = pb_newdets;
					pb_newdets = pb_temp;
					pb_temp = pb_swap;

					ppl_swap = ppl_occ;
					ppl_occ = ppl_temp;
					ppl_temp = ppl_swap;
				}
				*/
				
			}

			MPI_Barrier(MPI_COMM_WORLD); //TODO Necc??

			// pointer swapping
			pb_swap = pb_newdets;
			pb_newdets = pb_temp;
			pb_temp = pb_swap;

			ppl_swap = ppl_occ;
			ppl_occ = ppl_temp;
			ppl_temp = ppl_swap;
			
		}/*endvimfold*/
		/*
		 *
		 * Find the largest and smallest slaterdeterminant on each node and
		 * store them in the global arrays pb_firstdet and pb_lastdet.
		 *
		 */
		void findNodeLimits(int i_ndet_node, bitset<N>* pb_sd)
		{
			bitset<N> b_first, b_last; 
			
			if (i_ndet_node==0)
			{
				b_first = b_last = b_empty;
			}
			else 
			{
				b_first = pb_sd[0];
				b_last = pb_sd[i_ndet_node-1];
			}
#if 0
			if (i_myrank==0)
			{
				std::cout << std::endl;
				for (int i=0;i<i_nprocs;i++)
					std::cout << pb_firstdet[i] << " " << pb_lastdet[i] << " " << i << "-.,\n ";
				std::cout << std::endl;
			}
#endif
		//	MPI_Allgather(&b_first, t_bsbytes, MPI_CHAR, // TODO change name to nodelim
		//			pb_firstdet, t_bsbytes, MPI_CHAR, MPI_COMM_WORLD);
			//MPI_Allgather(&b_last, t_bsbytes, MPI_CHAR, // TODO delete. lastdet not necc
			//		pb_lastdet, t_bsbytes, MPI_CHAR, MPI_COMM_WORLD);
//			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgather(&b_last, 1, mpi_bitsetn, // TODO delete. lastdet not necc
					pb_lastdet, 1, mpi_bitsetn, MPI_COMM_WORLD);
		}

};
#endif





// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
