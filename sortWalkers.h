#ifndef SORTWALKERS_H
#define SORTWALKERS_H

#include <bitset>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <new>
#include "newmatrix.h"
//#include "last_run_definitions.h"

//#include <omp.h> //TODO NECC HERE?

#ifndef BITSET_BOUNDSCHECK
	#define BITSET_BOUNDSCHECK false
#endif

#if BITSET_BOUNDSCHECK
	#define SET set
	#define RESET reset
	#define TEST test
#else
	#define SET _Unchecked_set
	#define RESET _Unchecked_reset
	#define TEST _Unchecked_test
#endif

//#ifndef PARALLELSORT
#define PARALLELSORT false
//#endif

#if PARALLELSORT
	#include <parallel/algorithm>
	//using __gnu_parallel::sort;
	#define SORT __gnu_parallel::sort
#else
	#include <algorithm>
	//using std::sort;
	#define SORT std::sort
#endif

//#include <functional>

using std::bitset;
using std::cerr;

template <size_t N>
class sortWalkers
{
	private:

		bitset<N>* pb_newsdtemp;
		bitset<N>* pb_temp;
		unsigned long* ul_index;
		long** ppl_temp;
		long* pl_temp;
		
		bitset<N>* pb_sdswap;
		long** ppl_swap; 
		long* pl_swap;


		/*!
		 * Her. from std::binary_function<> functor, used to overload the sort function.
		 * <Less than> for two bitset's in the array pb_dets (pointer passed to the constructor) with
		 * the <const unsigned long indices> that  is passed to the functor.
		 *
		 */
		class LessThan 
		{
			private:
				bitset<N> b_tmp;
				unsigned int l_tmp;
				bitset<N>* pb_dets;
			public:
			LessThan (bitset<N>* pb_dets) { this->pb_dets = pb_dets; }

			bool operator()(const unsigned long& lhs, const unsigned long& rhs)
			{
				b_tmp = (pb_dets[rhs]^pb_dets[lhs]);
				l_tmp = b_tmp._Find_first();
				if (l_tmp==N) 
					return false;
				else if (pb_dets[rhs].TEST(l_tmp)) 
					return false;
				else 
					return true;
			}
		};


	public:

		/*!
		 * Constructor / Destructor.
		 */
		sortWalkers(const int n, const int m)
		{
			pb_newsdtemp = new bitset<N>[n];
			ul_index = new unsigned long[n];
			pb_temp = new bitset<N>[m];
			ppl_temp = (long**) matrix(n, 2, sizeof(long));
			pl_temp = new long[m];
		}
		/*!
		 *
		 */
		~sortWalkers()
		{
			delete [] pb_newsdtemp;
			delete [] pb_temp;
			delete [] ul_index;
			delete [] pl_temp;
			free_matrix((void**)ppl_temp);
		}

		/*!
		 *
		 * Sort the indices array s.t. to find the sorted ordering of the 
		 * array of the slater determinants.
		 *
		 */
		void sortIndices(bitset<N>* pb_dets, const unsigned long n)
		{
			for (unsigned long i=0; i<n; i++)
			{
				ul_index[i] = i;
			}

			SORT(&ul_index[0], &ul_index[n], LessThan(pb_dets) );
//			SORT(&ul_index[0], &ul_index[n],
//					//Sigma function (c++11), overloading <
//					[&k, &test, &pb_dets](const unsigned long& rhs, const unsigned long& lhs) -> bool
//					{
//						test = (pb_dets[lhs]^pb_dets[rhs]);
//						k = test._Find_first();
//						if ( (k==N) || (pb_dets[lhs].TEST(k)) ) 
//							return false;
//						else 
//							return true;
//					}
//				);
		}

		/*!
		 *
		 * Reorder the new walkers and merge equal determinants. Update the number of determinants. 
		 *
		 * The entire determinants array must be passed, Otherwise there will be problems
		 * when swapping the pointers.
		 *
		 * How do we parallelize this loop (openmp)
		 *
		 */
		void reorderWalkers(bitset<N>* &pb_dets, long** &ppl_occ, long &n)
		{
			long i = 0, j = 0;
			
			while (i<n)
			{
				pb_newsdtemp[j] = pb_dets[ul_index[i]];
				ppl_temp[j][0] = ppl_occ[ul_index[i]][0];
				ppl_temp[j][1] = ppl_occ[ul_index[i]][1];
				i++;

				while ((pb_dets[ul_index[i]]==pb_newsdtemp[j]) && (i<n))
				{
					ppl_temp[j][0] += ppl_occ[ul_index[i]][0];
					ppl_temp[j][1] += ppl_occ[ul_index[i]][1];
					i++;
				}
				
				if ( (ppl_temp[j][0]!=0) || (ppl_temp[j][1]!=0) ) 
					j++; //if pop = 0, remove the determinant
			}
			n = j; //update the number of determinants

			//Pointer swapping
			pb_sdswap = pb_dets;
			pb_dets = pb_newsdtemp;
			pb_newsdtemp = pb_sdswap;
			
			ppl_swap = ppl_occ;
			ppl_occ = ppl_temp;
			ppl_temp = ppl_swap;

		}

		void mergeNewAndOld( 
				bitset<N>* &pb_sd, long* &pl_occupancy, 
				long &l_ndet, bitset<N>* &pb_newdets, long** &ppl_occ, 
				long &l_nnew) 
		{
			long i = 0, j = 0, n = 0, l_temp;
			//Sort initiators and non initiators
			while ( (i<l_ndet) && (j<l_nnew) )
			{
				//pb_newdets[i] < pb_sd[j] ?
				const bitset<N> test = (pb_sd[i]^pb_newdets[j]);
				const int k = test._Find_first();
				if (k==N) //( pb_sd[i] == pb_newdets[j] )
				{

					pb_temp[n] = pb_sd[i];
					//old determinants
					pl_temp[n] = pl_occupancy[i];
					//initiators
					pl_temp[n] += ppl_occ[j][0];
					//non initiators can always spawn on previoulsy occupied walkers
					//if ( (pl_occupancy[i]!=0) || (abs(ppl_occ[j][1])>1) )
					pl_temp[n] += ppl_occ[j][1];
					
					i++;
					j++;
				}
				else if
					(!pb_sd[i].TEST(k)) //true
					{
						l_temp = 0;
						//initiators can always spawn everywhere
						l_temp += ppl_occ[j][0];
						//non-init can only spawn on non-occupied determinants 
						//if the |sum of non initiator spawning events| >1
						if ( abs(ppl_occ[j][1]) > 1 )
						{
							l_temp += ppl_occ[j][1];
						}

						pb_temp[n] = pb_newdets[j];
						pl_temp[n] = l_temp;
						j++;

					}
				else //false
				{
					//old determinants
					pb_temp[n] = pb_sd[i];
					pl_temp[n] = pl_occupancy[i];
					i++;
				}
				// if the population is larger than 0, accept the determiant 
				if ( pl_temp[n] != 0 )
				{
					n++;
				}
			}
			//all the rest are either initiators or non initiators XXX ??? or previously occupied!!!
			if (i==l_ndet)
			{
				while (j<l_nnew)
				{
					pb_temp[n] = pb_newdets[j];
					pl_temp[n] = ppl_occ[j][0];
					if ( abs(ppl_occ[j][1])>1 )
					{
						pl_temp[n] += ppl_occ[j][1];
					}
					j++;
					if ( pl_temp[n] != 0 )
					{
						n++;
					}
				}
			}
			else if ( j==l_nnew )
			{
				while ( i<l_ndet)
				{
					pb_temp[n] = pb_sd[i];
					pl_temp[n] = pl_occupancy[i];
					i++;
					if ( pl_temp[n] != 0 )
					{
						n++;
					}
				}
			}

			//reset number of determinants
			l_ndet = n;
			l_nnew = 0; //XXX NECC HERE?

			//Pointer swapping
			pb_sdswap = pb_sd;
			pb_sd = pb_temp;
			pb_temp = pb_sdswap;

			pl_swap = pl_occupancy;
			pl_occupancy = pl_temp;
			pl_temp = pl_swap;
		}
		/*!
		 * 
		 *
		 * If pb_newdets and pl_newocc must are sorted, the final list will be as well.
		 * XXX integrate in walkerDistribution class.
		 *
		 */
		void mergeRecvAndOld(
				long* &pl_occupancy, bitset<N>* &pb_sd, 
				long* pl_newocc, bitset<N>* pb_newdets, 
				long const l_nrecv_first, long const l_nrecv_last, long &l_ndet)
		{
			long j = 0, i, s;
			const long t = (l_nrecv_last >= 0) ? 0 : l_nrecv_last; // 0 or negative

			// first n walkers removed
			if (l_nrecv_first<=0)
			{
				for (i=abs(l_nrecv_first); i<l_ndet+t;i++)
				{
					pb_temp[j] = pb_sd[i];
					pl_temp[j] = pl_occupancy[i];
					++j;
				}
				s = 0;
			}
			// first n walkers recv
			else //if (l_nrecv_first>0)
			{
				s = l_nrecv_first; //>0
				for (i=0;i<s;i++)
				{
					pb_temp[j] = pb_newdets[i];
					pl_temp[j] = pl_newocc[i];
					++j;
				}
				for (i=0; i<l_ndet+t;i++)
				{
					pb_temp[j] = pb_sd[i];
					pl_temp[j] = pl_occupancy[i];
					++j;
				}
			}
			// last n walkers recv
			for (i=s;i<s+l_nrecv_last;i++) //only looping if l_nrecv_last>0
			{
				pb_temp[j] = pb_newdets[i];
				pl_temp[j] = pl_newocc[i];
				++j;
			}
			l_ndet = j;
			
			//Pointer swapping
			pb_sdswap = pb_sd;
			pb_sd = pb_temp;
			pb_temp = pb_sdswap;

			pl_swap = pl_occupancy;
			pl_occupancy = pl_temp;
			pl_temp = pl_swap;
		}
};
#endif
