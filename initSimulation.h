#ifndef INITSIMULATION_H
#define INITSIMULATION_H

#include <mpi.h>
#include <omp.h>
#include <bitset>
#include <iostream>
#include <string>
#include "imputVars.h"
#include <cmath>
#include "last_run_definitions.h"

using std::string;
using std::cout;
using std::bitset;

template <size_t N>
class initSimulation
{
	private:
		bitset<N> b_initial_state;
		bitset<N> b_frozen_orbs;
		imputVars* vars;
		int i_myrank;
	public:
		initSimulation(imputVars* vars, int i_myrank)
		{/*startvimfold*/
			this->i_myrank = i_myrank;
			this->vars = vars;
			b_initial_state.reset();
			for (int i=0;i<vars->iGETnumpart();i++)
			{
				b_initial_state.set(i);
			}
			// XXX nb: set initial occ orbs as frozen orbs
			b_frozen_orbs = b_initial_state;
		} /*endvimfold*/
		~initSimulation()
		{}
	
		bitset<N> bGETinitial_state() 
		{ /*startvimfold*/
			return b_initial_state; 
		}/*endvimfold*/

		bitset<N> bGETfrozen_orbs()
		{ /*startvimfold*/
			return b_frozen_orbs; 
		}/*endvimfold*/

		void init()
		{/*startvimfold*/

			// set random seed = time if set to < 0 in vars
			int i_ranseed = vars->iGETranseed();
			if (i_ranseed<0)
			{
				if (i_myrank==0)
				{
					i_ranseed = static_cast<int>(std::abs(time(NULL)));
				}
				MPI_Bcast(&i_ranseed, 1, MPI_INT, 0, MPI_COMM_WORLD);
				vars->iSETranseed(i_ranseed);
			}
			// print param. to screen
			if (i_myrank==0)
			{
				cout << "\n"
					<< "Simulation tag=" << vars->sGETofpath() 
					<< ", "
					<< "Random seed=" << vars->iGETranseed() 
					<< ".\n"
					<< "N=" << vars->iGETnumpart()
					<< ", M=" << INUM_ORBITALS
					<< ", r_s=" << vars->dGETrs()
					<< ", N_I=" << vars->iGETinitiatorlimit()
					<< ", d_t=" << vars->dGETdt() << "\n";
				// cout << walkers.pb_slater[0];
			}
			// set the number of threads to the value specified in imputVars
			// if the enviromental variable OMP_NUM_THREADS is not set
			omp_set_dynamic(0); //fixed number of threads
			int i_tmp = omp_thread_count();
			if (i_tmp!=0)
			{
				vars->iSETnumthreads(i_tmp);
				omp_set_num_threads(i_tmp);
			}
			else
			{
				omp_set_num_threads(vars->iGETnumthreads());
			}
		}/*endvimfold*/
	
	private:

		int omp_thread_count() // FIXME Not working on abel??
		{/*startvimfold*/
			int n = 0;
#pragma omp parallel reduction(+:n)
			n += 1;
			return n;
		}/*endvimfold*/
};
#endif

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
