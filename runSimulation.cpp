#ifndef RUNSIMULATION_CPP
#define RUNSIMULATION_CPP

#include <mpi.h>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "last_run_definitions.h"
#include "runSimulation.h"

// source must be included since these are template classes.
#include "walkerPropagator.cpp"
#include "walkerContainerClass.cpp"

//#include "eOutStream.h"

using namespace std;

template<size_t N>
runSimulation<N>::runSimulation(int i_myrank, int i_nprocs, imputVars* vars) 
	:
	initsim(vars, i_myrank),
	walkers(i_myrank, i_nprocs, *vars), 
	i_nprocs(i_nprocs),
	i_myrank(i_myrank),
	i_limit_nw(vars->iGETlimit_nw()),
	i_update_s_frequency(vars->iGETupdate_s_frequency()),
	d_xi(vars->dGETxi()),
	i_n_start_collecting_e(vars->iGETstart_collecting_e()),
	b_writeodata(vars->bGETwriteodata()),
	d_dt(vars->dGETdt())

{/*startvimfold*/
	this->vars = vars;
	walkers.init();
	initsim.init();
	int i_blocksize = 10000;

	i_numthreads = vars->iGETnumthreads();	
	pi_nnew_node = new unsigned int [i_numthreads];
	b_updates = false;
	i_nnew_node = 0;
	d_ecum_global = 0;
	d_projenergy_global = 0.0;
	d_projenergynow_node = 0.0;
	d_projenergynow_global = 0.0;
	ll_esamp_global = 0;
	ll_n0_global = 0;
	l_n0now_node = 0;//, i_nj = 0;
	l_n0now_global = 0;
	i_threadcount_node = 0;
	dS = vars->dGETdS();

}/*endvimfold*/

template<size_t N>
runSimulation<N>::~runSimulation()
{/*startvimfold*/
	delete [] pi_nnew_node;
}/*endvimfold*/

template<size_t N>
void runSimulation<N>::run()
{ /*startvimfold*/

	// init walkers
	setInitialState();
	walkers.loadBalance(); 

//	if (vars->bGETwriteodata() && (i_myrank==0) )
//		eOutStream::init(i_myrank, vars->sGETofpath());
	
	// Open filestream : blockingdata to file
	if (vars->bGETwriteodata() && (i_myrank==0) )
	{
		ostringstream ost;
		ost <<vars->sGETofpath()<<
			"mixed.dat";
		ofile.open(ost.str().c_str(), ios::out | ios::binary );
	}
	// Parallel section starts here.
#pragma omp parallel default(shared)
	{
		// Variables that are private to each thread
		const unsigned int i_thread = omp_get_thread_num();
#if 1
#pragma omp single
		cout << "threads (";
#pragma omp critical
		cout << i_thread << ",", flush(cout);
#pragma omp barrier

#pragma omp single
		cout << ") from task (" << i_myrank << ") sais hi!" << endl;
#endif
		int i_numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &i_numprocs);
		
		// Initiate the walkerPropagator objects
		walkerPropagator<N> o_walkerPropagator(*vars);
		o_walkerPropagator.init(
				i_thread,
				i_myrank);
		// FCIQMC Loop
#pragma omp barrier
		mainLoop(i_thread, o_walkerPropagator);
	} //end omp parallel

	if (i_myrank==0)
	{	
		cout << endl;
	}
	if ( (vars->bGETwriteodata()) && (i_myrank==0))
	{
		ofile.close();
	}
};/*endvimfold*/

template<size_t N>
void runSimulation<N>::mainLoop(const unsigned int i_thread, walkerPropagator<N> &o_walkerPropagator)
{/*startvimfold*/
	long** ppl_new_thread = (long**) matrix(vars->iGETmaxndets()/5, 2, sizeof(long)); // /5 OK??
	bitset<N>* pb_new_thread = new bitset<N>[vars->iGETmaxndets()/5];
	int j, i_nnew;
	unsigned int i_threadcount;

	for (int i_mainloop=0; i_mainloop<vars->iGETnum_loops(); i_mainloop++)
	{
		//if (l_numwalkers_global>i_limit_nw) break; 
		//#pragma omp barrier
		i_nnew = 0;
		// single FCIMC step (wdistthr : load balance on threads)
		double d_thread_projenergynow_node = 0;
		long l_thread_n0now = 0;
		o_walkerPropagator.singleStep(
				dS,
				&walkers.pb_slater[walkers.wdistthr[i_thread].first],
				pb_new_thread,
				&walkers.pl_population[walkers.wdistthr[i_thread].first],
				ppl_new_thread,
				i_nnew,
				d_thread_projenergynow_node,
				l_thread_n0now,
				walkers.wdistthr[i_thread].ndet,
				walkers.wdistthr[i_thread].nwfirst_thrd,
				walkers.wdistthr[i_thread].nwlast_thrd,
				labs(walkers.pl_population[walkers.wdistthr[i_thread].first]),
				labs(walkers.pl_population[walkers.wdistthr[i_thread].first]+walkers.wdistthr[i_thread].ndet-1)
				);
		// writing private variables to a global variable
#pragma omp atomic
		d_projenergynow_node += d_thread_projenergynow_node;
#pragma omp atomic
		l_n0now_node += l_thread_n0now;
#pragma omp atomic
		i_nnew_node += i_nnew;
		// move all the newly spawned walkers from all the threads to a shared array
#pragma omp critical (runSimulation)
		{
			i_threadcount = i_threadcount_node; ++i_threadcount_node;
			if (i_threadcount == 0)
				pi_nnew_node[i_threadcount] = i_nnew;
			else
				pi_nnew_node[i_threadcount] = pi_nnew_node[i_threadcount-1] + i_nnew;
		}
		// Store all elements on the shared arrays.
		j = 0;
		for (unsigned int ui = (i_threadcount==0) ? 0 : pi_nnew_node[i_threadcount-1]; 
				ui<pi_nnew_node[i_threadcount]; ui++)
		{
			walkers.ppl_new[ui][0] = ppl_new_thread[j][0]; //*(walkers.ppl_new+2*j)
			walkers.ppl_new[ui][1] = ppl_new_thread[j][1]; //*(walkers.ppl_new+2*j+1)
			walkers.pb_new[ui] = pb_new_thread[j]; //*(walkers.pb_new + j)
			++j;
		}
#pragma omp barrier

//#pragma omp single nowait // one thread writes energies to disc
//		{
//			if (i_myrank==0)
//				eOutStream::flushData();
//		}
#pragma omp single//the next thread manages sorting and distribution of walkers
		{
			walkers.SETnnew(i_nnew_node);
			walkers.reDistributeSortMerge();
			collectGlobalData(); 
			// if all walkers die out, reset to initial state!
			if (l_numwalkers_global==0) 
			{
				setInitialState();
			}
			walkers.findNumWalkers(); //FIXME NECC??
			walkers.loadBalance(); // load balance threads
			accumulateProjectedEnergy();
			writeMixedDataToFile();
			calcAndSaveShift(i_mainloop);
			printToScreen(i_mainloop);

			// Reset global vars
			l_n0now_node = 0;
			d_projenergynow_node = 0.0;
			i_threadcount_node = 0;
			i_nnew_node = 0; //XXX DELETE

		}
//#pragma omp single
//		{
//			if (i_myrank==0)
//				eOutStream::writeToScreen(
//						i_mainloop-1, vars->iGETnum_loops(), 
//						vars->iGETstart_collecting_e());
//		}
//#pragma omp single 
//		{
//			eOutStream::update();
//		}
		//END PRAGMA OMP SINGLE (Implicit barrier)
	} //END FCIQMC LOOP
	// Free memory 
	free_matrix((void**)ppl_new_thread);
	delete [] pb_new_thread;
}/*endvimfold*/	

template<size_t N>
void runSimulation<N>::setInitialState()
{
	unsigned int i_initnumwalkers = 1;
	if (i_myrank==i_nprocs-1)
	{
		walkers.initiateWalkers(
				i_initnumwalkers, initsim.bGETinitial_state());
	}
	else
	{
		walkers.initiateWalkers(0);
	}
	l_ndet_global = 1;
	l_numwalkers_global = i_initnumwalkers;
}

template<size_t N>
void runSimulation<N>::collectGlobalData()
{
	long l_numwalkers_node = walkers.findNumWalkers();
	long l_ndet_node = walkers.lGETndet();

	i_nnew_node = 0;

	l_lastnumwalkers_global = l_numwalkers_global;
	// mpi : gather i_numwalkers_node, ....
	MPI_Allreduce(&l_numwalkers_node, &l_numwalkers_global, 1, 
			MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&l_n0now_node, &l_n0now_global, 1, 
			MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&d_projenergynow_node, &d_projenergynow_global, 1, 
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&l_ndet_node, &l_ndet_global, 1, 
			MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
}

template<size_t N>
void runSimulation<N>::writeMixedDataToFile()
{
	if (b_writeodata)
	{
	   	//to find the corr. e. per particle in Ry
		double d_Stof = dS * 2./vars->iGETnumpart();
		double d_petof = d_projenergynow_global * 2./vars->iGETnumpart();

		//write to file
		ofile.write((char*)&d_Stof, sizeof (double));
		ofile.write((char*)&l_numwalkers_global, sizeof (long));
		ofile.write((char*)&l_ndet_global, sizeof (long));
		ofile.write((char*)&l_n0now_node, sizeof (long));
		ofile.write((char*)&d_petof, sizeof (double));
	}
}

template<size_t N>
void runSimulation<N>::calcAndSaveShift(const int i_mainloop)
{
	if (!b_updates)
	{
		if (l_numwalkers_global>i_limit_nw)
		{
			b_updates = true;
			//i_n_start_collecting_e = i_mainloop + i_n_start_collecting_e;
		}
	}
	else
	{
		if (i_mainloop%i_update_s_frequency==0)
		{
			dS = dS - d_xi /(static_cast<double>(i_update_s_frequency)*d_dt)
				* log(
						static_cast<double>(l_numwalkers_global)
						/ static_cast<double>(l_lastnumwalkers_global)
					 );
			d_ecum_global += dS; //TODO change name -_node
			ll_esamp_global++; //TODO change name - _node
		}
	}
}

template<size_t N>
void runSimulation<N>::printToScreen(const int i_mainloop)
{
	if ( (i_mainloop%50==0) && (i_myrank==0) ) //func printtoscreen
	{
		//write to screen
		cout << setprecision(8)
			<< "\rdS: " << dS
			<< " S_cum:" << d_ecum_global/static_cast<double>(ll_esamp_global)
			* 2./vars->iGETnumpart() //to find the corr. e. per particle in Ry
			<< " pe: " << //d_referenceenergy +
			d_projenergy_global/static_cast<double>(ll_n0_global)
			* 2./vars->iGETnumpart() //to find the corr. e. per particle in Ry
			<< " n0: " << l_n0now_global
			<< " N.Det's: " << l_ndet_global
			<< " N.W: " << l_numwalkers_global
			<< " L: " << i_mainloop << " of "
			<< vars->iGETnum_loops();
		fflush(stdout);
	}
	if ( i_mainloop < i_n_start_collecting_e )
	{
		d_ecum_global = d_projenergy_global = 0.0;
		ll_n0_global = ll_esamp_global = 0;
	}
}

template<size_t N>
void runSimulation<N>::accumulateProjectedEnergy()
{
	d_projenergy_global += d_projenergynow_global;
	ll_n0_global += l_n0now_global;
	// calc the curren t proj. energy
	if (l_n0now_global!=0)
	{
		d_projenergynow_global = d_projenergynow_global 
			/ static_cast<double>(l_n0now_global);
	}
	else 
	{
		d_projenergynow_global = 0;
	}
}

#endif

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
