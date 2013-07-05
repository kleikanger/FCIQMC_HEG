#ifndef RUNSIMULATION_H
#define RUNSIMULATION_H

#include "imputVars.h"
#include "initSimulation.h"
#include "walkerContainerClass.h"
#include "walkerPropagator.h"
#include <fstream>



template <size_t N>
class runSimulation
{
	private:
		imputVars* vars;
		initSimulation<N> initsim;
		walkerContainerClass<N> walkers;
		std::ofstream ofile;
		
		unsigned int* pi_nnew_node;
		
		unsigned int i_threadcount_node;
		int i_numthreads;
		int i_nnew_node;
		const int i_nprocs;
		const int i_myrank;
		const int i_limit_nw;
		const int i_update_s_frequency;
		const int i_n_start_collecting_e;
		long l_n0now_node;
		long l_n0now_global;
		long l_ndet_global;
		long l_numwalkers_global;
		long l_lastnumwalkers_global;
		long long ll_esamp_global;
		long long ll_n0_global;
		bool b_updates;
		const bool b_writeodata;
		double d_ecum_global;
		double d_projenergy_global;
		double d_projenergynow_node;
		double d_projenergynow_global;
		double dS;
		const double d_dt;
		const double d_xi;

		void setInitialState();
		void collectGlobalData();
		void writeMixedDataToFile();
		void calcAndSaveShift(const int);
		void printToScreen(const int);
		void accumulateProjectedEnergy();
		void mainLoop(const unsigned int, walkerPropagator<N>&);

	public:
		runSimulation(int, int, imputVars*); 
		~runSimulation(); 
		void run();
};

#endif
