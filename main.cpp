#include <mpi.h>
#include <iostream>
#include "last_run_definitions.h"
#include "imputVars.h"
#include "runSimulation.h"
#include "runSimulation.cpp"

using std::cerr;
/*!
 *
 * Main
 *
 */
int main(int argc, char** argv)
{
	// initiate MPI
	int i_required = MPI_THREAD_SERIALIZED;
	int i_provided;
	int i_myrank;
	int i_nprocs;
	MPI_Comm MPI_COMM_WORLD;
	MPI_Init_thread(&argc, &argv, i_required, &i_provided);
	//MPI_Init(&argc, &argv);//, i_required, &i_provided);
	MPI_Comm_size(MPI_COMM_WORLD, &i_nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &i_myrank);
	if (i_required<i_provided)
	{
		cerr << "\nMPI warning: MPI implementation has no support for muitiple threads\n";
	}
	MPI_Barrier(MPI_COMM_WORLD); //XXX delete

	// init run time variables
	imputVars vars;
	vars.init();

	runSimulation<INUM_ORBITALS> ors(i_myrank, i_nprocs, &vars);
	ors.run();

	MPI_Finalize();
}
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
