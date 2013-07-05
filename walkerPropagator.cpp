#ifndef WALKERPROPAGATOR_CPP
#define WALKERPROPAGATOR_CPP

#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "walkerPropagator.h"
#include "last_run_definitions.h"

//instantiations of the abstract classes basis and hamiltonianElements
#include "hamiltonian2DHEG.h"
#include "hamiltonian3DHEG.h" 
#include "basis3DHEG.h" 
#include "basis3DHEG.cpp" //must be included since this is a template class 
#include "basis2DHEG.h" 
#include "basis2DHEG.cpp" //must be included since this is a template class 

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

using namespace std;

template <size_t INORB>
walkerPropagator<INORB>::walkerPropagator(imputVars &vars) : o_random(0)
{
	this->vars = vars;
}

template <size_t INORB>
walkerPropagator<INORB>::~walkerPropagator()
{
	delete [] pi_ss;
	delete [] pi_ss_reference;
}

template <size_t INORB>
void walkerPropagator<INORB>::init(
		int i_thread_id,
		int i_myrank)
{
	inum_part = vars.iGETnumpart();
	d_dt = vars.dGETdt();
	i_initiatorlimit = vars.iGETinitiatorlimit();
	d_pexone = vars.dGETpexone();
	dnum_partC2 = static_cast<double>((inum_part-1)*inum_part/2.);

	pi_ss = new int[inum_part];
	
	// initialize random number generator
	const int i_seed_a = vars.iGETranseed();
	const int i_seed_b = i_thread_id;
	const int i_seed_c = i_myrank;
	// array of seeds
	const int pi_seeds[] = 
		{
			i_seed_a,
			i_seed_b,
			i_seed_c
		};
	// generate random seeds
	o_random.RandomInitByArray(pi_seeds, 3); // Initialize RNG

	if (vars.sGETsystem()=="2DHEG")
	{
	 	o_hamiltonianmatrix = new hamiltonian2DHEG;
		o_basis = new basis2DHEG<INORB/2> (vars.dGETkc()); //XXX program this class
	}
	else if (vars.sGETsystem()=="3DHEG")
	{
		o_hamiltonianmatrix = new hamiltonian3DHEG;
		o_basis = new basis3DHEG<INORB/2> (vars.dGETkc());;
	}
	else
	{
	 	cout << "error: set the physical system (param: s_system)." ;
	}
	int *pi_k, *pi_l, *pi_m;
	o_basis->getPointersToMomVecs(pi_k, pi_l, pi_m);
	o_hamiltonianmatrix->init(inum_part, INORB, &vars, pi_k, pi_l, pi_m);

	// set the reference state to calculate the projected energy
	pi_ss_reference = new int[inum_part];
	b_referenceslater.reset(); // = o_basis->getFullShells(vars.iGETnumfullshells());
	for (int i=0;i<vars.iGETnumpart();i++)
		b_referenceslater.set(i);

	pi_ss_reference[0] = b_referenceslater._Find_first();
	for (int i=1;i<inum_part;i++)
	{
		pi_ss_reference[i] = b_referenceslater._Find_next(pi_ss_reference[i-1]);
	}
	
	d_HF = o_hamiltonianmatrix->getDiagonalElement(pi_ss_reference);
	if ( (i_myrank==0) && (i_thread_id==0) )
	cout << "\nHF energy: " << d_HF*2./vars.iGETnumpart() << "\n";

	// the following vars defines the CAS
	b_frozen_orbs = b_referenceslater; //TODO we only include the reference slater in the CAS
	b_empty.reset();
//	i_active_orbs = vars.iGETactive_orbs(); //TODO remove
   	
}

template <size_t INORB>
void walkerPropagator<INORB>::singleStep(
		double dS,
		bitset<INORB>* pbslater,
		bitset<INORB>* pbnewdets,
		long* pi_occupancy,
		long** ppl_occ,
		int &i_nnewdets,
		double &d_projectedenergynow,
		long &i_n0now,
		const long i_ndetnow,
		const long i_nw_first,
		const long i_nw_last,
		const long i_nwglob_first,
		const long i_nwglob_last
		)
{
	long i_njnow = 0;
	double d_pd;
	long i, j, i_loopw, i_numactw;
	int ia1,ia2,ic1,ic2, i_nspawn; //annihilated & created orbitals
	double d_pgen; //generational probability.
	bool b_sign_actw;

	dS += d_HF;

	//i_ndetnow populated determinants
	for	(int iact_w=0;iact_w<i_ndetnow;iact_w++)
	{
		//the active walker
		b_slater_active = pbslater[iact_w];

		if (iact_w==i_ndetnow-1) 
		{
			//num w.s
			i_loopw = abs(i_nw_last);
			//the sign of the active walker
			b_sign_actw = (i_nw_last>0);
			//the total unsigned population on this det
			i_numactw = i_nwglob_last;
		}
		else if (iact_w==0) 
		{
			i_loopw = abs(i_nw_first);
			b_sign_actw = (i_nw_first>0); 
			i_numactw = i_nwglob_first; 
		}
		else
		{	
			i_loopw = abs(pi_occupancy[iact_w]);
			i_numactw = i_loopw; 
			b_sign_actw = (pi_occupancy[iact_w]>0); 
		}

		//calculate the diagonal hamiltonian element
		i = b_slater_active._Find_first();
		pi_ss[0] = i;
		for (j=1;j<inum_part;j++)
		{
			i = b_slater_active._Find_next(i);
			pi_ss[j] = i;
		}
		d_pd = o_hamiltonianmatrix->getDiagonalElement(pi_ss);

		//accumulate the projected energy
		d_projectedenergynow += updateProjectedEnergy(
				iact_w, d_pd, i_n0now, 
				i_njnow, i_loopw * ((b_sign_actw) ? 1 : -1) );

		//calculate the killing/cloning probability
		d_pd = d_dt*(d_pd-dS);

		//loop over all walkers on the given determinant. loop_w is the number of walkers 
		//processed by this thread and can be different from pi_occupancy if the walkers 
		//on one determinant are split between two threads.
		for (; i_loopw>0; i_loopw--)
		{
			// TODO remove test, The probability to do a single ex is always 
			// 0 for the Homogenous Electron Gas
			if (false) // (o_random.Random())<d_pexone) //try to do a single excitation
			{
				//Find dp_gen. returns negative value the excitation failed
				d_pgen = sampleProjector(iact_w,ia1,ic1);
				if(d_pgen>0)
				{
					d_pgen /= d_pexone;
					trySingleExcitation(
							pbnewdets, ppl_occ, 
							d_pgen, i_nnewdets, iact_w, 
							ia1, ic1, i_numactw, b_sign_actw);
				}
			}
			else //try to do a double excitation.
			{
				//find p_gen. Returns a negative number if the exitation failed 
				d_pgen = sampleProjector(iact_w,ia1,ia2,ic1,ic2);
				if(d_pgen>0)
				{
					d_pgen /= (1.-d_pexone); // TODO not necc for HEG
					tryDoubleExcitation(
							pbnewdets, ppl_occ, 
							d_pgen, i_nnewdets, iact_w,  
							ia1, ia2, ic1, ic2, i_numactw, b_sign_actw);
				}
			}
			//Spawn or die on ii
			if (d_pd>0) //
			{
				if (d_pd>(o_random.Random()))
				{
					//kill walker(s?)
					i_nspawn = static_cast<int>(abs(d_pd)+1); 
#pragma omp atomic // FIXME ONLY NECC IF iact_w = 0 or n_det-1 !
					pi_occupancy[iact_w] +=	(b_sign_actw) ? -i_nspawn : i_nspawn;
				}
			}
			else if (fabs(d_pd)>(o_random.Random()))
			{
				//spawn walker(s?)
				i_nspawn = static_cast<int>(abs(d_pd)+1); 
#pragma omp atomic //FIXME ONLY NECC IF iact_w = 0 or n_det-1 !
				pi_occupancy[iact_w] +=	(b_sign_actw) ? i_nspawn : -i_nspawn;
			}
		}// END WHILELOOP OVER WALKERS
	}// END FORLOOP OVER DETERMINANT
}

template <size_t INORB>
double walkerPropagator<INORB>::updateProjectedEnergy(
		const int iactw,
		double d_eactw,
		long& i_n0, 
		long& i_nj, 
		const int i_niactw
		)
{
	double ret = 0.0;
	
	b_temp = (b_slater_active&~b_referenceslater);
	const int n = b_temp.count();

	if (n>2)
	{
		return 0.0;
	}
	else if (n==1)
	{
		const int a = b_temp._Find_first();
		b_temp = (b_referenceslater&~b_slater_active);
		const int b = b_temp._Find_first();
		//calculate single ex. ampl.
		ret = o_hamiltonianmatrix->getSingleExcitationAmplitude(b, a, pi_ss_reference);
	}
	else if (n==2)
	{
		const int a = b_temp._Find_first();
		const int b = b_temp._Find_next(a);
		b_temp = (b_referenceslater&~b_slater_active);
		const int c = b_temp._Find_first();
		const int d = b_temp._Find_next(c);
		//calculate double ex. ampl.
		ret = o_hamiltonianmatrix->getDoubleExcitationAmplitude(b, a, d, c);
		ret *= transSign(b, d, a, c, iactw);	
	}
	else if (n==0)
	{
		i_n0 += i_niactw;
		//ret the diagonal elem
		ret = d_eactw-d_HF;
		// cout << "\n" << d_eactw << "\n" ;
	}
	return ret * i_niactw;
}

template <size_t INORB>
double walkerPropagator<INORB>::transSign(
		const int ia1, 
		const int ic1, 
		const int ia2, 
		const int ic2, 
		const int iact_w)
{
	int ip = 0;
	int i;

	if (ic1<ic2)
	{
		i = ic1;
		while (true)
		{	
			i = b_slater_active._Find_next(i);
			if (i<ic2)
			{
				if ( (i!=ia1) && (i!=ia2) ) ++ip;
			}
			else break;
		}
	}
	else
	{
		i = ic2;
		while (true)
		{
			i = b_slater_active._Find_next(i);
			if (i<ic1)
			{
				if ( (i!=ia1) && (i!=ia2) ) ++ip;
			}
			else break;
		}
	}
	if (ia1<ia2)
	{
		i = ia1;
		while (true)
		{
			i = b_slater_active._Find_next(i);
			if (i<ia2) ++ip;
			else break;
		}
	}
	else
	{
		i = ia2;
		while (true)
		{
			i = b_slater_active._Find_next(i);
			if (i<ia1) ++ip;
			else break;
		}
	}

	return (ip%2==0) ? 1.0 : -1.0; //static_cast<double>( 1-2*(ip%2) );
}

template <size_t INORB>
void walkerPropagator<INORB>::trySingleExcitation(
		bitset<INORB>* pbnewdets,
		long** ppl_occ,
		const double d_pgen, 
		int &i_nnewdets, 
		const int iact_w, 
		const int ia1, 
		const int ic1,
		const long i_numactw,	
		const bool b_sign_actw)
{

	//Find off diagonal matrix element and the spawning amplitude
	double dp = o_hamiltonianmatrix->getSingleExcitationAmplitude(ia1, ic1, pi_ss);
	dp *= d_dt*d_pgen;
	const bool bsign = (dp<0);
	dp = fabs(dp);

	//Count how many new walkers to spawn
	int i_nspawn=0;
	if (dp>1)
	{
		while (dp>1)
		{
		   i_nspawn++;
		   dp-=1.0;
		}
	}
	if (dp>(o_random.Random())) 
		i_nspawn++;
	//spawn walkers
	if (i_nspawn>0)
	{
		// check if the walker is an initiator
		bool b_initiator = true;
		if (i_numactw<i_initiatorlimit)
		{
			// Check if the determinant is inside the cas space
			// test if all set bits in b_frozen_orbs are set in b_slater_active
			if ( (b_frozen_orbs&b_slater_active) != b_frozen_orbs)
			{
				b_initiator = false;
			}
//TODO we only use the determinant with the frozen orbs (= the reference determinant) as the CAS
//			// test if any particles are outside the active orbs
//			else if (b_slater_active._Find_next(i_active_orbs-1) < INORB)
//			{
//				b_initiator = false;
//			}
		}

		//if initiator spawn from initiator 
		if (b_initiator)
		{
			ppl_occ[i_nnewdets][0] = 
				-i_nspawn*(2.*(int)(bsign!=(b_sign_actw))-1);
			ppl_occ[i_nnewdets][1] = 0; 
		}
		//else spawn walker from non initiator
		else
		{
			ppl_occ[i_nnewdets][0] = 0; 
			ppl_occ[i_nnewdets][1] = 
				-i_nspawn*(2.*(int)(bsign!=(b_sign_actw))-1);
		}

		//spawn Walker
		pbnewdets[i_nnewdets] = b_slater_active;
		pbnewdets[i_nnewdets].RESET(ia1);
		pbnewdets[i_nnewdets].SET(ic1);

		i_nnewdets++;
	}
}

template <size_t INORB>
void walkerPropagator<INORB>::tryDoubleExcitation(
		bitset<INORB>* pbnewdets,
		long** ppl_occ,
		const double d_pgen, int &i_nnewdets, 
		const int iact_w,  
		const int ia1, const int ia2, 
		const int ic1, const int ic2,
		const long i_numactw,	
		const bool b_sign_actw)
{

	//Find off diagonal matrix element
	double dp = o_hamiltonianmatrix->getDoubleExcitationAmplitude(ia1, ia2, ic1, ic2);
	dp *= transSign(ia1, ic1, ia2, ic2, iact_w);	
	dp *= d_dt*d_pgen;
	
	const bool bsign = (dp<0);
	dp = fabs(dp);

	//Count how many new walkers to spawn
	int i_nspawn=0;
	if (dp>1)
	{
		while (dp>1)
		{
		   i_nspawn++;
		   dp-=1.0;
		}
	}
	if (dp>(o_random.Random())) 
		i_nspawn++;

	//spawn walkers
	if (i_nspawn>0)
	{
		// check if the walker is an initiator
		bool b_initiator = true;
		if (i_numactw<i_initiatorlimit)
		{
			// Check if the determinant is inside the cas space
			// test if all set bits in b_frozen_orbs are set in b_slater_active
			if ( (b_frozen_orbs&b_slater_active) != b_frozen_orbs)
			{
				b_initiator = false;
			}
//TODO we only use the (frozen orbs = the reference determinant) 
// = the lowest energy determinant as the cas space
//			// test if any particles are outside the active orbs
//			else if (b_slater_active._Find_next(i_active_orbs-1) < INORB)
//			{
//				b_initiator = false;
//			}
		}

		//if initiator spawn from initiator 
		if (b_initiator)
		{
			ppl_occ[i_nnewdets][0] = 
				-i_nspawn*(2.*(int)(bsign!=(b_sign_actw))-1);
			ppl_occ[i_nnewdets][1] = 0; 
		}
		//else spawn walker from non initiator
		else
		{
			ppl_occ[i_nnewdets][0] = 0; 
			ppl_occ[i_nnewdets][1] = 
				-i_nspawn*(2.*(int)(bsign!=(b_sign_actw))-1);
		}

		//spawn a new walker;
		pbnewdets[i_nnewdets] = b_slater_active;
		pbnewdets[i_nnewdets].RESET(ia1);
		pbnewdets[i_nnewdets].RESET(ia2);
		pbnewdets[i_nnewdets].SET(ic1);
		pbnewdets[i_nnewdets].SET(ic2);

		++i_nnewdets;
	}
}

template <size_t INORB> //XXX single ex not possible for HEG
double walkerPropagator<INORB>::sampleProjector(
		const int iact_w, 
		int &ia1, 
		int &ic1)
{
	cerr << "\nerror : zero probability for single excitations\n";
	return static_cast<double>(-1);
}

template <size_t INORB>
double walkerPropagator<INORB>::sampleProjector(
		const int iact_w, 
		int &ia1, 
		int &ia2, 
		int &ic1, 
		int &ic2)
{
	int //ia1,ia2,ic1,ic2,
		i,k;

	int i_inverseSProb;
//old method
#if 0	
	//pick random walkers to be anihilated
	ia1 = o_random.IRandomX(1,inum_part);
	ia2 = o_random.IRandomX(1,inum_part-1);
	//find the position of the ai1'th and ai2'th set bit	
	if (ia1>ia2)
	{	
		swap(ia1, ia2);
	}
	else //if (iia1<=iia2)
	{
		ia2++;
	}
	i = b_slater_active._Find_first();
	for (k=1;k<ia1;k++)
		i = b_slater_active._Find_next(i);
	ia1 = i;
	for (;k<ia2;k++)
		i = b_slater_active._Find_next(i);
	ia2 = i;

//this method should be much faster for large s.p. bases
#else

	//pick random walkers to be anihilated
	ia1 = o_random.IRandomX(1,inum_part);
	ia2 = o_random.IRandomX(1,inum_part-1);
	if (ia1>ia2)
	{	
		swap(ia1, ia2);
	}
	else //if (iia1<=iia2)
	{
		ia2++;
	}
	ia1 = pi_ss[ia1-1];
	ia2 = pi_ss[ia2-1];
#endif
	
//	const int i_spin1 = ia1%2;
//	const int i_spin2 = ia2%2;

	if (ia1%2==ia2%2)
	{
		const double d_ran1 = o_random.Random();
		o_basis->sampleInteractionParallelSpins
			(ic1, ic2, i_inverseSProb, ia1, ia2, ia1%2, d_ran1);
	}
	else
	{
		const double d_ran1 = o_random.Random();
		const double d_ran2 = o_random.Random();
		o_basis->sampleInteractionAntiparallelSpins
			(ic1, ic2, i_inverseSProb, ia1, ia2, d_ran1, d_ran2);
	}

	//zero prob if either ic1 or ic2 was previously occupied
	//b_temp.reset(); //TODO faster way??
	//b_temp.set(ic1);
	//b_temp.set(ic2);
	//if ( (b_slater_active&b_temp) != b_empty)
	if (b_slater_active[ic1]||b_slater_active[ic2])
	{
		return -1.0;
	}
	if (ic1>ic2)
	{
		const int i_temp = ic1;
		ic1 = ic2;
		ic2 = i_temp;
		//icc2 = ic1;
		//icc1 = ic2;
	}
//	else 
//	{
//		icc1 = ic1;
//		icc2 = ic2;
//	}
	// TODO not possible	
//	if (ia1>ia2)
//	{
//		iaa2 = ia1;
//		iaa1 = ia2;
//	}
//	else
//	{
//		iaa1 = ia1;
//		iaa2 = ia2;
//	}
	return dnum_partC2 * static_cast<double>(i_inverseSProb);

}

#endif
