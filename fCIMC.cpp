#include <bitset>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "fCIMC.h"
#include "last_run_definitions.h"
#include <exception>

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
/*
 *
 * Constructor
 *
 */
template <size_t INORB>
fCIMC<INORB>::fCIMC(imputVars &vars) : o_random(0)
{/*startvimfold*/
	this->vars = vars;
}/*endvimfold*/
/*
 *
 * Destructor
 *
 */
template <size_t INORB>
fCIMC<INORB>::~fCIMC()
{/*startvimfold*/
	delete [] pbamc_up;
	delete [] pbamc_dwn;
	delete [] pi_m;
	delete [] pi_n;
	delete [] pi_nn;
	delete [] pi_mm;	
	delete [] pi_ss;
	delete [] pi_ss_reference;
}/*endvimfold*/
/*
 *
 * init
 *
 */
template <size_t INORB>
void fCIMC<INORB>::init(
		int *pi_mARG,
		int *pi_nARG,
		bitset<INORB*2> b_referenceslaterARG,
		bitset<INORB*2> b_frozen_orbsARG,
		libGRIE* ogrie,
		int i_thread_id,
		int i_myrank)
{/*startvimfold*/
	ir = vars.iGETr();	
	inum_part = vars.iGETnumpart();
	d_dt = vars.dGETdt();
	i_initiatorlimit = vars.iGETinitiatorlimit();

	//Set the limit for the initiators in iFCIMC

	//
	// each object holds a copy of this small array to reduce 
	// overhead using openmpi
	//
	pi_m = new int[INORB];
	pi_n = new int[INORB];
	for (unsigned int i=0; i<INORB; i++)
		pi_m[i] = pi_mARG[i];
	for (unsigned int i=0; i<INORB; i++)
		pi_n[i] = pi_nARG[i];
	
	pi_mm = new int[inum_part];
	pi_nn = new int[inum_part];
	pi_ss = new int[inum_part];
	
	pbamc_up = new bitset<INORB*2>[ir*2+1];
	pbamc_dwn = new bitset<INORB*2>[ir*2+1];
	
	//the number of ways to draw a random occupied orbital
	dnum_partC2 = static_cast<double>((inum_part-1)*inum_part/2.);
	
	//Set the probability to do a single excitation instead if a double
	d_pexone = .90;
	
	//
	// initialize random number generator
	//
	int i_seed_a = vars.iGETranseed();
	int i_seed_b = i_thread_id+1;
	int i_seed_c = i_myrank+1;
	//int i_seed_c = i_nodenr;
	// array of seeds
	const int pi_seeds[] = 
		{
			i_seed_a,
			i_seed_b,
			i_seed_c
		};
	// generate random seeds
	o_random.RandomInitByArray(pi_seeds, 3); // Initialize RNG

	//
	// This is as array of bitsets with all determinants with a given angular momentum.
	// It is used to optimize the function sampleProjector()
	//
	initializepbamcArrays(vars.bGETenergy_cut());

	//Initiate the class returning the hamiltonian elements
	o_hamiltonianelem.init(inum_part, ogrie); 

	// Set the reference state to calculate the projected energy.
	// TODO make better method. Should automatically find 
	// the sd with the largest population
	// b_referenceslater.RESET();
	pi_ss_reference = new int[inum_part];
	b_referenceslater = b_referenceslaterARG;
	pi_ss_reference[0] = b_referenceslater._Find_first();
	for (int i=1;i<inum_part;i++)
	{
		pi_ss_reference[i] = b_referenceslater._Find_next(pi_ss_reference[i-1]);
	}
	
	// the following vars defines the CAS
	b_frozen_orbs = b_frozen_orbsARG;
	b_empty.reset();
	i_active_orbs = vars.iGETactive_orbs();

}/*endvimfold*/
/*
 *
 * singleIteration
 *
 */
template <size_t INORB>
void fCIMC<INORB>::singleIteration(
		double dS,
		bitset<INORB*2>* pbslater,
		bitset<INORB*2>* pbnewdets,
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
{/*startvimfold*/

	//reset some of the class variables
	//pbslater = pbslaterARG;	
	//pi_occupancy = pi_occupancyARG;

	long i_njnow = 0;

	double d_pd;
	long i, j, i_loopw, i_numactw;
	int ia1,ia2,ic1,ic2, i_nspawn; //annihilated & created orbitals

	double d_pgen; //generational probability.
	bool b_sign_actw;

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
		pi_mm[0] = pi_m[i/2];
		pi_nn[0] = pi_n[i/2];
		pi_ss[0] = i;
		for (j=1;j<inum_part;j++)
		{
			i = b_slater_active._Find_next(i);
			pi_mm[j] = pi_m[i/2];
			pi_nn[j] = pi_n[i/2];
			pi_ss[j] = i;
		}
		d_pd = o_hamiltonianelem.getHii(pi_mm, pi_nn, pi_ss);

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
			if (( (o_random.Random())<d_pexone)) //try to do a single excitation
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
					d_pgen /= (1.-d_pexone);
					tryDoubleExcitation(
							pbnewdets, ppl_occ, 
							d_pgen, i_nnewdets, iact_w,  
							ia1, ia2, ic1, ic2, i_numactw, b_sign_actw);
				}
			}
			//Spawn or die on ii
			if (d_pd>0)
			{
				//cerr << "c" ;
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
				//cerr << "d" ;
				//spawn walker(s?)
				i_nspawn = static_cast<int>(abs(d_pd)+1); 
#pragma omp atomic //FIXME ONLY NECC IF iact_w = 0 or n_det-1 !
				pi_occupancy[iact_w] +=	(b_sign_actw) ? i_nspawn : -i_nspawn;
			}
		}// END WHILELOOP OVER WALKERS
	}// END FORLOOP OVER DETERMINANT
}/*endvimfold*/
/*
 * Find <0|H|iactw>
 *
 * XXX TESTING XXX
 * XXX NB NOT NECC TO COLLECT i_nj XXX
 *
 */
template <size_t INORB>
double fCIMC<INORB>::updateProjectedEnergy(
		const int iactw,
		double d_eactw,
		long& i_n0, 
		long& i_nj, 
		const int i_niactw
		)
{/*startvimfold*/
	double ret = 0.0;
	
	//b_temp = (pbslater[iactw]&~b_referenceslater);
	b_temp = (b_slater_active&~b_referenceslater);
	const int n = b_temp.count();

	if (n>2)
	{
		return 0.0;
	}
	else if (n==1)
	{
		//number of samples +
		//i_nj += abs(i_niactw);
		
		//XXX Not optimal. Find smarter bitoperations 
		const int a = b_temp._Find_first();
		//b_temp = (b_referenceslater&~pbslater[iactw]);
		b_temp = (b_referenceslater&~b_slater_active);
		const int b = b_temp._Find_first();
		//calculate single ex. ampl.
		//a<b??
		ret = o_hamiltonianelem.getSingleExcitationAmplitude(b, a, pi_ss_reference);
	}
	else if (n==2)
	{
		//number of samples +
		//i_nj += abs(i_niactw);
		
		//XXX Not optimal. Find smarter bitoperations 
		const int a = b_temp._Find_first();
		const int b = b_temp._Find_next(a);
		//b_temp = (b_referenceslater&~pbslater[iactw]);
		b_temp = (b_referenceslater&~b_slater_active);
		const int c = b_temp._Find_first();
		const int d = b_temp._Find_next(c);

		//a<b? ?
		ret = o_hamiltonianelem.getDoubleExcitationAmplitude(b, a, d, c);
		ret *= transSign(b, d, a, c, iactw);	
	}
	else if (n==0)
	{
		//number of samples +i_n0 += abs(i_niactw);
		i_n0 += i_niactw;// abs(i_niactw);
		//cout << d_referenceenergy << endl;
		ret = d_eactw;
	}
	return ret * static_cast<double>(i_niactw);
}/*endvimfold*/
/*
 *
 * transSign
 *
 * Find the sign of the transition. Count the number of occupied states (ip) 
 * between the two created(new) states. Return (-1)^ip
 *
 * TODO : This method can probably be speeded up a bit.
 * Eg. if ia1 in(ic1,ic2) then ip++ (fewer if tests for many particles)
 */
template <size_t INORB>
double fCIMC<INORB>::transSign(
		const int ia1, 
		const int ic1, 
		const int ia2, 
		const int ic2, 
		const int iact_w)
{/*startvimfold*/
	int ip = 0;
	int i;

	if (ic1<ic2)
	{
		i = ic1;
		while (true)
		//for (i=ic1+1; i<ic2; i++)
		{	
			i = b_slater_active._Find_next(i);
			if (i<ic2)
			{
				if ( (i!=ia1) && (i!=ia2) ) ++ip;
			}
			else break;
		//	if ( (i!=ia1) && (i!=ia2) )
		//		{
		//			if (b_slater_active.TEST(i)) ++ip;
		//		}
		}
	}
	else
	{
		i = ic2;
		while (true)
		//for (i=ic2+1; i<ic1; i++)
		{
			i = b_slater_active._Find_next(i);
			if (i<ic1)
			{
				if ( (i!=ia1) && (i!=ia2) ) ++ip;
			}
			else break;
		//	if ( (i!=ia1) && (i!=ia2) )
		//	{
		//		if (b_slater_active.TEST(i)) ++ip;
		//	}
		}
	}
	if (ia1<ia2)
	{
		i = ia1;
		while (true)
	//	for (i=ia1+1; i<ia2; i++)
		{
	//		if (b_slater_active.TEST(i)) ++ip;

			i = b_slater_active._Find_next(i);
			if (i<ia2) ++ip;
			else break;
		}
	}
	else
	{
		i = ia2;
		while (true)
		//for (i=ia2+1; i<ia1; i++)
		{
		//	if (b_slater_active.TEST(i)) ++ip;
		
			i = b_slater_active._Find_next(i);
			if (i<ia1) ++ip;
			else break;
		}
	}

	return (ip%2==0) ? 1.0 : -1.0; //static_cast<double>( 1-2*(ip%2) );
}
/*endvimfold*/
/*
   Try to do a single excitation 
   d_pgen - generation probability
   i_ndet - number of determinants
   iact_w - active walker
   ia1,ic1 - annihi. + created orbs
   pb_sign_w[i] - signs of walker 1 = +
   returns number of walkers spawned on new determinants
   */
template <size_t INORB>
void fCIMC<INORB>::trySingleExcitation(
		bitset<INORB*2>* pbnewdets,
		long** ppl_occ,
		const double d_pgen, 
		int &i_nnewdets, 
		const int iact_w, 
		const int ia1, 
		const int ic1,
		const long i_numactw,	
		const bool b_sign_actw)
{/*startvimfold*/

	//Find off diagonal matrix element and the spawning amplitude
	double dp = o_hamiltonianelem.getSingleExcitationAmplitude(ia1, ic1, pi_ss);
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
		// i-FCI-QMC step
		//if the determinant is unpopulated, spawn on new determinant only if 
		//the current determinant has a large enough population.
		
	/*	// check if the walker is an initiator 
		if (i_numactw>=i_initiatorlimit)
		{
			ppl_occ[i_nnewdets][0] = 
				-i_nspawn*(2.*(int)(bsign!=(b_sign_actw))-1);
			ppl_occ[i_nnewdets][1] = 0; 
		}
		else
		{
			// Check if the determinant is inside the cas space
			bool b_incas = true;
			// test if all set bits in b_frozen_orbs are set in b_slater_active
			if ( (b_frozen_orbs&b_slater_active) != b_frozen_orbs )
			{
				b_incas = false;
			}
			// test if any particles are outside the active orbs
			else if (b_slater_active._Find_next(i_active_orbs-1) < INORB*2)
			{
				b_incas = false;
			}
			//if walker is in CAS, spawn walker from initiator 
			if (b_incas)
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
		}
		*/
		
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
			// test if any particles are outside the active orbs
			else if (b_slater_active._Find_next(i_active_orbs-1) < INORB*2)
			{
				b_initiator = false;
			}
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
}/*endvimfold*/
/*
   Try to do a double excitation 
   d_pgen - generation probability
   i_ndet - number of determinants
   iact_w - active walker
   ia1,ia2,ic1,ic2 - annihi. + created orbs
   pb_sign_w[i] - array of signs 1 = +
   returns number of walkers spawned on new determinants
   */
template <size_t INORB>
void fCIMC<INORB>::tryDoubleExcitation(
		bitset<INORB*2>* pbnewdets,
		long** ppl_occ,
		const double d_pgen, int &i_nnewdets, 
		const int iact_w,  
		const int ia1, const int ia2, 
		const int ic1, const int ic2,
		const long i_numactw,	
		const bool b_sign_actw)
{//startvimfold

	//Find off diagonal matrix element
	double dp = o_hamiltonianelem.getDoubleExcitationAmplitude(ia1, ia2, ic1, ic2);
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
			// test if any particles are outside the active orbs
			else if (b_slater_active._Find_next(i_active_orbs-1) < INORB*2)
			{
				b_initiator = false;
			}
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
}//endvimfold
/*
 *
 * sampleProjector
 *
 * find exited orbital.
 * double excitation.
 * iact_w - active walker
 * ia1 - anihilated walker
 * ic1 - created walker
 * returns p_gen
 *
 */
template <size_t INORB>
double fCIMC<INORB>::sampleProjector(
		const int iact_w, 
		int &ia1, 
		int &ic1)
{/*startvimfold*/
	int i,j;
	int ima;
	int idim_S5; 
	bitset<INORB*2> b_slatertemp;

	//pick random particle to be annihilated
	//and store it as ai1	
	ia1 = o_random.IRandomX(1,inum_part);
	
//	i=-1;
//	while (ia1>0)
//	{
//		i++;
//		if (b_slater_active[i])
//		{
//			ia1--;
//			//if (ia1==0)
//			//	break;
//		}
//	}
//	ia1=i;
	
	//find the position of the i'th set bit	
	i = b_slater_active._Find_first();
	for (j=1;j<ia1;j++)
		i = b_slater_active._Find_next(i);
	ia1 = i;

	/*
	//reset b_slatertemp and find s2
	b_slatertemp.RESET();
	*/
	
	//find the angular momentum of the new state
	ima = pi_m[ia1/2];
	//the spin must be the same as the annihil. orb
	
	
	/*
	if (ia1%2) //ic1 is spin down (1)
		for (i=0;i<INORB;i++)
			{
			if ( (pi_m[i]==ima) )
				b_slatertemp.SET(2*i+1);
			}
	else //ic1 is spin up (0)
		for (i=0;i<INORB;i++)
			{
			if ( (pi_m[i]==ima) )
				b_slatertemp.SET(2*i);
			}
	 */
	if (ia1%2) //ic1 is spin down (1)
		b_slatertemp = pbamc_dwn[ima+ir]&~b_slater_active;
	else //ic1 is spin up (0)
		b_slatertemp = pbamc_up[ima+ir]&~b_slater_active;
	
	/*
	//remove states that is part of the old slater determinant	
	b_slatertemp = b_slatertemp&~b_slater_active; //~returns the inverse
	//find dimension of S5
	*/

	idim_S5 = b_slatertemp.count();
	// Return -1 if dim  S2 = 0
	if (idim_S5==0) {return -1;}

	//pick random from S2
	ic1 = //static_cast<int>((o_random.Random())*idim_S5)+1;//
		o_random.IRandomX(1,idim_S5);
//	i=-1;
//	while (ic1>0)
//	{
//		i++;
//		if (b_slatertemp[i])
//			ic1--;
//	}
//	ic1=i;
	
	//find the position of the i'th set bit	
	i = b_slatertemp._Find_first();
	for (j=1;j<ic1;j++)
		i = b_slatertemp._Find_next(i);
	ic1 = i;
	
	return static_cast<double>(inum_part*idim_S5);
}/*endvimfold*/
/*
 *
 * sampleProjector
 *
 * suggest exited orbital.
 * double excitation.
 * iact_w - active walker
 * ia1 - first annihilated walker
 * ia2 - second annihilated walker
 * ic1 - first created walker
 * ic2 - second created walker
 * returns p_gen
 *
 *
   */
template <size_t INORB>
double fCIMC<INORB>::sampleProjector(
		const int iact_w, 
		int &iaa1, 
		int &iaa2, 
		int &icc1, 
		int &icc2)
{	/*startvimfold*/
	//active walker (loop through all)
	bool ba1, ba2;
	int i,j,k, iang_mom1, iang_mom2, ima, imb;
	int idim_S1, idim_S2, idim_S3, idim_S4;
	bitset<INORB*2> b_slatertemp;
	int ia1,ia2,ic1,ic2;

	//pick random walkers to be anihilated
	ia1 = o_random.IRandomX(1,inum_part);
	ia2 = o_random.IRandomX(1,inum_part-1);
	
//	i=-1;
//	while (ia1>0)
//	{
//		i++;
//		if (b_slater_active[i])
//			ia1--;
//	}
//	ia1=i; 
	
	//find the position of the ai1'th and ai2'th set bit	
	if (ia1>ia2)
	{	
		swap(ia1, ia2);
	}
	else //if (iia1<iia2)
	{
		ia2++;
	}
	i = b_slater_active._Find_first();
	for (k=1;k<ia1;k++)
		i = b_slater_active._Find_next(i);
	//Actually not all draws will give S1>0 here! eg two ang mom both = ir 
	ia1 = i;
	for (;k<ia2;k++)
		i = b_slater_active._Find_next(i);
	ia2 = i;

//	i=-1;
//	while ((ia2>0)) 
//	{
//		i++;
//		if (b_slater_active[i])
//			if (i!=ia1)
//				ia2--;
//	}
//	ia2=i;

	//find spins
	ba1=ia1%2; //1=spin up
	ba2=ia2%2;
	//find angular momenta
	iang_mom1=pi_m[ia1/2]; //integer division 
	iang_mom2=pi_m[ia2/2];
	
	//Find ima,imb (min/max spin of new states) using formula Eq. (...).
	//Removing most cases of S1=0 and improving the importance sampling.
	if (ba1!=ba2)
	{
		ima=iang_mom1+iang_mom2-ir;
		if (ima<-ir) ima=-ir;
		imb=iang_mom1+iang_mom2+ir;
		if (imb>ir) imb=ir;

		//checking if S1=0.
	   	//( speeding up the code ? )	
		//(More cases leading to dim(S1)=0 possible!?)
//		if ((ima==imb) && (abs(ima)==ir))
//			//check if immax / immin is occupied
//			if (ia1/2==i_maxm 
//					|| ia2/2==i_maxm
//					|| ia1/2==i_minm
//					|| ia2/2==i_minm) //integer division
//			{
//				return -1;
//			}
	}
	else //removing one case of dim(S1)=0 improving the importance sampling
	{
		ima=iang_mom1+iang_mom2-ir;
		if (ima<=-ir) ima=-ir+1;
		imb=iang_mom1+iang_mom2+ir;
		if (imb>=ir) imb=ir-1;
	}
	//calculate S1
	b_slatertemp.reset();
	if (ba1==ba2)//the spins of the new particles are also equal
		if (ba1) //(spin up) new particle has spin down
		{
			for (i=0;i<INORB;i++)
			{
				if ( (pi_m[i]>=ima) && (pi_m[i]<=imb) )
					b_slatertemp.SET(2*i+1);
			}
		} 
		else //new particle spin down
		{
			for (i=0;i<INORB;i++)
			{
				if ( (pi_m[i]>=ima) && (pi_m[i]<=imb) )
					b_slatertemp.SET(2*i);
			}
		}
	else //both spins allowed
	{
		for (i=0;i<INORB;i++)
			if ( (pi_m[i]>=ima) && (pi_m[i]<=imb) )
			{
				b_slatertemp.SET(2*i);
				b_slatertemp.SET(2*i+1);
			}
	}
	//remove states that is part of the old slater determinant	
	b_slatertemp=b_slatertemp&~b_slater_active; //~returns the inverse
	//find dim(S1), return -1 if S1 empty
	idim_S1=b_slatertemp.count();
	if (idim_S1==0) return -1;
	//Pick a random orbital from S1
	ic1 = //static_cast<int>((o_random.Random())*idim_S1)+1; //
		o_random.IRandomX(1,idim_S1);
	
//	i=-1;
//	while (ic1>0)
//	{
//		i++;
//		if (b_slatertemp[i])
//			ic1--;
//	}
//	ic1=i;

	//fint the index of the ic1'th orbital/set bit	
	i = b_slatertemp._Find_first();
	for (j=1;j<ic1;j++)
		i = b_slatertemp._Find_next(i);
	ic1 = i;

	//reset b_slatertemp and find s2
	b_slatertemp.reset();
	//find the angular momentum of the new state
	ima = iang_mom1+iang_mom2-pi_m[ic1/2];
	//find ic2 and dim (S2)
	if (ba1==ba2) //then ic2 spin must be the same as of ba1
	{
		if (ba1) //ic1 is spin down (1)
			b_slatertemp = pbamc_dwn[ima+ir]&~b_slater_active;
		else //ic1 is spin up (0)
			b_slatertemp = pbamc_up[ima+ir]&~b_slater_active;
			/*for (i=0;i<inorb;i++)
				{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i+1);
				}
		else //ic1 is spin up (0)
			for (i=0;i<inorb;i++)
				{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i);
				}*/
	}
	else //then ic2 spin must be opposite of ic1
	{
		if (ic1%2==0)
			b_slatertemp = pbamc_dwn[ima+ir]&~b_slater_active;
		else //ic1 is spin up (0)
			b_slatertemp = pbamc_up[ima+ir]&~b_slater_active;
		/*if (ic1%2==0)
		{
			for (i=0;i<inorb;i++)
			{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i+1);
			}
		}
		else
		{
			for (i=0;i<inorb;i++)
			{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i);
			}
		}*/
	}
	/*
	//remove states that is part of the old slater determinant	
	b_slatertemp=b_slatertemp&~b_slater_active; //~returns the inverse
	*/
	b_slatertemp.RESET(ic1); //already picked
	//find dimension	
	idim_S2 = b_slatertemp.count();

	// Return -2 if dim  S2 = 0
	if (idim_S2==0) {return -2;}
	
	//pick random from S2
	ic2 = //static_cast<int>((o_random.Random())*idim_S2)+1; //
		o_random.IRandomX(1,idim_S2);

//	i=-1;
//	while (ic2>0)
//	{
//		i++;
//		if (b_slatertemp[i])
//			ic2--;
//	}
//	ic2=i;

	//Find the ic2'th set bit	
	i = b_slatertemp._Find_first();
	for (j=1;j<ic2;j++)
		i = b_slatertemp._Find_next(i);
	ic2 = i;
	//now, find the frequency that ic2 is picked first and ic1 second
	idim_S3=idim_S1; //XXX always true?
	//reset b_slatertemp and use to find s4
	/*
	b_slatertemp.RESET();
	*/
	//find angular momentum of new state
	ima = iang_mom1+iang_mom2-pi_m[ic2/2];
	//find dim (S4)
	if (ba1==ba2) //then ic1 spin must be opposite of ba1
	{
		if (ba1) //ic2 is spin up (0)
			b_slatertemp = pbamc_dwn[ima+ir]&~b_slater_active;
		else //ic1 is spin up (0)
			b_slatertemp = pbamc_up[ima+ir]&~b_slater_active;
		/*	for (i=0;i<inorb;i++)
				{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i+1);
				}
		else //ic1 is spin down (1)
			for (i=0;i<inorb;i++)
				{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i);
				}*/
	}
	else //then ic1 spin must be opposite of ic2
	{
		if (ic2%2==0)
			b_slatertemp = pbamc_dwn[ima+ir]&~b_slater_active;
		else //ic1 is spin up (0)
			b_slatertemp = pbamc_up[ima+ir]&~b_slater_active;
		/*if (ic2%2==0)
		{
			for (i=0;i<inorb;i++)
			{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i+1);
			}
		}
		else
		{
			for (i=0;i<inorb;i++)
			{
				if ( (pi_m[i]==ima) )
					b_slatertemp.SET(2*i);
			}
		}*/
	}
	/*
	//remove states that is part of the old slater determinant	
	b_slatertemp=b_slatertemp&~b_slater_active; //~returns the inverse
	*/
	b_slatertemp.RESET(ic2); //already picked

	//find dimension	
	idim_S4=b_slatertemp.count();
	//if (idim_S4==0) cout << " !#!!"  << endl; should not be possible!
	
//XXX necc??	
	if (ia1>ia2)
	{
		iaa2=ia1;
		iaa1=ia2;
	}
	else
	{
		iaa1=ia1;
		iaa2=ia2;
	}
	if (ic1>ic2)
	{
		icc2=ic1;
		icc1=ic2;
	}
	else
	{
		icc1=ic1;
		icc2=ic2;
	}

	return dnum_partC2 * 
	   	static_cast<double>(idim_S1*idim_S2*idim_S3*idim_S4)
		/static_cast<double>(idim_S1*idim_S2 + idim_S3*idim_S4); 

}/*endvimfold*/
/*
 *
 * Init bitsets pbamc___[i] with all states with an angulat momentum i set.
 * pbamc_up  : only spin up states.
 * pbamc_dwn : only spin down states. 
 *
 * These arrays are used to optimize the functions sampleProjector(...) 
 *
 * TODO : b_energy_cut not needed? remove.
 *
 */
template <size_t INORB>
void fCIMC<INORB>::initializepbamcArrays(const bool b_energy_cut)
{/*startvimfold*/
	int i;
	unsigned int j;

	//set all to 0
	for (i=0; i<ir*2+1; i++)
		pbamc_up[i].reset();
	for (i=0; i<ir*2+1; i++)
		pbamc_dwn[i].reset();

	//init pbamc_up, pbamc_dwn arrays m-ir
	for (i=-ir; i<=ir; i++)
		for (j=0; j<INORB; j++)
			if (pi_m[j]==i)
			{
				pbamc_up[i+ir].SET(2*j);
				pbamc_dwn[i+ir].SET(2*j+1);
			}
}/*endvimfold*/


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=startvimfold,endvimfold
