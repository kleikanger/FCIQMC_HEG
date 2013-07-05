#ifndef FCIMC_H
#define FCIMC_H

#include <bitset>
#include "hamiltonianElements.h" //TODO move to source
#include "imputVars.h" //TODO move to source
#include "basis.h" //TODO move to source

//class hamiltonianElements;
//class basis;
//class imputVars;

//
// Random number generators:
// Choose between Mersenne Twister and SFMT generator
//
#ifndef STOC_BASE
#if 1 
	#include "rng/randomc.h"
 	#define  STOC_BASE CRandomMersenne // define random number generator base class
#else
	#include "rng/sfmt.h"
	#define  STOC_BASE CRandomSFMT
#endif
#endif

//using std::ofstream;
//using std::string;
using std::bitset;

/*!
 * This class performs the kill/clone and spawn steps of a single excitation 
 * of the FCIMQC algorithm for a number of walkers.
 * @INORB the number of spin orbitals in the basis.
 */
template <size_t INORB>
class walkerPropagator 
{
	private:
		
		int i_active_orbs, inum_part;
		unsigned int i_initiatorlimit;
		double d_dt, d_pexone, d_referenceenergy, dnum_partC2, d_HF;
		STOC_BASE o_random;
		imputVars vars;
		hamiltonianElements* o_hamiltonianmatrix;
		basis<INORB/2>* o_basis;
		int *pi_ss, *pi_ss_reference;
		bitset<INORB> b_temp;
		bitset<INORB> b_referenceslater;
		bitset<INORB> b_slater_active;
		bitset<INORB> b_frozen_orbs;
		bitset<INORB> b_empty;

	public:
		/*!
		 * Constructor.
		 * @vars instance of the class imputVars contains runtime variables
		 */
		walkerPropagator(imputVars& vars);
		/*!
		 * Destructor.
		 */
		~walkerPropagator();
		/*!
		 * Initiate class.
		 * @i_thread_id number of the current MP thread.
		 * @i_myrank number of current MPI task. 
		 */
		void init(
				int i_thread_id,
				int i_myrank);
		/*!
		 * Do a single FCIMQC iteration.
		 * @dS the energy shift.
		 * @pbslater array of the occupied determinants (dets).
		 * @pbnewdets array to store the spawned (new) dets.
		 * @pi_occupancy array with the populations of the dets.
		 * @ppl_occ the populations of the newly spawned dets.
		 * @i_nnewdets are changed to the number of spawned dets.
		 * @d_projectedenergynow changed to the sum of the projected energies. 
		 * @i_n0now is changed to the number of dets on the reference det |D_0>.
		 * @i_ndetnow the number of dets distribution to this thread.
		 * @i_nw_first the number of walkers on the first det (this thread).
		 * @i_nw_last the number of walkers on the last det (this thread).
		 * @i_nwglob_first the  number of walkers on the first det (all threads).
		 * @i_nwglob_last the  number of walkers on the last det (all threads).
		 */
		void singleStep(
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
				);

	private:
		/*!
		 * Find <0|H|iactw>. Sample the projected (correlation) energy and add samples to 
		 * d_projectedenergynow and i_n0now. 
		 * @iactw the index of the "active" det.
		 * @d_eactw the energy of the active walker.
		 * @i_n0 changed to the sum of the number of samples of <0|H-E_HF|0>.
		 * @i_nj changed to the sum of the number of samples of <0|H|iactw>, iactw!=0.
		 * @i_niactw the (signed) number of walker on the active det.
		 * @return the summed proj. of the walkers on this determinant.
		 */
		double updateProjectedEnergy(
				const int iactw, 
				double d_eactw,
				long &i_n0,
				long &i_nj,
				const int i_niactw);
		/*
		 * Sample a connected determinant which is a double excitation of the 
		 * active determinant.
		 * @iact_w the index of the active walker
		 * @ia1 changed to the index of the 1. annihilated orbital.
		 * @ia2 changed to the index of the 2. annihilated orbital.
		 * @ic1 changed to the index of the 1. created orbital.
		 * @ic2 changed to the index of the 2. created orbital.
		 * @return the inverse sampling probability.
		 */
		double sampleProjector(
				const int iactw, 
				int &ia1, int &ia2, 
				int &ic1, int &ic2);
		/*
		 * Sample a connected determinant which is a single excitation of the 
		 * active determinant.
		 * @iact_w the index of the active walker
		 * @ia1 changed to the index of the 1. annihilated orbital.
		 * @ic1 changed to the index of the 1. created orbital.
		 * @return the inverse sampling probability.
		 */
		double sampleProjector(
				const int, 
				int &, int &);
		/*!
		 * Try to perform a double excitation.
		 * @pbnewdets
		 * @ppl_occ 
		 * @dp_gen
		 * @i_nnewdets
		 * @iact_w
		 * @ia1
		 * @ia2
		 * @ic1
		 * @ic2
		 * @i_numactw
		 * @b_sign_actw
		 */
		void tryDoubleExcitation(
				bitset<INORB>* pbnewdets,
				long** ppl_occ,
				const double d_pgen, int &i_nnewdets, 
				const int iact_w,  
				const int ia1, const int ia2, 
				const int ic1, const int ic2,
				const long i_numactw,	
				const bool b_sign_actw);
		/*!
		 * Try to perform a single excitation.
		 * @pbnewdets
		 * @ppl_occ 
		 * @dp_gen
		 * @i_nnewdets
		 * @iact_w
		 * @ia1
		 * @ic1
		 * @i_numactw
		 * @b_sign_actw
		 */
		void trySingleExcitation(
				bitset<INORB>* pbnewdets,
				long** ppl_occ,
				const double d_pgen, 
				int &i_nnewdets, 
				const int iact_w, 
				const int ia1, 
				const int ic1,
				const long i_numactw,	
				const bool b_sign_actw);
		/*!
		 * Find the sign of the determinant |D_j> = +/- a_c1 a_c2 a_a2 a_a1|D_i>
		 * @ia1 the index of the 1. annihilated spin orbital. 
		 * @ia2 the index of the 2. annihilated spin orbital. 
		 * @ic1 the index of the 1. excited spin orbital. 
		 * @ic1 the index of the 2. excited spin orbital. 
		 * @iact_w the index of the current determinant.
		 * @return the sign +/-1.
		 */
		double transSign(
				const int ia1, const int ic1, 
				const int ia2, const int ic2, 
				const int iact_w);
};

#endif
