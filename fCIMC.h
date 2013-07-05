#ifndef FCIMC_H
#define FCIMC_H

#include <bitset>
#include "qDotHamiltonianElement_OFCI.h"
#include "imputVars.h"

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

template <size_t INORB>
class walkerPropagator 
{
	private:
		
		int ir;
		int i_active_orbs;
		int inum_part;
		double dnum_partC2;
		double d_dt;
		unsigned int i_initiatorlimit;
		double d_pexone;
		double d_referenceenergy;
	
		STOC_BASE o_random;
		imputVars vars;
		qDotHamiltonianElement o_hamiltonianelem;	

		int *pi_n;
		int *pi_m;
		int *pi_ss;
		int* pi_mm;
		int* pi_nn;

		bitset<INORB*2> b_temp;
		bitset<INORB*2>* pbamc_up;
		bitset<INORB*2>* pbamc_dwn;
	
		int* pi_ss_reference;
		bitset<INORB*2> b_referenceslater;
		bitset<INORB*2> b_slater_active;
		bitset<INORB*2> b_frozen_orbs;
		bitset<INORB*2> b_empty;


	public:
		walkerPropagator(imputVars&);
		~walkerPropagator();
		void init(
			int*,
			int*,
			bitset<INORB*2>,
			bitset<INORB*2>,
			libGRIE*,
			int,
			int);
		void singleIteration(
			const double, 
			bitset<INORB*2>*,
			bitset<INORB*2>*,
			long*,
			long**,
			int&,
			double&,
			long&,
        	const long,
			const long,
			const long,
			const long,
			const long
			);
		double updateProjectedEnergy(
				const int, 
				double,
				long &,
				long &,
				const int);
		void initializepbamcArrays(
				const bool);
		void setInitialState(
				const std::string);
		int totalSpin(
				const int ) const;
		int totalAngularMomentum(
				const int ) const;
		double sampleProjector(
				const int, 
				int &, int &, 
				int &, int &);
		double sampleProjector(
				const int, 
				int &, int &);
		void tryDoubleExcitation( 
				bitset<INORB*2>*,
				long**,
				const double, int &, 
				const int, 
				const int, const int, 
				const int, const int, 
				const long,	
				const bool);
		void trySingleExcitation(
				bitset<INORB*2>*,
				long**, const double, int &, 
				const int, const int, const int,
				const long,	
			   	const bool);
		double transSign(
				const int ia1, const int ic1, 
				const int ia2, const int ic2, 
				const int iact_w);
};

#endif
