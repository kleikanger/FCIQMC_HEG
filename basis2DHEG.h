#ifndef basis2DHEG_H
#define basis2DHEG_H

#include <cstdlib>
#include <bitset>
#include "basis.h"

using std::bitset;

/*!
 * Class to organize the single particle basis, sample connected 
 * determinants and other tasks related to the single particle basis.
 * @N is thetemplate parameter half the size of the sp basis N=M/2
 */
template<size_t N>
class basis2DHEG : public basis<N>
{
	private:
		double d_kc, d_kc2;
		double d_l, d_omg;
		int i_kc;
		int i_ms, i_ns;
		int *pi_k, *pi_l, *pi_m;
		int *pi_pprobs, *pi_aprobs;
		int *pi_pnumstates, *pi_anumstates;
		int *pi_m2i;
		bitset<N> *pb_asets, *pb_psets;
		
		/*!
		 * functor to overload lessthan in sort
		 */
		class LessThan
		{
			private:
				int* pi_mom;
			public:
				/*! 
				 * Constructor
				 * @pi_mom squared mom. array.
				 */
				LessThan (int* pi_mom) { this->pi_mom = pi_mom; }

				bool operator()(const int& lhs, const int& rhs)
				{
					if (pi_mom[lhs]<pi_mom[rhs]) 
						return true;
					else
						return false;
				}
		};

	public:
		/*!
		 * Constructor. Initiate all arrays and parameters.
		 * @d_kc The cutoff momentum.
		 */
		basis2DHEG(double d_kc);
		/*!
		 * Destructor.
		 */
		~basis2DHEG();
		/*!
		 * @a component of the angular momentum.
		 * @b component of the angular momentum
		 * @c component of the angular momentum
		 * Get the index of the orbital with angular momenta
		 * (a,b,c)
		 */
		int getOrb(const int a, const int b, const int c);
		/*!
		 * Sample a random pair of orbitals to be created i_c1,i_c2 
		 * with probability -(i_invprob)^-1 given that the orbitals 
		 * i_a1,i_a2 are annihilated and that i_a1 and i_a2 have parallel
		 * spins.Only pairs where the momenta and the 
		 * spin is conserved are sampled

		 *
		 * @i_c1 is changed to the index of the first created orbital.
		 * @i_c2 is changed to the index of the second created orbital.
		 * @i_invprob is changed to the inverse prob. of drawing i_c1,i_c2.
		 * @i_a1 is the first annihilated orbital.
		 * @i_a2 is the second annihilated orbital.
		 * @i_spin is the spin of the annihilated orbitals (i_a1%2)
		 * @i_dran is a random uniform in (0,1)
		 */
		void sampleInteractionParallelSpins(
				int &i_c1, int &i_c2, int &i_invprob, 
				int i_a1, int i_a2, int i_spin, double d_ran);
		/*!
		 * Sample a random pair of orbitals to be created i_c1,i_c2 
		 * with probability -(i_invprob)^-1 given that the orbitals 
		 * i_a1,i_a2 are annihilated and that i_a1 and i_a2 have 
		 * anti-parallel spins. Only pairs where the momenta and the 
		 * spin is conserved are sampled
		 *
		 * @i_c1 is changed to the index of the first created orbital.
		 * @i_c2 is changed to the index of the second created orbital.
		 * @i_invprob is changed to the inverse prob. of drawing i_c1,i_c2.
		 * @i_a1 is the first annihilated orbital.
		 * @i_a2 is the second annihilated orbital.
		 * @i_dran1 is a random uniform in (0,1)
		 * @i_dran2 is a random uniform in (0,1)
		 */
		void sampleInteractionAntiparallelSpins(
				int &i_c1, int &i_c2, int &i_invprob, 
				int i_a1, int i_a2, 
				double d_ran1, double d_ran2);
		/*!
		 * Returns bitset with all bits set that represents orbitals in the
		 * first i_R+1 full shells. f.ex: use to set frozen orbs or activeorbs.
		 * @i_R is the shellnumber i_R = shells -1.
		 */
		bitset<2*N> getFullShells(int i_R);
		/*!
		 * @pi_k_ are changed to pointer to x-comp mom vecs
		 * @pi_l_ are changed to pointer to y-comp mom vecs
		 * @pi_m_ are changed to pointer to z-comp mom vecs
		 */
		void getPointersToMomVecs(int* &pi_k_, int* &pi_l_, int* &pi_m_);

	private:
		/*!
		* Initialize the single particle basis with all spin orbitals 
		* with a momentum that has smaller magnitude than d_kc. The 
		* arrays with the momenta (pi_k, pi_l, pi_n) are initialized. 
		*/
		void initBasis();
		/*
		* Initiate bitsets with all pairs of s.p. orbitals with a given 
		* total momentum (a,b,c). bitsets pb_asets and pb_psets contains pairs 
		* of anti-parallel or parallel spin. pi_pprobs and pi_aprobs contains 
		* the inverse probability of drawing one particular pair and 
		* pi_pnumstates and pi_anumstates contains the number of set bits in
		* the corresponding bitsets.
		*/
		void initMBS();
		/*
		 * The int arrays pi_pnumstates and pi_anumstates contains 
		 * the number of set bits in the corresponding bitsets pb_asets 
		 * and pb_psets.
		 */
		void initNumOrbsPerM();
		/*
		 * Initiate array with the index of the orbital with mom
		 * (a,b,c). stored as the a+b*i_ns+b*i_ns^2 'th element. 
		 */
		void initMom2Indx();
};
#endif
