#ifndef basis_H
#define basis_H

#include <cstdlib>
#include <bitset>

using std::bitset;

/*!
 * Class to organize the single particle basis, sample connected 
 * determinants and other tasks related to the single particle basis.
 * @N is thetemplate parameter half the size of the sp basis N=M/2
 */
template<size_t N>
class basis
{
	public:
		/*!
		 * Constructor. Initiate all arrays and parameters.
		 * @d_kc The cutoff momentum.
		 */
//		init();
		/*!
		 * Destructor.
		 */
//		~basis();
		/*!
		 * @a component of the angular momentum.
		 * @b component of the angular momentum
		 * @c component of the angular momentum
		 * Get the index of the orbital with angular momenta
		 * (a,b,c)
		 */
		virtual int getOrb(const int a, const int b, const int c) =0;
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
		virtual void sampleInteractionParallelSpins(
				int &i_c1, int &i_c2, int &i_invprob, 
				int i_a1, int i_a2, int i_spin, double d_ran) =0;
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
		virtual void sampleInteractionAntiparallelSpins(
				int &i_c1, int &i_c2, int &i_invprob, 
				int i_a1, int i_a2, 
				double d_ran1, double d_ran2) =0;
		/*!
		 * Returns bitset with all bits set that represents orbitals in the
		 * first i_R+1 full shells. f.ex: use to set frozen orbs or activeorbs.
		 * @i_R is the shellnumber i_R = shells -1.
		 */
		virtual void getPointersToMomVecs(int* &pi_k_, int* &pi_l_, int* &pi_m_) =0;
};
#endif
