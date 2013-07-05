#ifndef HAMILTONIAN2DHEG_H
#define HAMILTONIAN2DHEG_H

class imputVars;
#include "hamiltonianElements.h"

/*!
 * This class calculates the Hamiltonian matrix elements of the 2D Homogenous Electron Gas.
 */
class hamiltonian2DHEG : public hamiltonianElements
{
	private:
		int i_numpart;
		double d_k1;
		double d_k2;
		double d_kk;
		double d_mu2;
		double d_l;
		double *pd_spe;
		int* pi_k;
		int* pi_l;
	public:
		/*!
		 * Initiate the class.
		 * @i_numpart the number of particles
		 * @i_norb the number of orbitals in the single particle basis
		 * @vars pointer to instance of the imputVars class
		 * @pi_k pointer to array with x-components of orbitals (i_norb/2 elems)  
		 * @pi_l pointer to array with y-components of orbitals (i_norb/2 elems)  
		 * @pi_m empty pointer in the 2D case 
		 * @d_r0 the Wigner-Seitz radius
		 */
		void init(
				int i_numpart,
				int i_norb,
				imputVars* vars,
				int* pi_k,
				int* pi_l,
				int* pi_m
				);
		/*!
		 * Get the diagonal matrix element <D_i|H|D_i>.
		 * @pi_s array with the indices of the occupied arrays in |D_i>
		 * @return the diagonal matrix element <D_i|H|D_i>.
		 */
		double getDiagonalElement(
			const int *pi_s) const;
		/*!
		 * Get the off diagonal matrix element <D_i|H{a_k^+ a_p}|D_i>.
		 * @i_sp the annihilated orbital
		 * @i_sk the created orbital
		 * @pi_s array with the indices of the occupied arrays in |D_i>
		 * @return the matrix element <D_i|H{a_k^+ a_p}|D_i>.
		 */
		double getSingleExcitationAmplitude(
			const int i_sp, const int i_sk, const int *pi_s) const;
		/*!
		 * Get the off diagonal matrix element <D_i|H{a_p^+ a_q^+ a_l a_k}|D_i>.
		 * @i_sk the first annihilated orbital
		 * @i_sl the second annihilated orbital
		 * @i_sp the first created orbital
		 * @i_sq the second created orbital
		 * @return the off diagonal matrix element <D_i|H{a_p^+ a_q^+ a_l a_k}|D_i>.
		 */
		double getDoubleExcitationAmplitude(
			const int iak, const int ial, const int iap, const int iaq) const;
};

#endif
