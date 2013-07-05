#ifndef HAMILTONIANELEMENT_H
#define HAMILTONIANELEMENT_H

class imputVars;

/*!
 * This abstract class calculates the Hamiltonian matrix elements.
 */
class hamiltonianElements
{
	public:
		/*!
		 * Initiate the class.
		 * @i_numpart the number of particles
		 * @i_norb the number of orbitals in the single particle basis
		 * @vars pointer to instance of the imputVars class
		 * @pi_k pointer to array with x-components of orbitals (i_norb/2 elems)  
		 * @pi_l pointer to array with y-components of orbitals (i_norb/2 elems)  
		 * @pi_m pointer to array with z-components of orbitals (i_norb/2 elems)  
		 * @d_r0 the Wigner-Seitz radius
		 */
		virtual void init(
				int i_numpart,
				int i_norb,
				imputVars* vars,
				int* pi_k,
				int* pi_l,
				int* pi_m
				) =0;
		/*!
		 * Get the diagonal matrix element <D_i|H|D_i>.
		 * @pi_s array with the indices of the occupied arrays in |D_i>
		 * @return the diagonal matrix element <D_i|H|D_i>.
		 */
		virtual double getDiagonalElement(
			const int *pi_s) const =0;
		/*!
		 * Get the off diagonal matrix element <D_i|H{a_k^+ a_p}|D_i>.
		 * @i_sp the annihilated orbital
		 * @i_sk the created orbital
		 * @pi_s array with the indices of the occupied arrays in |D_i>
		 * @return the matrix element <D_i|H{a_k^+ a_p}|D_i>.
		 */
		virtual double getSingleExcitationAmplitude(
			const int i_sp, const int i_sk, const int *pi_s) const =0;
		/*!
		 * Get the off diagonal matrix element <D_i|H{a_p^+ a_q^+ a_l a_k}|D_i>.
		 * @i_sk the first annihilated orbital
		 * @i_sl the second annihilated orbital
		 * @i_sp the first created orbital
		 * @i_sq the second created orbital
		 * @return the off diagonal matrix element <D_i|H{a_p^+ a_q^+ a_l a_k}|D_i>.
		 */
		virtual double getDoubleExcitationAmplitude(
			const int iak, const int ial, const int iap, const int iaq) const =0;
};

#endif
