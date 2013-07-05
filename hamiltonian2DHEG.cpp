#ifndef HAMILTONIAN2DHEG_CPP
#define HAMILTONIAN2DHEGT_CPP


#include <iostream>
#include <cmath>
#include <cstdlib>
#include "hamiltonian2DHEG.h"
#include "imputVars.h"

#define PI 3.141592653589793

using namespace std;

void hamiltonian2DHEG::init(
		int i_numpart,
		int i_norb,
	    imputVars* vars,
		int* pi_k,
		int* pi_l,
		int* pi_m
		)
{
	d_l = pow(PI*vars->iGETnumpart()*vars->dGETrs()*vars->dGETrs(), 1./2.); //2D
	d_k1 = 2.*PI*PI/(d_l*d_l);
	d_k2 = 2*PI/(d_l*d_l); // 2D
	d_mu2 = pow(vars->dGETmu(), 2.0);
	d_kk = pow(2.*PI/d_l, 2.0);

	this->i_numpart = i_numpart;
	this->pi_k = pi_k;
	this->pi_l = pi_l;
	
	// initiate single particle energies
	pd_spe = new double[i_norb];
	for (int i=0;i<i_norb;i++)
	{
		pd_spe[i] = d_k1 *
			( 	
				pi_k[i/2]*pi_k[i/2] + 
				pi_l[i/2]*pi_l[i/2] 
			);
	}
}

double hamiltonian2DHEG::getDiagonalElement(
		const int *pi_s) const
{

	int i,j,a,b;
	double d_eo = 0.0, d_et = 0.0;

	//calculate the single particle energies
	for (i=0;i<i_numpart;i++)
	{
		d_eo += pd_spe[pi_s[i]];
	}

	//XXX Correct to include here??
	//calculate the two body part of the expression 
	for (i=0;i<i_numpart-1;i++)
		for (j=i+1;j<i_numpart;j++)
		{
			// if (i%2==i%2) // && (j%2==j%2) osv
			// the direct term is zero because the total momentum is not conserved.

			// the exchange term
			if (pi_s[i]%2==pi_s[j]%2) // && (ial%2==iaq%2) 
				//second test alw true if first is because of sampling method
			{
				a=pi_s[i]/2;
				b=pi_s[j]/2;

				d_et -= 1./
					sqrt( d_mu2 +  // 2D
							d_kk*(pow(pi_k[a]-pi_k[b], 2) + 
							pow(pi_l[a]-pi_l[b], 2) 
							) );
			}
		}

	return d_eo + d_et*d_k2;
}

double hamiltonian2DHEG::getSingleExcitationAmplitude (
		const int i_sk,
		const int i_sp,
		const int *pi_s) const
{
	cout << "\n error : single excitation amplitudes always zero for the HEG \n";
	//always 0 for the HEG
	return 0;
}

double hamiltonian2DHEG::getDoubleExcitationAmplitude(
			const int iak, 
			const int ial,
			const int iap,
			const int iaq) const
{
	// note that the total momenta always are conserved becauce of the sampling function
	const int a = iak/2;
	double ret = 0.0;
	if (iak%2==iap%2) // && (ial%2==iaq%2) 
	{
		const int b = iap/2;
		//second test alw true if first is because of sampling method
		ret += 1./	
			sqrt( d_mu2 + // 2D
					d_kk *(
						pow(pi_k[a]-pi_k[b], 2) + 
						pow(pi_l[a]-pi_l[b], 2) 
						) ); 
	}
	// else ok here
	if (iak%2==iaq%2) // && (ial%2==iap%2) 
	{
		const int c = iaq/2;
		//second test alw true if first is because of sampling method
		ret -= 1./
			sqrt( d_mu2 + // 2D
					d_kk*(
						pow(pi_k[a]-pi_k[c], 2) + 
						pow(pi_l[a]-pi_l[c], 2) 
						) );
	}
	/*
	if ((pi_k[iak/2]+pi_k[ial/2])!=(pi_k[iap/2]+pi_k[iaq/2]))
		cout << "\nerr\n";
	if ((pi_l[iak/2]+pi_l[ial/2])!=(pi_l[iap/2]+pi_l[iaq/2]))
		cout << "\nerr\n";
	*/
	return ret*d_k2;
}
#endif
