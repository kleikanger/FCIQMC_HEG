#ifndef basis3DHEG_CPP
#define basis3DHEG_CPP

#include <cstdlib>
#include <bitset>
#include <iostream>
#include <algorithm>
#include "basis3DHEG.h"

#define PI 3.1415923585

using std::bitset;
using std::cout;

template <size_t N>
basis3DHEG<N>::basis3DHEG(double d_kc)
{
	this->d_kc = d_kc;
	d_kc2 = d_kc*d_kc;
	i_kc = static_cast<int>(d_kc);

	i_ms = 4*i_kc+1;
	i_ms2 = i_ms*i_ms;
	i_ns = 2*i_kc+1;
	i_ns2 = i_ns*i_ns;

	pi_k = new int[N];
	pi_l = new int[N];
	pi_m = new int[N];
	pb_asets = new bitset<N>[i_ms*i_ms*i_ms];
	pb_psets = new bitset<N>[i_ms*i_ms*i_ms];
	pi_aprobs = new int[i_ms*i_ms*i_ms];
	pi_pprobs = new int[i_ms*i_ms*i_ms];
	pi_anumstates = new int[i_ms*i_ms*i_ms];
	pi_pnumstates = new int[i_ms*i_ms*i_ms];
	pi_m2i = new int[i_ns*i_ns*i_ns];

	initBasis();
	initMBS();
	initMom2Indx();
	initNumOrbsPerM();
}

template <size_t N>
basis3DHEG<N>::~basis3DHEG()
{
	delete [] pi_k;
	delete [] pi_l;
	delete [] pi_m;
	delete [] pb_asets;
	delete [] pb_psets;
	delete [] pi_aprobs;
	delete [] pi_pprobs;
	delete [] pi_anumstates;
	delete [] pi_pnumstates;
	delete [] pi_m2i;
}

template <size_t N>
void basis3DHEG<N>::initBasis()
{
	int n = 0;

	//initiate the basis
	for (int i=-i_kc;i<=i_kc;i++)
		for (int j=-i_kc;j<=i_kc;j++)
			for (int k=-i_kc;k<=i_kc;k++)
				if ( ( i*i+j*j+k*k )<=d_kc2 ) //XXX <=
				{
					pi_k[n]=i;
					pi_l[n]=j;
					pi_m[n]=k;

					n+=1;
				}


	int *pi_mom = new int[N];
	int *pi_indx = new int[N];
	int *pi_tk = new int[N];
	int *pi_tl = new int[N];
	int *pi_tm = new int[N];

	for (int i=0;i<N;i++)
		pi_mom[i] = pi_k[i]*pi_k[i]+
					pi_l[i]*pi_l[i]+
				 	pi_m[i]*pi_m[i];
	for (int i=0;i<N;i++)
		pi_indx[i] = i;

	// sort indx array
	std::sort(&pi_indx[0], &pi_indx[N], LessThan(pi_mom) );
	// reorder
	for (int i=0;i<N;i++)
	{
		pi_tk[i] = pi_k[pi_indx[i]];
		pi_tl[i] = pi_l[pi_indx[i]];
		pi_tm[i] = pi_m[pi_indx[i]];
	}
	
	// swap pointers
	int* pi_swap = pi_tk;
	pi_tk = pi_k;
	pi_k = pi_swap;
	
	pi_swap = pi_tl;
	pi_tl = pi_l;
	pi_l = pi_swap;
	
	pi_swap = pi_tm;
	pi_tm = pi_m;
	pi_m = pi_swap;
			
	delete [] pi_mom;
	delete [] pi_indx;
	delete [] pi_tk;
	delete [] pi_tl;
	delete [] pi_tm;
}

template<size_t N>
void basis3DHEG<N>::initMBS()
{
	int i_indx, a, b, c;

	//reset all (should be done by constructor)
	for (int i=0;i<i_ms*i_ms*i_ms;i++)
		pb_psets[i].reset();
	for (int i=0;i<i_ms*i_ms*i_ms;i++)
		pb_asets[i].reset();
	for (int i=0;i<i_ms*i_ms*i_ms;i++)
		pi_aprobs[i] = 0;
	for (int i=0;i<i_ms*i_ms*i_ms;i++)
		pi_pprobs[i] = 0;


	for (int i=0;i<N;i++)
		for (int j=i;j<N;j++)
		{
			a = pi_k[i]+pi_k[j]+2*i_kc;
			b = pi_l[i]+pi_l[j]+2*i_kc; 
			c = pi_m[i]+pi_m[j]+2*i_kc; 
			i_indx = a+b*i_ms+c*i_ms2;

			if (i==j)
			{
				pb_asets[i_indx].set(i);

				//TODO correct ?
				pi_aprobs[i_indx] += 1;
			}
			else
			{
				pb_asets[i_indx].set(i);
				pb_asets[i_indx].set(j);

				pb_psets[i_indx].set(i);
				pb_psets[i_indx].set(j);

				//TODO correct ?
				pi_aprobs[i_indx] += 2;
				pi_pprobs[i_indx] += 1;
			}
		}
	/*
	   p_gen^-1: 
	   paralell spins: pprobs: multiply with 2
	   anti-parallel spins: aprobs: multiply with 4 if i!=j and 2 if i==j
	   (because of spin symmetries)
	   */
}

template <size_t N>
void basis3DHEG<N>::initNumOrbsPerM()
{
	int i_indx, a, b, c;

	for (int i=0;i<N;i++)
		for (int j=i;j<N;j++)
		{
			a=pi_k[i]+pi_k[j]+2*i_kc;
			b=pi_l[i]+pi_l[j]+2*i_kc; 
			c=pi_m[i]+pi_m[j]+2*i_kc; 
			i_indx = a+b*i_ms+c*i_ms2;

			pi_anumstates[i_indx] = pb_asets[i_indx].count();
			pi_pnumstates[i_indx] = pb_psets[i_indx].count();
		}
}

template <size_t N>
void basis3DHEG<N>::initMom2Indx()
{
	int a, b, c;

	for (int i=0;i<i_ns*i_ns*i_ns;i++)
		pi_m2i[i] = -1;

	for (int i=0;i<N;i++)
	{
		a = pi_k[i]+i_kc;
		b = i_ns*(pi_l[i]+i_kc);
		c = i_ns2*(pi_m[i]+i_kc);
		pi_m2i[a+b+c] = i;
	}
}

template <size_t N>
int basis3DHEG<N>::getOrb(const int a, const int b, const int c)
{	
	const int p = a+i_kc;
	const int q = (b+i_kc)*i_ns;
	const int r = (c+i_kc)*i_ns2;

	return pi_m2i[p+q+r];
}

template <size_t N>
void basis3DHEG<N>::sampleInteractionParallelSpins(
		int &i_c1, int &i_c2, int &i_invprob, 
		int i_a1, int i_a2, int i_spin, double d_ran)
{
	// i_a2 >= i_a1
	i_a1/=2;
	i_a2/=2;

	const int a = pi_k[i_a1]+pi_k[i_a2];
	const int b = pi_l[i_a1]+pi_l[i_a2];
	const int c = pi_m[i_a1]+pi_m[i_a2];

	int p = a+2*i_kc;
	int q = (b+2*i_kc)*i_ms;
	int r = (c+2*i_kc)*i_ms2;

	i_invprob = pi_pprobs[p+q+r];

	//TODO Precount
	//const int i_dimstatespace = pb_psets[p+q+r].count();
	const int i_dimstatespace = pi_pnumstates[p+q+r];

	i_c2 = static_cast<int>(d_ran*i_dimstatespace); //betw 0 and dim-1
	i_c1 = pb_psets[p+q+r]._Find_first();
	for (int i=0;i<i_c2;i++)
		i_c1 = pb_psets[p+q+r]._Find_next(i_c1);

	p = a+i_kc-pi_k[i_c1];
	q = (b+i_kc-pi_l[i_c1])*i_ns;
	r = (c+i_kc-pi_m[i_c1])*i_ns2;

	i_c2 = pi_m2i[p+q+r];

	i_c1 = 2*i_c1+i_spin;
	i_c2 = 2*i_c2+i_spin;
#if 0		
	cout 
		<< "\ns"
		<< "\n" << pi_k[i_a1]+pi_k[i_a2] << ", " << pi_k[i_c1/2]+pi_k[i_c2/2] 
		<< "\n" << pi_l[i_a1]+pi_l[i_a2] << ", " << pi_l[i_c1/2]+pi_l[i_c2/2] 
		<< "\n" << pi_m[i_a1]+pi_m[i_a2] << ", " << pi_m[i_c1/2]+pi_m[i_c2/2] 
		<< "\n";
#endif
}

template <size_t N>
void basis3DHEG<N>::sampleInteractionAntiparallelSpins(
		int &i_c1, int &i_c2, int &i_invprob, 
		int i_a1, int i_a2, 
		double d_ran1, double d_ran2)
{
	// i_a2 >= i_a1
	i_a1/=2;
	i_a2/=2;

	const int a = pi_k[i_a1]+pi_k[i_a2];
	const int b = pi_l[i_a1]+pi_l[i_a2];
	const int c = pi_m[i_a1]+pi_m[i_a2];

	int p = a+2*i_kc;
	int q = (b+2*i_kc)*i_ms;
	int r = (c+2*i_kc)*i_ms2;

	i_invprob = pi_aprobs[p+q+r];

	//TODO Precount
	//const int i_dimstatespace = pb_asets[p+q+r].count();
	const int i_dimstatespace = pi_anumstates[p+q+r];
	//XXX what if statestpace has size 0. then break and return -1

	i_c2 = static_cast<int>(d_ran1*i_dimstatespace); //betw 0 and dim-1
	i_c1 = pb_asets[p+q+r]._Find_first();
	for (int i=0;i<i_c2;i++)
		i_c1 = pb_asets[p+q+r]._Find_next(i_c1);

	p = a+i_kc-pi_k[i_c1];
	q = (b+i_kc-pi_l[i_c1])*i_ns;
	r = (c+i_kc-pi_m[i_c1])*i_ns2;

	i_c2 = pi_m2i[p+q+r];

	if (d_ran2<.5)
	{
		i_c1 = 2*i_c1+1;
		i_c2 = 2*i_c2;
	}
	else
	{
		i_c1 = 2*i_c1;
		i_c2 = 2*i_c2+1;
	}
#if 0
	cout
		<< "\nas"	
		<< "\n" << pi_k[i_a1]+pi_k[i_a2] << ", " << pi_k[i_c1/2]+pi_k[i_c2/2] 
		<< "\n" << pi_l[i_a1]+pi_l[i_a2] << ", " << pi_l[i_c1/2]+pi_l[i_c2/2] 
		<< "\n" << pi_m[i_a1]+pi_m[i_a2] << ", " << pi_m[i_c1/2]+pi_m[i_c2/2] 
		<< "\n";
#endif
}
		
template <size_t N>
void basis3DHEG<N>::getPointersToMomVecs(int* &pi_k_, int* &pi_l_, int* &pi_m_)
{
	pi_k_ = pi_k;
	pi_l_ = pi_l;
	pi_m_ = pi_m;
}

#endif
