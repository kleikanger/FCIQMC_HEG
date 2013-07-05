#include "imputVars.h"
void imputVars::init()
{
	
	d_dt = 0.00500;
 	dS = 5.00;
	d_rs =  0.5;
	d_mu = 0.0;

#if 0
	i_numpart = 10;

	s_system = "2DHEG";
	//d_kc=1; //N 10 //2 shells
	d_kc=1.5; //N 18 //3 shells
	//d_kc=2; //N 26 //4 shells
	//d_kc=2.5; //N 42 //5 shells
#else
	s_system = "3DHEG";
 	i_numpart = 14; // 2,14,38,54
	
	// change d_kc, INUM_ORBITALS sim.	
	//d_kc=1;
	//N=14
	d_kc = 1.45; //N=38
	//d_kc = 2.; //N = 66
	//d_kc =  2.4; //N = 114
	//d_kc =  2.9; //N = 186
#endif

	b_writeodata=1;
 	i_numthreads=2;
 	i_initiatorlimit=3;

 	i_limit_nw=10000;
 	i_num_loops=100000;
 	i_n_start_collecting_e=5000;
 	i_update_s_frequency=1;
	s_ofpath="outputdata/3D_pt10_sh3_rs05_i3_NW10k_dt0005"; 
	
	d_pexone=0.0;//probability of a double excitation
	d_xi=0.3;//see:umrigars update formula. S(n+1) = S(n)-(d_xi/d_dt)*(N_W(n)/N_W(n+1)) 
 	
	//i_active_orbs=12;
 	
 	d_detweight=1.00000;
 	i_maxndets=1000000;
 	i_ranseed=1;
 	d_redistlimit=0.040000;
}
