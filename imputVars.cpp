#include "imputVars.h"
void imputVars::init()
{
 	d_dt=0.001557;
 	dS=5.000000;
 	d_rs=1.000000;
 	d_mu=0.000000;
 	i_numpart=10;
 	s_system="2DHEG";
 	d_kc=2.001000;
 	b_writeodata=1;
 	i_numthreads=2;
 	i_initiatorlimit=3;
 	i_limit_nw=2000;
 	i_num_loops=100000;
 	i_n_start_collecting_e=5000;
 	i_update_s_frequency=1;
 	s_ofpath="outputdata/10pt_2D_rs1_5shell_i4_nw2k";
 	d_pexone=0.000000;
 	d_xi=0.300000;
 	d_detweight=1.300000;
 	i_maxndets=1000000;
 	i_ranseed=-1;
 	d_redistlimit=0.040000;
}