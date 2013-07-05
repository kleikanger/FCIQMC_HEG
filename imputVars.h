#ifndef IMPUTVARS_H
#define IMPUTVARS_H
/*
 * Container for all runtime parameters with get methods
 */

#include <iostream>
using std::string;

class imputVars
{/*startvimfold*/
	private:
		int i_numpart;
		int i_limit_nw;
		int i_num_loops; 
		int i_n_start_collecting_e;
		int i_update_s_frequency;
		int i_initiatorlimit;
		int i_numthreads;
		int i_active_orbs;
		int i_maxndets;
		int i_ranseed;
		double d_dt;
		double dS;
		double d_detweight;
		double d_pexone;
		double d_redistlimit;
		double d_xi;
		double d_rs;
		double d_kc;
		double d_mu;
		bool b_writeodata;  
		string s_ofpath;
		string s_system;
	public:
		void init();
		int iGETnumpart(){ return i_numpart; }
		int iGETlimit_nw(){ return i_limit_nw; }
		int iGETnum_loops(){ return i_num_loops; }
		int iGETstart_collecting_e(){ return i_n_start_collecting_e; }
		int iGETupdate_s_frequency(){ return i_update_s_frequency; }
		int iGETinitiatorlimit(){ return i_initiatorlimit; }
		int iGETnumthreads(){ return i_numthreads; }
		int iGETactive_orbs(){ return i_active_orbs;}
		int iGETmaxndets(){ return i_maxndets;}
		int iGETranseed(){ return i_ranseed; }
		double dGETdt(){ return d_dt; }
		double dGETdS(){ return dS; }
		double dGETdetweight(){ return d_detweight; }
		double dGETredistlimit(){ return d_redistlimit; }
		double dGETpexone(){ return d_pexone; }
		double dGETxi(){ return d_xi; }
		double dGETkc(){ return d_kc; }
		double dGETrs(){ return d_rs; }
		double dGETmu(){ return d_mu; } 
		bool bGETwriteodata(){ return b_writeodata; }
		string sGETofpath(){ return s_ofpath; }
		string sGETsystem(){ return s_system; }
		
		void iSETranseed(int i_ranseed){ this->i_ranseed = i_ranseed; }
		void iSETnumthreads(int i_numthreads){ this->i_numthreads = i_numthreads; };
};
#endif
