// *** grid.h ***  for the hop_cart project -----------------  db 09/98
#ifndef TIMEGRID_H
#define TIMEGRID_H
#include<iostream>
#include<math.h>
#include<vector> 
#include "constant.h"
using namespace std;


class timegrid
	{
	public:
		
		//string name;
		int N;
		
		double dt;
		double t0;
		
		
		vector<double> t;
		
		
		/********************************************/
		// Useful data
		
		inline double sizet_au() const { return dt*N; }		
		inline double vol_elem() const { return dt; }
		
		inline long tot_pts() const { return N; }
		inline double aprox_size_Mb() const{ return N*16.*1e-6; }
		
		
		/********************************************/
		// Initialize
		inline void set_dt(double _t0, double _dt)
		  {
			dt=_dt;
			t0=_t0;
		  }
		
		
		
		inline void set_n(int _n)
		  {
			N=_n;
		  }
		
		
		
		inline void set_grid(int _n, double _dt, double _t0)
		  {			  
			  set_n(  _n);
			  set_dt(_t0, _dt);
			  
			  t.resize(N, 0.);			
			  
			  for(int i=0;i<N;i++)
			  	  t[i] = t0 + i*dt;	//t[i]	=-n*dt/2+i*dt;
			  
		  }	
}; //End of the grid.h 
#endif
