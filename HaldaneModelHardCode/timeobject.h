#ifndef TIMEOBJECT_H
#define TIMEOBJECT_H
#include<iostream>
#include<vector>



class timeobject: public timegrid
{
public:
	

	int something;
	double fR_double;
	double fI_double;
	complex *f;
	
	void initialize()
    {
		f =  (complex*)malloc(N*sizeof(complex));
        memset( f,    0,    sizeof(complex)*N );

		fR_double=0.0;
		fI_double=0.0;
		
    }
	
	
	void put_on_grid(timegrid g)
    {
		set_grid(g.N,g.dt,g.t0 );
		initialize();
    }
    
	/*** RungeKuttaIntegrator  ***/
	void integrateRK4()
	{
		//Double time grid and electric field
		complex alpha0 = complex(0.0,0.0);
		complex alpha1 = complex(0.0,0.0);
		complex alpha2 = complex(0.0,0.0);
		complex alpha3 = complex(0.0,0.0);
		
		vector<complex> v(N, 0.0);
		
		for (int ktime=0; ktime<N; ktime++)
            v[ktime] = f[ktime];
		//{
			//real(v[ktime])=real(f[ktime]);
			//imag(v[ktime])=imag(f[ktime]);
		//}
		
		//real(f[0])=0.0;
		//imag(f[0])=0.0;
		
		//ImplementatioN RK4th		
		for (int ktime=0; ktime<N-1 ; ktime++) 
		{
			
			alpha0 = v[ktime];
			 
			
			fR_double = real(v[ktime]) + ( real(v[ktime+1]) - real(v[ktime]) )/2. ;
			fI_double = imag(v[ktime]) + ( imag(v[ktime+1]) - imag(v[ktime]) )/2. ;
			
			
			alpha1 = complex( fR_double, fI_double );
			
			
			alpha2 = alpha1;
			
			
			alpha3 = v[ktime+1];	
			
            f[ktime+1] =  f[ktime] + dt/6.0*(alpha0 + 2.0*(alpha1 + alpha2) + alpha3);
			//real(f[ktime+1]) = ( real(f[ktime]) + dt/6.0*(real(alpha0 + 2.0*(alpha1 + alpha2) + alpha3) ) );

			
			//imag(f[ktime+1]) = ( imag(f[ktime]) + dt/6.0*(imag(alpha0 + 2.0*(alpha1 + alpha2) + alpha3) ) );
		
		}
			
		
	}//End rungekuta4th integrator
	
	/*~timeobject(){
		free(f);
	}*/
	
	
};
#endif
