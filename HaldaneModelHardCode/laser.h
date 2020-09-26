/*
 *  laser.cpp laser pulses train differents characteristics by pulse
 *
 *  Created by Alexis Chacon
 *
 */

#ifndef LASER_H
#define LASER_H
#include <stdlib.h>
#include <string>
#include "timegrid.h"
#include "timeobject.h"
#include "constant.h"
#include <vector>
#include <algorithm>
using namespace std;

const double gaussian_width_multiply = 2.0;
const double gaussian_factor = gaussian_width_multiply*gaussian_width_multiply* 4*log(2);

enum Envelope
{
    ENVELOPE_GAUSSIAN,
    //ENVELOPE_COS2,
    ENVELOPE_COS4,
    ENVELOPE_REC
};

struct LaserParam
{
public:
    double I0;			                    //Intensity per pulse
	double e;			                    //Ellipticity per pulse
    double E0;                              //Electric Field Magnitude


    double E0x;			                    //Electric Field Component in the x-direction
	double E0y;			                    //Electric Field Component in the y-direction

    double w0;			                    //Frequency per pulse
	double period0;		                    //Period per pulse
	double cycles0;		                    //Cycles number per pulse


	double cep0;		                    //Carrier envelope phase per pulse
	double phi_x;		                    //Additional phase in the x-direction given by elliptical polarization 
    double phi_y;                           //Additional phase in the y-direction given by elliptical polarization 

    double twidth;		                    //Time bandwidth per pulse
    double t0;                              // Mid-time of pulses
    double theta0;
    
	Envelope envelope;	//Envelope type used

    // Essential initial values - I0, e, w0, cycles0, cep0, phi_rel, t0, theta0
    void Initialize(double _E0, double _e, double _w0, double _ncycle, double _cep,  double _t0, double _theta0, string _env_name)
    {
        E0 = _E0;   e = _e; w0 = _w0;   cycles0 = _ncycle;  cep0 = _cep;  t0 = _t0;   theta0 = _theta0;
        /*if (_env_name == "COS2" || _env_name == "sin2")
        {
            envelope = ENVELOPE_COS2;
        }*/
        if (_env_name == "cos4")
        {
            envelope = ENVELOPE_COS4;
        }
        else if (_env_name == "gauss")
        {
            envelope = ENVELOPE_GAUSSIAN;
        }
        else
        {
            cout << "Unknown envelope\n";
            exit(1);
        }
    }
    void PreProcessing()
    {
        
        I0 = E0*E0 * 3.5e16;

        period0    =  dospi/w0; // laser period or cycle
        
        // laser time duration
        switch (envelope)
        {
        case ENVELOPE_GAUSSIAN:
            twidth = gaussian_width_multiply * cycles0*period0;
            break;
        
        default:
            twidth = cycles0*period0/2;
            break;
        }
        
        // The ellipticity is defined by e = b/a. and -1 <= e <= 1
        // theta0 is degree between major axis and x-axis
        // If there is no additional phase, E(t=t0) = {E0 cos(theta0), E0 sin(theta0)}
        // i.e. additional phase is adjusted to electric field is at major axis at t=t0
        // Also note that theta0 and pi-theta0 represents same major axis, but phase of electric field is changed
        // i.e. E(theta0 = pi-at0) = -E(theta0 = at0)
        // If e = 0, the laser is circularly polarized independent of value of theta0
        // however, theta0 give additional CEP. additional CEP = theta0 in this case
        // e < 0 --> right helicity, e > 0 --> lefet helicity, e = 0 --> linear polarization

        // Note that here we assume that laser filed is propagating in +z direction
        // i.e. E is proportional to cos(kz - wt + phase), so cos(wt - phase) is used in code.

        // For e= 0 and phi_rel = 0, the laser is linearly polarized along x-direction
        // For e= 0 and phi_rel = 0, and theta != 0 (different to zero), the laser is linearly
        //polarized along a stringht line with theta angle with respect to positive x-direction
        // For e= 1 and phi_rel = pi/2., the laser is circularly polarized
        // Between e = (0, 1) and phi_rel = pi/2., the laser is ellipticaly polarized with major axis along x
        
        // not fully tested yet!

        double temp_a = 1.0 / sqrt(1.0 + e*e) * E0;
        double temp_b = temp_a * e;

        E0x = sqrt(pow(temp_a*cos(theta0), 2) + pow(temp_b*sin(theta0), 2));
        E0y = sqrt(pow(temp_a*sin(theta0), 2) + pow(temp_b*cos(theta0), 2));

        phi_x = atan2(-temp_b*sin(theta0), temp_a*cos(theta0));
        phi_y = atan2( temp_b*cos(theta0), temp_a*sin(theta0));
    };
};

class laser 
{

public:

	
    
    int Npulses;				//Total pulses number
    int Nt;

    vector<LaserParam> PulseParam; // paramters for pulses

	double dt;					//Time step
	double blaser;				//Time before of the pulses
	double alaser;				//Time after of the pulses

    double ncycles2;
    
	double major0;				//Major value
	double minus0;				//Minus value
    double width;
    complex a_ef0,a_ef1;
    complex a_af0,a_af1;
    double atmin, atmax;
    
    
	vector<double> clock0;      //Initial time per pulse
	vector<double> clock1;      //The time until half per pulse
	vector<double> clock2;		//The time until final per pulse
	
	
    
    
	/*==========================*/
	/*      MAIN FUNCTIONS
	/*==========================*/
    laser(){};
    laser(int _Npulses);
	void init_laser(int _Npulses);      							                                // Creator Object
	//~laser();														//Destructor
    void laser_pulses(double _dt, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
    
    
	/*=========================*/
	/*   SECUNDARY FUNCTIONS
	/*=========================*/
	inline void Initialize_Parameters();			// Initialize paramters like ellipticity and duration
	inline void Set_StartTime_EndTime();		        // "Start" and "end" time per pulse
	
    void Set_Time_Domain();								// set time domain
	

    double f_envelope(double _t, LaserParam const *_param1, bool _isDt);    // Generate enevelop function
    
    
    complex avlaser( double const *t );
    complex elaser( double const *t );
    inline complex avlaser( double const *t , LaserParam const *_param1);
    inline complex elaser( double const *t , LaserParam const *_param1);

    void Print_LaserInfo();
    
    
};


//==============================================================================//
			/*=== MAIN FUNCTIONS ===*/

/*===========================================================          		
		/*=== OBJECT'S LASER CONSTRUCTOR  ===*/
laser::laser(int _Npulses)
{
    init_laser(_Npulses);
}
void laser::init_laser(int _Npulses)
{
    Npulses = _Npulses;

    clock0.resize( Npulses, 0.0 );
    clock1.resize( Npulses, 0.0 );
    clock2.resize( Npulses, 0.0 );
    
    PulseParam.resize(Npulses);
    
    a_ef0 = 0.;
    a_af0 = 0.;
    a_ef1 = 0.;
    a_af1 = 0.;
}//End initialize variable  



/*===========================================================          		
 /*===  LASER PULSES FUNCTION  ===*/
void laser::laser_pulses(double _dt, double _blaser, double _alaser)
{
    
    
   	/*======= Pulse's Parameter ======*/
   	dt      = abs( _dt );			                    //	Time step
    
    
    
    blaser  = abs(_blaser);	                        //	Time before the laser or the pulse train
   	alaser  = abs(_alaser);	                        //	Time after the laser or the pulse train
	
    
    
    
   	//================================//
	Initialize_Parameters();					//	Initialize max amplitude and period        
   	Set_StartTime_EndTime();	 	     	        //	"Start" and "end" time by each pulse
	
    
    
    
   	//major0	= qmajor(clock2);	     	                //	Maximum time to all pulse
   	//minus0	= qminus(clock0);          	            //	Minimun time to all pulse     
    major0 = *max_element(clock2.begin(), clock2.end());
    minus0 = *min_element(clock0.begin(), clock0.end());
    
    Set_Time_Domain();       	     	                //	Objet time axis


}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/***********************************************************/
		/*========= SECUNDARY FUNCTIONS ===========*/
/***********************************************************/

//== Function initialize max amplitude and period ==// 
inline void laser::Initialize_Parameters()
{
	for( int kpulse = 0; kpulse < Npulses; kpulse++ )
	{
        PulseParam[kpulse].PreProcessing();
	}

}


/***********************************************************/
  //== Function "start" and "end" time by each pulse ==//
/***********************************************************/
inline void laser::Set_StartTime_EndTime()
{
    
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
        clock0[kpulse] = PulseParam[kpulse].t0 - PulseParam[kpulse].twidth;
        clock1[kpulse] = PulseParam[kpulse].t0;
        clock2[kpulse] = PulseParam[kpulse].t0 + PulseParam[kpulse].twidth;
        
	}//End loop
 
    
}//End function



/***********************************************************/
                //== Time axis and field set  ==//
/***********************************************************/
void laser::Set_Time_Domain()   
{
    atmin = minus0 - blaser;
    atmax = major0 + alaser;

    a_ef1 = elaser(  &(atmin) );
    a_af1 = avlaser( &(atmin) );
    
    Nt   = floor( (atmax - atmin)/dt );

    if ( Nt%2 !=0 ) Nt+=1; 
    
}
//== End of Time axis and field set  ==//
/***********************************************************/

// A(t) = A0 f(t) cos(w0t - cep) = E0/w0 f(t) cos(w0t - cep)
// E(t) = -dA(t)/dt = - E0/w0 f'(t) cos(w0t - cep) - E0 f(t) sin(w0t - cep)

// _t : actual time value, isDt : true for df(t)/dt, false for f(t)
double laser::f_envelope(double _t, LaserParam const * p1, bool _isDt)
{
    switch (p1->envelope)
    {
    /*case ENVELOPE_COS2:
        if (abs(_t) > p1->twidth) return 0;
        if (_isDt) return -p1->w0/(2.*p1->cycles0)*sin(_t* p1->w0 /(p1->cycles0));
        else return pow( cos(_t* p1->w0/(2*p1->cycles0) ), 2. );
        break;*/
    
    case ENVELOPE_COS4:
        if (abs(_t) > p1->twidth) return 0;
        if (_isDt) return -4.*p1->w0/(2.*p1->cycles0) *sin(_t* p1->w0 /(2*p1->cycles0)) * pow( cos(_t* p1->w0/(2*p1->cycles0) ), 3. );
        else return pow( cos(_t* p1->w0/(2*p1->cycles0) ), 4. );
        break;
    
    default:
        if (_isDt) return -2*gaussian_factor* _t*exp( -gaussian_factor* _t*_t/p1->twidth/p1->twidth )/(p1->twidth*p1->twidth);
        else return exp( -gaussian_factor* _t*_t/p1->twidth/p1->twidth );
        break;
    }
}

complex laser::elaser( double const *t )
{
    //if (*t < minus0 || *t > major0) return 0;
    a_ef0=0.;
    LaserParam *p1;
    for (int kpulse = 0; kpulse < Npulses; ++kpulse)
    {
        p1= &(PulseParam[kpulse]);
        a_ef0 += elaser(t, p1);
    }
    return a_ef0-  a_ef1;
}
inline complex laser::elaser( double const *t, LaserParam const *p1)
{
    double kt = *t - p1->t0;
    return complex(p1->E0x*cos(kt* p1->w0 - p1->cep0 - p1->phi_x) * f_envelope(kt, p1, false)
            + p1->E0x/p1->w0 * sin(kt * p1->w0 - p1->cep0 - p1->phi_x) * f_envelope(kt, p1, true),
            p1->E0y * cos(kt * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope(kt, p1, false)
            + p1->E0y/p1->w0 * sin(kt * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope(kt, p1, true));
}
// Value of phase is adjust to E0x = E0(t=t0, cep = 0)
// So E(cos, sin) --> A(sin, -cos)
complex laser::avlaser( double const *t )
{
    //if (*t < minus0 || *t > major0) return 0;
    a_af0=0.;
    LaserParam *p1;
    for (int kpulse = 0; kpulse < Npulses; ++kpulse)
    {
        p1= &(PulseParam[kpulse]);
        a_af0 += avlaser(t, p1);
    }
    return a_af0-  a_af1;
}
inline complex laser::avlaser( double const *t, LaserParam const *p1)
{
    double kt = *t - p1->t0;
    return complex(-p1->E0x/p1->w0 * sin(kt * p1->w0 - p1->cep0 - p1->phi_x) * f_envelope(kt, p1, false),
            -p1->E0y/p1->w0 * sin(kt * p1->w0 - p1->cep0 - p1->phi_y) * f_envelope(kt, p1, false));
    /*return complex(p1->E0x/sqrt(1. + p1->e*p1->e) * (cos(p1->w0 * kt) * (-1. + p1->cycles0*p1->cycles0 + p1->cycles0*p1->cycles0 * cos((kt * p1->w0)/p1->cycles0) 
        + p1->cycles0 * sin(kt * p1->w0) * sin((kt * p1->w0) / p1->cycles0) )) / (2. * (-1. + p1->cycles0) * (1. + p1->cycles0) * p1->w0), 0. );*/
}

void laser::Print_LaserInfo()
{
    cout << "\n\n*************************";
    cout << "\nNtime = " << Nt;
    cout << "\n*************************** \n\n";
}

#endif
//END
