//
//  solidstructure.h
//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
// Structural information defined as follow:

/***********************************************
 
 (1) Energy dispersions,
 (2) Dipole matrix elements,
 (3) Berry connection and
 (4) Berry curvature,
 (5) Group velocities,
 (6) Anomalous velocities
 (7) Band gap
 (8) Chern No.
 
**********************************/

#ifndef solidstructure_h
#define solidstructure_h

#include <stdlib.h>
#include <string>
//#include "mkl.h"
#include "constant.h"
#include "momaxis.h"

using namespace std;

class solidstructure
{
  
private:
    
    
    //Indexes on loops for
    int i, j, l, n;
    int ai, aj, al;
    
    
    
public:
    
    /**********************************************
     
     INPUTS PARAMETERS, VARIABLES AND OBJECTS
     
     */
    
    int rflag;
    
    //Hopping parameters
    double t1, t2, phi0, M0;
    
    //Lattice constant
    //double lattice_a0;
    
    momaxis *g; //Momentum grid input
    

    
    
    /***********************
     VARIABLES, for OUTPUTS
     ******/
    double *energy_gap;            // Energy gap difference
    double ks00[Ngrad];
    double ks01[Ngrad];
    
    double ks10[Ngrad];
    double ks11[Ngrad];
    double ks12[Ngrad];
    double ks13[Ngrad];
    

    double ks20[Ngrad];
    double ks21[Ngrad];
    double ks22[Ngrad];
    double ks23[Ngrad];
    
    
    //Variables
    double *chi_gap[Ngrad];        // Berry Connection difference
    double *groupvel_c[Ngrad];
    double *groupvel_v[Ngrad];
    
    
    double *zbcurva_v;
    double *zbcurva_c;
    double T1,T2;

    
    
    /************************************
     VARIABLES, var., for PROCESSES
            ***************************/
    //Basic vectors for Honneycomb Lattice
    double vec_a[Ngrad][Nvects];
    double vec_b[Ngrad][Nvects];
    
    
    
    
    //Variables for B0, B1, B2, B3 evaluations and their gradient along x and y directions
    double b0res0, b1res0, b2res0, b3res0, bnorm ;
    double xb0der0, xb1der0, xb2der0, xb3der0, xbnorm0;
    double yb0der0, yb1der0, yb2der0, yb3der0, ybnorm0;
    
    int gauge;
    
    double wave_norm, wave_norm2;
    double gauge_b_set[Nvects];
    double gauge_b_set_grad[Ngrad][Nvects];
    double w_norm_grad[Ngrad];
    double grad_bnorm[Ngrad];
    double xtemp[Ngrad*2];
    double ytemp[Ngrad*2];
    
    
    //Variables for energy dispersion
    double energyv0, energyc0, energyg0, eg0;
    
    
    double zBcurva_c0, zBcurva_v0;
    double xcberryconn0, xvberryconn0, xchig0;
    double ycberryconn0, yvberryconn0, ychig0;
    double chig0[Ngrad]; 
    double MinEGap, MaxEGap; 
    
    double re_width;
    double reg0;
    
    
    //Anomalous velocity variables
    double xanomalous_v0, yanomalous_v0;
    double xanomalous_c0, yanomalous_c0;
    
    
    
    //Group velocity vars.
    double group_v0[Ngrad];
    double group_c0[Ngrad];
    

    
    
    //Dipole matrix element from val. -> cond. band vars.
    complex curvatemp1, curvatemp2;
    
    complex wv_grad0[Ngrad][Ngrad];
    complex wc_grad0[Ngrad][Ngrad];
    complex velocity_cv0[Ngrad];
    
    complex dip_cv0[Ngrad];
    complex tempc[Ngrad];
    
    complex *dipole_cv[Ngrad];
    
    double chi_c0[Ngrad];
    double chi_v0[Ngrad];
    
    
    //Smoothed or regularization parameter
    double eps_phi0;
    
    
    //K and K' points
    double K[Ngrad];
    double Kprime[Ngrad];
    
    
    //Maxima momentum
    double honeycomb_lattice_kmax[Ngrad];
    
    double flag;
    double angle_a0, angle_b0;
    
    
    complex uwave_val0[Ngrad];
    complex uwave_cond0[Ngrad];
    
    
    
    /********************************
         OWN METHODS OR FUNCTIONS FOR
         INPUTS, PROCESSES and OUTPUTS
        *******************************/

    solidstructure(  const momaxis *_g ,const double *_a0, const double *_dephasing, const double *_dephasing2 );
    ~solidstructure( );
    void set_memory();
    
    
    
    
    void haldane_model_params( double const *_t1
                              ,double const *_t2
                              ,double const *_phi0
                              ,double const *_M0
                              ,int const _gauge );
    
    
    
    void bases_vectors( );
    
    void Set_Of_B_BGrad( double const *kx, double const *ky );
    
    
    void Ec( );
    void Ev( );
    void Eg( );
    

    void BerryConnectionCV( );
    void dipoleCV( );
    void zBerryCurvaCV( );
    void group_velCV( );

    double min( double *v, int Nsize );
    double max( double *v, int Nsize );
    

    double myRegularization( double const *kx, double const *ky, double *skx, double *sky  );
    
    void lattice_structure( );
    
    void anomalous_velCV( complex const efield, double const zBerryCurvature );
    
    void cdot( const complex *v1, const complex *v2, const int Nv, complex *res);
    
    double ChernNumber( );
    
    /*
     
     DATA FORM .TXT OR .DAT OUTPUTS
     
     */
    //Outputs
    void energy_vc_output(      FILE *output, int skip0, int skip1, int skip2 );
    void dipole_cv_output(      FILE *output, int skip0, int skip1, int skip2 );
    void connection_c_output(   FILE *output, int skip0, int skip1, int skip2 );
    void curvature_c_output(    FILE *output, int skip0, int skip1, int skip2 );
    void group_vel_vc_output(   FILE *output, int skip0, int skip1, int skip2 );
    
    
};


//######################################
solidstructure::solidstructure(  const momaxis *_g
                                ,const double *_a0
                                ,const double *_dephasing
                                ,const double *_dephasing2
                               )
{
    

    //##################################
    //lattice_a0  = *_a0;
    T1          = *_dephasing;
    T2          = *_dephasing2;
    flag        = 0;
    
    
    //Regularization parameters by defaul
    reg0        = 0.;
    re_width    = 0.11;  //5% of total ky-momentum-length
    rflag       = 0;
    
    
    //Momentum axis and grid
    g = new momaxis( _g->N, _g->kmaxs, &(_g->lattice_a0) );
    g->integral_method( _g->imethod );
    

    //Bravais vector-basis
    bases_vectors( );
    set_memory( );
    
}




//######################################
solidstructure::~solidstructure( )
{
   
    
    delete g;
    free(energy_gap);
    
    
    for (n = 0; n < Ngrad; n++ )
    {
        
        
        free( chi_gap[n] );
        free( groupvel_c[n] );
        free( dipole_cv[n] );
        
        
    }
    
    
    
    free(zbcurva_c);
    free(zbcurva_v);
    
    
}


double solidstructure::myRegularization( double const *kx, double const *ky, double *skx, double *sky  )
{
    
    
    return reg0*exp( - ( ( *kx - *skx)*( *kx - *skx) + ( *ky - *sky)*( *ky - *sky) )/re_width/re_width );
    
    
}


//######################################
void  solidstructure::haldane_model_params( double const *_t1
                                           ,double const *_t2
                                           ,double const *_phi0
                                           ,double const *_M0
                                           ,int const _gauge )
{
    
    t1          = *_t1;
    t2          = *_t2;
    phi0        = *_phi0;
    M0          = *_M0;
    eps_phi0    = t1/10.;
    gauge       = _gauge;

    

    
}

//######################################
void solidstructure::lattice_structure()
{
  
    flag = 1;
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++ )
            for ( i = 0; i < g->N[0]; i++ )
            {
                
                Set_Of_B_BGrad( &g->k[0][i]
                         ,&g->k[1][j]
                         );
                
                zBerryCurvaCV();
                group_velCV( );
                
                energy_gap[ g->index(&i,&j,&l) ]       = eg0;
                
                chi_gap[0][ g->index(&i,&j,&l) ]       = chig0[0];
                chi_gap[1][ g->index(&i,&j,&l) ]       = chig0[1];
                
                dipole_cv[0][ g->index(&i,&j,&l) ]     = dip_cv0[0];
                dipole_cv[1][ g->index(&i,&j,&l) ]     = dip_cv0[1];
                
                zbcurva_c[ g->index(&i,&j,&l) ]        = zBcurva_c0;
                zbcurva_v[ g->index(&i,&j,&l) ]        = zBcurva_v0;
                
                groupvel_c[0][ g->index(&i,&j,&l) ]    = group_c0[0];
                groupvel_c[1][ g->index(&i,&j,&l) ]    = group_c0[1];
                
                groupvel_v[0][ g->index(&i,&j,&l) ]    = group_v0[0];
                groupvel_v[1][ g->index(&i,&j,&l) ]    = group_v0[1];
                
        }
    

       MinEGap = min( energy_gap, g->N[2]*g->N[1]*g->N[0] ); 
       MaxEGap = max( energy_gap, g->N[2]*g->N[1]*g->N[0] ); 

}



//######################################
double solidstructure::min( double *v, int Nsize )
{

	double min0 = v[0];
	
	for ( n =1 ; n<Nsize; n++)
		if (v[n]<min0) min0=v[n];	
	
	return min0;

}//End minus



//######################################
double solidstructure::max( double *v, int Nsize)
{

	double max0 = v[0];
	
	for (n=1;n<Nsize;n++)   
		if (v[n]>max0) max0=v[n];
	
	return max0;

}//End major


//######################################
void solidstructure::set_memory()
{

    energy_gap         = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;   // Energy gap difference
    memset( energy_gap, 0, (g->Ntotal)*sizeof( double ) );
    
    
    zbcurva_v    = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;
    memset( zbcurva_v,   0, (g->Ntotal)*sizeof( double ) );
    
    
    zbcurva_c    = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;
    memset( zbcurva_c,   0, (g->Ntotal)*sizeof( double ) );
    
    
    for (n = 0; n < Ngrad; n++ )
    {
    
        
        chi_gap[n]     = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;
        groupvel_c[n]  = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;
        groupvel_v[n]  = (double*) malloc( (g->Ntotal) * sizeof( double ) ) ;
        
        dipole_cv[n]   = (complex*) malloc( (g->Ntotal) * sizeof( complex ) ) ;
        
        memset( chi_gap[n],    0, (g->Ntotal)*sizeof( double ) );
        memset( groupvel_c[n], 0, (g->Ntotal)*sizeof( double ) );
        memset( groupvel_v[n], 0, (g->Ntotal)*sizeof( double ) );

        memset( dipole_cv[n],  0, (g->Ntotal)*sizeof( complex ) );
        

    }
    

}


//######################################
void solidstructure::bases_vectors()
{
    //First vector base: position atom
    vec_a[0][0] = 0.;
    vec_a[1][0] = g->lattice_a0;
    
    vec_a[0][1] = -sqrt(3.)/2.*(g->lattice_a0);
    vec_a[1][1] = -1./2.*(g->lattice_a0);
    
    vec_a[0][2] =  sqrt(3.)/2.*(g->lattice_a0);
    vec_a[1][2] = -1./2.*(g->lattice_a0);
    
    
    
    //b-vectors of the hexagonal lattice
    vec_b[0][0] = sqrt( 3. )*(g->lattice_a0);
    vec_b[1][0] = 0.;
    
    vec_b[0][1] = -sqrt(3.)/2.*(g->lattice_a0);
    vec_b[1][1] = +3./2.*(g->lattice_a0);
    
    vec_b[0][2] = -sqrt(3.)/2.*(g->lattice_a0);
    vec_b[1][2] = -3./2.*(g->lattice_a0);
    
    
    //K and K' points
    K[0] = 2.*dospi/sqrt(3.)/3./(g->lattice_a0) ;
    K[1] = 0. ;
    
    
    Kprime[0] = -2.*dospi/sqrt(3.)/3./(g->lattice_a0) ;
    Kprime[1] = 0. ;
    
    
    //Honeycomb lattice maxima momentum on x and y directions definition
    honeycomb_lattice_kmax[0] = dospi/sqrt(3.)/(g->lattice_a0);
    honeycomb_lattice_kmax[1] = dospi/3./(g->lattice_a0);
    

    ks00[0]= -1.7;
    ks00[1]= 0;

    ks01[0] = -1.7 + 2.*abs(g->k[0][0]);
    ks01[1] = 0;
    

    ks10[0] = 0.80 ;
    ks10[1] = 1.45;

    ks11[0] = 0.80 + 2.*g->k[0][0];
    ks11[1] = 1.45;
    
    ks12[0] = +0.80 ;
    ks12[1] = -1.45;
    
    ks13[0] = +0.80 + 2.*g->k[0][0];
    ks13[1] = -1.45;
    
    /*
     
     ks20 = {-1.1, 1.85};
     ks21 = {-1.1 + 2 kxMax, 1.85};
     ks22 = {-1.1, -1.85};
     ks23 = {-1.1 + 2 kxMax, -1.85};
    */
    

    ks20[0] = -1.10;
    ks20[1] =  1.85;
    
    ks21[0] =  -1.10 + 2.*abs(g->k[0][0]);
    ks21[1] =  1.85;
    
    ks22[0] = -1.10;
    ks22[1] = -1.85;
    
    ks23[0] = -1.10 + 2.*abs(g->k[0][0]);
    ks23[1] = -1.85;
    
    
    
}
//End lattice constant




//######################################
//New function
void solidstructure::Set_Of_B_BGrad( double const *kx, double const *ky )
{
    
    b0res0  = 0.; b1res0  = 0.;
    b2res0  = 0.; b3res0  = 0.;
    
    xb0der0 = 0.; xb1der0 = 0.;
    xb2der0 = 0.; xb3der0 = 0.;
    
    yb0der0 = 0.; yb1der0 = 0.;
    yb2der0 = 0.; yb3der0 = 0.;
    
    
    for ( n = 0; n < Nvects; n++ )
    {
        angle_a0 = (*kx)*vec_a[0][n] + (*ky)*vec_a[1][n];
        angle_b0 = (*kx)*vec_b[0][n] + (*ky)*vec_b[1][n];
        
        b0res0+= cos( angle_b0 );
        b1res0+= cos( angle_a0 );
        b2res0+= sin( angle_a0 );
        b3res0+= sin( angle_b0 );
        
        xb0der0+=  vec_b[0][n]*sin( angle_b0 );
        xb1der0+=  vec_a[0][n]*sin( angle_a0 );
        xb2der0+=  vec_a[0][n]*cos( angle_a0 );
        xb3der0+=  vec_b[0][n]*cos( angle_b0 );
        
        yb0der0+=  vec_b[1][n]*sin( angle_b0 );
        yb1der0+=  vec_a[1][n]*sin( angle_a0 );
        yb2der0+=  vec_a[1][n]*cos( angle_a0 );
        yb3der0+=  vec_b[1][n]*cos( angle_b0 );
        
    }
    
    b0res0*= 2. * t2 * cos( phi0 );
    b1res0*= t1;
    b2res0*= t1;
    b3res0*= -2.0 * t2 * sin( phi0 );
    b3res0+= M0;
    

    
    bnorm = sqrt(
                  b1res0 * b1res0
                  +b2res0 * b2res0
                  +b3res0 * b3res0
                  );
    
    
    xb0der0*= -2.*t2*cos( phi0 ) ;
    xb1der0*= -t1;
    xb2der0*= t1;
    xb3der0*= -2.*t2*sin( phi0 );
    
    
    
    yb0der0*= -2.*t2*cos( phi0 );
    yb1der0*= -t1;
    yb2der0*= t1;
    yb3der0*= -2.*t2*sin( phi0 );
    
    
    //Gauge choise
    if( gauge == 0 )
    {
        
        gauge_b_set[0] = b1res0; gauge_b_set[1] = b2res0; gauge_b_set[2] = b3res0;
        gauge_b_set_grad[0][0] = xb1der0; gauge_b_set_grad[0][1] = xb2der0; gauge_b_set_grad[0][2]= xb3der0;
        gauge_b_set_grad[1][0] = yb1der0; gauge_b_set_grad[1][1] = yb2der0; gauge_b_set_grad[1][2]= yb3der0;
        
    }
    
    if( gauge == 1 )
    {
        
        gauge_b_set[0] = b2res0; gauge_b_set[1] = b3res0; gauge_b_set[2] = b1res0;
        gauge_b_set_grad[0][0] = xb2der0; gauge_b_set_grad[0][1] = xb3der0; gauge_b_set_grad[0][2]= xb1der0;
        gauge_b_set_grad[1][0] = yb2der0; gauge_b_set_grad[1][1] = yb3der0; gauge_b_set_grad[1][2]= yb1der0;
        
    }
    
    if( gauge == 2 )
    {
        
        gauge_b_set[0] = b3res0; gauge_b_set[1] = b1res0; gauge_b_set[2] = b2res0;
        gauge_b_set_grad[0][0] = xb3der0; gauge_b_set_grad[0][1] = xb1der0; gauge_b_set_grad[0][2]= xb2der0;
        gauge_b_set_grad[1][0] = yb3der0; gauge_b_set_grad[1][1] = yb1der0; gauge_b_set_grad[1][2]= yb2der0;
        
    }
    
    
    
    //######################################
    //## Wavefunction normalization factor
    if ( rflag == 0 )
    {
        
        wave_norm  = sqrt( 2.*( bnorm * bnorm +  bnorm * gauge_b_set[2] ) );
        
    }
    
    if ( rflag == 1 )
    {
        
        wave_norm = 2./sqrt( 2. ) * ( bnorm + gauge_b_set[2]/2. - gauge_b_set[2] * gauge_b_set[2] / bnorm / 8. );
        
    }
    
    if ( rflag == 2 )
    {
        
        wave_norm  =      sqrt( 2.*( bnorm * bnorm + bnorm * gauge_b_set[2] ) ) +
        myRegularization( kx, ky, &ks00[0], &ks00[1]  ) + myRegularization( kx, ky, &ks01[0], &ks01[1]  )
        + myRegularization( kx, ky, &ks10[0], &ks10[1]  ) + myRegularization( kx, ky, &ks11[0], &ks11[1]  )
        + myRegularization( kx, ky, &ks12[0], &ks12[1]  ) + myRegularization( kx, ky, &ks13[0], &ks13[1]  )
        + myRegularization( kx, ky, &ks20[0], &ks20[1]  ) + myRegularization( kx, ky, &ks21[0], &ks21[1]  )
        + myRegularization( kx, ky, &ks22[0], &ks22[1]  ) + myRegularization( kx, ky, &ks23[0], &ks23[1]  );
        
    }
    

    
 

    //############################
    //##Wave-functions
    uwave_val0[0] = ( -gauge_b_set[0] + I*gauge_b_set[1] )/wave_norm;
    uwave_val0[1] = (  gauge_b_set[2] + bnorm  )/wave_norm;
    
    
    uwave_cond0[0] = ( gauge_b_set[2] + bnorm )/wave_norm;
    uwave_cond0[1] = ( gauge_b_set[0] + I*gauge_b_set[1] )/wave_norm;
    
    
    
    //#####################################
    // Gradient of BNorm = [ B1^2 + B2^2 +B3^2 ]^1/2
    // and
    // Gradient on the  wavefunction normalization factor
    for (n = 0; n<Ngrad; n++)
    {
        
        
        grad_bnorm[n]  = (gauge_b_set_grad[n][0] * gauge_b_set[0] + gauge_b_set_grad[n][1]*gauge_b_set[1] + gauge_b_set_grad[n][2] * gauge_b_set[2] )/bnorm;
        
        
       
        if ( rflag == 0 || rflag==2 )
        {
            
            w_norm_grad[n] = ( grad_bnorm[n] * ( 2. * bnorm + gauge_b_set[2] ) + bnorm * gauge_b_set_grad[n][2] )/wave_norm;
            
        }
        else
        {
            
            w_norm_grad[n] = 2./sqrt( 2. ) * ( grad_bnorm[n] + gauge_b_set_grad[n][2]/2.
                                          -  ( 2. *  bnorm * gauge_b_set[2] * gauge_b_set_grad[n][2]  - grad_bnorm[n] * gauge_b_set[2] * gauge_b_set[2] )/ bnorm/bnorm/8. );
        }
        
        //*/
    }
    //###########################################################


    
    //#####################################
    //###### WAVEFUNCTION GRADIENTs
    //#x-dir-grad
    xtemp[0]    = gauge_b_set_grad[0][0] * wave_norm - w_norm_grad[0] * gauge_b_set[0];
    xtemp[1]    = gauge_b_set_grad[0][1] * wave_norm - w_norm_grad[0] * gauge_b_set[1];
    xtemp[2]    = ( gauge_b_set_grad[0][2] + grad_bnorm[0] ) * wave_norm;
    xtemp[3]    = ( gauge_b_set[2] + bnorm ) * w_norm_grad[0];
    
    
    //#y-dir-grad
    ytemp[0]    = gauge_b_set_grad[1][0] * wave_norm - w_norm_grad[1] * gauge_b_set[0];
    ytemp[1]    = gauge_b_set_grad[1][1] * wave_norm - w_norm_grad[1] * gauge_b_set[1];
    ytemp[2]    = ( gauge_b_set_grad[1][2] + grad_bnorm[1] ) * wave_norm;
    ytemp[3]    = ( gauge_b_set[2] + bnorm ) * w_norm_grad[1];
    
    
    
    
    //wavefunction norm, its square
    wave_norm2  = wave_norm * wave_norm;
    
    
    
    
    //##################################
    //## W.F. VALENCE BAND GRAD
    wv_grad0[0][0]  = ( -xtemp[0] +  I*xtemp[1]  )/wave_norm2;
    wv_grad0[0][1]  = ( +xtemp[2] -   xtemp[3]   )/wave_norm2;
    
    wv_grad0[1][0]   = ( -ytemp[0] + I*ytemp[1]  )/wave_norm2;
    wv_grad0[1][1]   = ( +ytemp[2] -  ytemp[3]   )/wave_norm2;
    
    
    
    //#################################
    //## W.F. CONDUCTION BAND GRAD
    wc_grad0[0][0] = (  xtemp[2] -   xtemp[3]  )/wave_norm2;
    wc_grad0[0][1] = (  xtemp[0] + I*xtemp[1]  )/wave_norm2;
    
    wc_grad0[1][0] = (  ytemp[2]  - ytemp[3]   )/wave_norm2;
    wc_grad0[1][1] = (  ytemp[0] + I*ytemp[1]  )/wave_norm2;
    
    
    
    
    //##########################
    //## Structure Energy, Berry Connection
    //## and Dipole Matrix Element
    eg0         = 2. * bnorm;
    
    
    
    
    //Berry Connection difference: Xig0 = Xi_c - Xi_v = 2*Xi_c
    cdot( uwave_cond0, wc_grad0[0], Ngrad, &tempc[0] ) ;
    cdot( uwave_cond0, wc_grad0[1], Ngrad, &tempc[1] ) ;
    
    
    //Connection Difference
    chig0[0]    = 2.* real( I * tempc[0] ) ;
    chig0[1]    = 2.* real( I * tempc[1] ) ;
    
    
    
    
    
    
    //Dipole transiton val. -> cond. band
    cdot( uwave_cond0, wv_grad0[0], Ngrad, &tempc[0] );
    cdot( uwave_cond0, wv_grad0[1], Ngrad, &tempc[1] );
    
    
    dip_cv0[0]      =  I*tempc[0];
    dip_cv0[1]      =  I*tempc[1];
    
    
    velocity_cv0[0] = I * eg0 * dip_cv0[0];
    velocity_cv0[1] = I * eg0 * dip_cv0[1];
    
    zBcurva_c0      =  2. * imag( conj(dip_cv0[0])*dip_cv0[1] );
    
    
    
}



//Conduction Band energy/spectrum
void solidstructure::Ec( )
{
    
    energyc0 = b0res0 + bnorm;
    
}



//Valence Band energy/spectrum
void solidstructure::Ev( )
{
    
    energyv0 = b0res0 - bnorm;
    
}


void solidstructure::BerryConnectionCV()
{
    
    //Connection difference Xig0 = Xi_c - Xi_v = 2*Xi_c
    cdot( uwave_cond0, wc_grad0[0], Ngrad, &tempc[0] ) ;
    cdot( uwave_cond0, wc_grad0[1], Ngrad, &tempc[1] ) ;
    
    
    chi_c0[0]    = + real( I * tempc[0] ) ;
    chi_c0[1]    = + real( I * tempc[1] ) ;
    
    chi_v0[0]    = - real( I * tempc[0] ) ;
    chi_v0[1]    = - real( I * tempc[1] ) ;
    
    
}


void solidstructure::dipoleCV( )
{
    
    //Dipole transiton val. -> cond. band, dip_cv = i<u_c_k|\nabla_k|u_v_k>
    cdot( uwave_cond0, wv_grad0[0], Ngrad, &tempc[0] );
    cdot( uwave_cond0, wv_grad0[1], Ngrad, &tempc[1] );
    
    dip_cv0[0]  =  I*tempc[0];
    dip_cv0[1]  =  I*tempc[1];
    
}


//Cond. Band Berry Curvature: -((xGradtheta*yGradphi - xGradphi*yGradtheta)*Sin(theta))/2.
void solidstructure::zBerryCurvaCV( )
{
    
/*****************
 Before using this fungtion we should called:
 *******/
    
    cdot(  wc_grad0[0] , wc_grad0[1], Ngrad, &curvatemp1 );
    cdot(  wc_grad0[1] , wc_grad0[0], Ngrad, &curvatemp2 );
    
    zBcurva_c0 = real( I*(curvatemp1 - curvatemp2) );
    
    
    
}//*/



//##########################
//Group velocity
void solidstructure::group_velCV( )
{
    
    
    group_c0[0] =  xb0der0 + grad_bnorm[0] ;
    group_c0[1] =  yb0der0 + grad_bnorm[1] ;
    
    group_v0[0] =  xb0der0 - grad_bnorm[0] ;
    group_v0[1] =  yb0der0 - grad_bnorm[1] ;
    
    
}


//##########################
//Anomalous velocity
void solidstructure::anomalous_velCV( complex const efield, double const zBerryCurvature )
{
    
    
    xanomalous_c0   =  -imag( efield ) * zBerryCurvature;
    yanomalous_c0   =   real( efield ) * zBerryCurvature;
    
    
    xanomalous_v0   =  -imag( efield ) * (-zBerryCurvature);
    yanomalous_v0   =  +real( efield ) * (-zBerryCurvature);
    
    
}


//###########################
double solidstructure::ChernNumber( )
{
    
    double chern0=0.;
    int ktemp=0;

    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++ )
        {
            
            for ( i = 0; i < g->N[0]; i++ )
            {
            
                Set_Of_B_BGrad( &g->k[0][i], &g->k[1][j] );
               
                chern0+=zBcurva_c0*g->weight[ g->index(&i,&j,&ktemp) ]*g->dk[0]*g->dk[1]/4./pi;
            
            }
        
        }
    
    return chern0;
}



//##########################
void solidstructure::cdot( const complex *v1, const complex *v2, const int Nv, complex *res)
{
    
    *res = 0.;
    for ( int n = 0; n<Nv; n++)
        *res+= conj(v1[n])*v2[n];
    
}


//##########################
void solidstructure::energy_vc_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{

    if (flag==0)
    {
        
        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
        exit (EXIT_FAILURE);
        
    }
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
        for ( j = 0; j < g->N[1]/skip1; j++ )
        {
            fprintf( output, "\n");
            for ( i = 0; i < g->N[0]/skip0; i++ )
            {
                
                Set_Of_B_BGrad( &g->k[0][i*skip0], &g->k[1][j*skip1] );
                
                Ev();
                Ec();
                
                fprintf( output, "%e %e %.16e %.16e %.16e \n"
                        ,g->k[0][i*skip0]
                        ,g->k[1][j*skip1]
                        ,energyv0
                        ,energyc0
                        ,energyc0 - energyv0
                        );
                
            }

        }
    
    fflush(output);
    
}



//##########################
void solidstructure::dipole_cv_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    if (flag==0)
    {
        
        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
        exit (EXIT_FAILURE);
        
    }
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        
        fprintf( output, "\n");
        
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            ai = i*skip0;
            aj = j*skip1;
            al = l*skip2;
            
            
            fprintf( output, "%e %e %.16e %.16e %.16e %.16e \n"
                    ,g->k[0][i*skip0]
                    ,g->k[1][j*skip1]
                    ,real( dipole_cv[0][g->index( &ai, &aj, &al ) ] )
                    ,imag( dipole_cv[0][g->index( &ai, &aj, &al ) ] )
                    ,real( dipole_cv[1][g->index( &ai, &aj, &al ) ] )
                    ,imag( dipole_cv[1][g->index( &ai, &aj, &al ) ] )
                    );
            
        }
        
    }
    
    fflush(output);
}


//##########################
void solidstructure::connection_c_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
    if (flag==0)
    {
        
        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
        exit (EXIT_FAILURE);
        
    }
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        
        fprintf( output, "\n");
        
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            
            ai = i*skip0;
            aj = j*skip1;
            al = l*skip2;
            
            fprintf( output, "%e %e %e %e \n"
                    ,g->k[0][i*skip0]
                    ,g->k[1][j*skip1]
                    ,chi_gap[0][  g->index( &ai, &aj, &al ) ]/2.
                    ,chi_gap[1][  g->index( &ai, &aj, &al ) ]/2.
                    );
            
        }
        
    }
    
    fflush(output);
}


//##########################
void solidstructure::curvature_c_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    if (flag==0)
    {
        
        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
        exit (EXIT_FAILURE);
        
    }
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        
        fprintf( output, "\n");
        
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            ai = i*skip0;
            aj = j*skip1;
            al = l*skip2;
            
            
            fprintf( output, "%e %e %e \n"
                    ,g->k[0][i*skip0]
                    ,g->k[1][j*skip1]
                    ,zbcurva_c[ g->index( &ai, &aj, &al ) ]
                    );
            
            
        }
        
    }
    
    fflush(output);
    
}


//##########################
void solidstructure::group_vel_vc_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
    
    if (flag==0)
    {
        
        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
        exit (EXIT_FAILURE);
        
    }
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        
        fprintf( output, "\n");
        
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            //groupvel_c
            ai = i*skip0;
            aj = j*skip1;
            al = l*skip2;
            
            
            fprintf( output, "%e %e %e %e %e %e\n"
                    ,g->k[0][i*skip0]
                    ,g->k[1][j*skip1]
                    ,groupvel_v[0][  g->index(&ai, &aj, &al) ]
                    ,groupvel_v[1][  g->index(&ai, &aj, &al) ]
                    ,groupvel_c[0][  g->index(&ai, &aj, &al) ]
                    ,groupvel_c[1][  g->index(&ai, &aj, &al) ]
                    );
            
            
            
        }
        
    }
    
    fflush(output);
}






#endif /* solidstructure_h */
