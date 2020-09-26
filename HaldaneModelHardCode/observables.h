//
//  observables.h
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
//

#ifndef observables_h
#define observables_h

#include <stdlib.h>
#include <string>
//#include "mkl.h"
#include "constant.h"
#include "momaxis.h"
#include "solidstructure.h"
#include "laser.h"
#include <complex>
#define complex complex<double>

using namespace std;

class observables
{
    
private:
    int i, j, l, n;
    int it, jt, lt;
    
public:
    /**********************************************
     
     INPUT PARAMETERS, VARIABLES AND OBJECTS
     
     */
    
    momaxis *g;
    //laser *flaser;
    
    
    /**********************************************
     
     PROCESS PARAMETERS, VARIABLES AND OBJECTS
     
     */
    
    complex *coherence_pi0, *coherence_pi, *inhomo_pi, *homo_pi;
    complex *rabbi_omega;
    
    complex *occup_nv, *occup_nc, *inhomo_nc, *homo_nc;
    
    
    complex *occup_source;
    complex initial_nv0, initial_nc0, initial_pi0 ;
    complex temp;
    
    

    
    long int Ntime;
    
    
    /**********************************************
     
     OUTPUT PARAMETERS, VARIABLES AND OBJECTS
     
     */

    
    double interband_dipole[Ndim];     //x,y and z component, interband current
    double intraband_currents[Ndim];     //x,y and z component, intraband current
    double occup_nv0, occup_nc0;
    double occup_nv01, occup_nc01;
    
    double intg_nc0, intg_nv0;
    complex intg_pi0;



/*
 
 Outputs functions for,
 (1) pi(k,t) coherence
 (2) nc(k,t) occupation
 (3) time-dependent interband current
 (4) time-dependent intraband current
 (5) computing expectation values on for inter-dipole moment
    and intra-current
 
 */
    
    
    observables( momaxis *_g,  laser *_flaser );
    ~observables();
    void set_memory( );
    
    
    
    void cinitials_conditions( complex *f, complex init );
    void dinitials_conditions( double *f, double init );
    void set_inititals( complex init_nv0, complex init_nc0, complex init_pi0 );
    
    
    void RabbiOmega( const solidstructure *s, complex efield );
    void RabbiOmegaOccup( const complex *coh_pi0, const solidstructure *s, complex efield );
    
    
    
    void inter_dipole_expectation_value( solidstructure *s );
    void classical_velocity_expectation_value( solidstructure *s, complex efield );
    
    
    void integrate_occupations();
    void integrate_coherence();
    
    
    void coherence_output( FILE *coh_out, int skip0, int skip1, int skip2 );
    void gnuplot_coherence_output( FILE *coh_out, int skip0, int skip1, int skip2 );
    
    
    void cb_occup_output( FILE *coh_out, int skip0, int skip1, int skip2 );
    void gnuplot_cb_occup_output( FILE *coh_out, int skip0, int skip1, int skip2 );
    
    
    
    void intreband_dipole_moment_output( FILE *cfile, solidstructure *s, int *ktime, double *t, complex efield   );
    void intraband_current_output( FILE *cfile, solidstructure *s, int *ktime, double *t, complex efield  );
    
    
};



observables::observables( momaxis *_g,  laser *_flaser )
{
    
    g       = new momaxis( _g -> N, _g -> kmaxs );
    
    //flaser  = new laser( _flaser -> Npulses );
    //flaser -> copy_laser( _flaser );
    
    Ntime   = _flaser->Nmaxt;
    
    
    set_memory();
    
    
}



observables::~observables()
{
    
    
    delete g;
    //delete flaser;
    
    
    free( coherence_pi0 );
    free( coherence_pi );
    free( inhomo_pi );
    free( homo_pi );
    
    
    
    free( occup_nc );
    free( inhomo_nc );
    free( homo_nc );
    
    free( occup_source );
    
    free(rabbi_omega);
    //free( occup_nv );
    
    /*
    for ( n = 0; n < Ndim; n++ )
    {
        
        free(interband_dipole[n]);
        free(intraband_currents[n]);
        
    }*/
    
    
}




void  observables::set_inititals( complex init_nv0, complex init_nc0, complex init_pi0 )
{
    
    
    
    initial_nv0 = init_nv0;
    initial_nc0 = init_nc0;
    initial_pi0 = init_pi0;
    
    
    
    //cinitials_conditions( occup_nv, &initial_nv0 );
    cinitials_conditions( occup_nc, initial_nc0 );
    cinitials_conditions( coherence_pi, initial_pi0 );
    
    
    
}




void observables::set_memory( )
{
    
    
    
    coherence_pi0 = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( coherence_pi0,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    
    coherence_pi = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( coherence_pi,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    inhomo_pi = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( inhomo_pi,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    
    homo_pi = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( homo_pi,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    
    rabbi_omega = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( rabbi_omega,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    
    //occup_nv    = ( double* )malloc( ( g->Ntotal )*sizeof( double ) );
    //memset( occup_nv,     0,    sizeof( double )*(g->Ntotal) );
    
    
    
    occup_nc    = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( occup_nc,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    inhomo_nc    = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( inhomo_nc,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    homo_nc    = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( homo_nc,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    occup_source = ( complex* )malloc( ( g->Ntotal )*sizeof( complex ) );
    memset( occup_source,     0,    sizeof( complex )*(g->Ntotal) );
    
    
    interband_dipole[0] = 0.;
    interband_dipole[1] = 0.;
    
   // for ( n = 0; n < Ndim; n++ )
    //{
        
        
        //interband_dipole[n] = ( double* )malloc( ( Ntime )*sizeof( double ) );
        //memset( interband_dipole[n],     0,    sizeof( double )*( Ntime ) );
        
        
        //intraband_currents[n] = ( double* )malloc( ( Ntime )*sizeof( double ) );
        //memset( intraband_currents[n],     0,    sizeof( double )*( Ntime ) );
        
        
    //}
    
    
    
}



void observables::RabbiOmega( const solidstructure *s, complex efield )
{
    
    
    
    for ( n = 0; n < g->Ntotal; n++ )
    {
        

        rabbi_omega[ n ] =    s->dipole_cv[0][ n ]*real( efield )
                            + s->dipole_cv[1][ n ]*imag( efield );
        
        
    }
    
    
}





void observables::RabbiOmegaOccup( const complex *coh_pi0, const solidstructure *s, complex efield )
{
    
    
    for ( n = 0; n < g->Ntotal; n++ )
    {
        
    
        
      
        
        
        occup_source[ n ] =  -2.*imag( ( conj( s->dipole_cv[0][ n ] )*real( efield )
                                        +  conj( s->dipole_cv[1][ n ] )*imag( efield )
                                        )*coh_pi0[ n ]
                                      );
        

        
        
        
    }
    
    
}








void observables::cinitials_conditions( complex *f, complex init )
{
    
    
    for( i = 0; i < g -> Ntotal; i++ )
            f[i] = init;
    // g->index(&i,&j,&l)
    
    
}




void observables::dinitials_conditions( double *f, double init )
{


    for( i = 0; i < g -> Ntotal; i++ )
        f[i] = init;
    // g->index(&i,&j,&l)


}



void observables::integrate_occupations()
{
    
    
   intg_nc0 = 0.;
   //intg_nv0 = 0.;
    
   for ( n = 0; n < g->Ntotal; n++ )
    {
            
        intg_nc0+= real( occup_nc[n] );
        //intg_nv0+= 1. - real( occup_nc[n] );
        
    }
    
    intg_nc0*= ( (g->dV)/double(g->Ntotal) );
    intg_nv0 = 1. - intg_nc0;
    
    //intg_nv0*= ( (g->dV)/double(g->Ntotal) );
    
}


void observables::integrate_coherence()
{
    
    
    intg_pi0 = complex( 0., 0. );
    
    for ( n = 0; n < g->Ntotal; n++ )
    {
        
        intg_pi0+= coherence_pi[n] ;
        
    }
    
    intg_pi0*= (g->dV/double(g->Ntotal));

    
}




void observables::inter_dipole_expectation_value( solidstructure *s  )
{
    
    
    
    interband_dipole[0] = 0.;
    interband_dipole[1] = 0.;
    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++ )
            for ( i = 0; i < g->N[0]; i++)
            {
                
                
                //x-direction dipole expectation value
                interband_dipole[0]+=  real( coherence_pi[ g->index( &i, &j, &l ) ]
                                                    * conj( s->dipole_cv[0][ g->index( &i, &j, &l ) ] ) );
                
                
                
                //y-direction dipole expectation value
                interband_dipole[1]+=  real( coherence_pi[ g->index( &i, &j, &l ) ]
                                                * conj( s->dipole_cv[1][ g->index( &i, &j, &l ) ] ) );
                
                
                
            }
    
    
    
    interband_dipole[0]*= 2.*charge_electron_au * g->dV;
    interband_dipole[1]*= 2.*charge_electron_au * g->dV;
    
    
    
    //interband_dipole[0];
    //interband_dipole[1];
    
    
    
    
}



void observables::classical_velocity_expectation_value( solidstructure *s, complex efield )
{
    
    /********************************************
     
        vclassical = vg_c * n_c   + vg_v * n_v
                     + va_c * n_c + va_v * n_v
     
     *****/
    
    
    intraband_currents[0]=0.;
    intraband_currents[1]=0.;
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++ )
            for ( i = 0; i < g->N[0]; i++ )
            {
                
                
                occup_nc0 = real( occup_nc[ g->index( &i, &j, &l ) ] );
                occup_nv0 = 1. - occup_nc0;
                
                
                s->anomalous_velV( efield, s->zbcurva_v[ g->index(&i,&j,&l) ] );
                s->anomalous_velC( efield, s->zbcurva_c[ g->index(&i,&j,&l) ] );
                
                
                
                intraband_currents[0]+=    occup_nv0 * ( s->groupvel_v[0][ g->index(&i,&j,&l) ] + s->xanomalous_v0 )
                                                 + occup_nc0 * ( s->groupvel_c[0][ g->index(&i,&j,&l) ] + s->xanomalous_c0 );
                
                
                intraband_currents[1]+=    occup_nv0 * ( s->groupvel_v[1][ g->index(&i,&j,&l) ] + s->yanomalous_v0 )
                                                 + occup_nc0 * ( s->groupvel_c[1][ g->index(&i,&j,&l) ] + s->yanomalous_c0 );
                
                
            }
    
    
    intraband_currents[0]*= charge_electron_au * g->dV;
    intraband_currents[1]*= charge_electron_au * g->dV;
    
    
    //intraband_currents[0]++;
    //intraband_currents[1]++;
    
    
}




void observables::intraband_current_output( FILE *cfile, solidstructure *s, int *ktime, double *t, complex efield  )
{
    

    
    classical_velocity_expectation_value( s, efield );
    
    
    
    fprintf( cfile, " %.5d    %e    %e    %e    %.16e     %.16e  \n",
            *ktime, *t, real( efield ), imag( efield )
            ,intraband_currents[0]
            ,intraband_currents[1]
            );
    
    
    
    fflush( cfile );
                                     
                                     

}


void observables::intreband_dipole_moment_output( FILE *cfile, solidstructure *s, int *ktime, double *t, complex efield  )
{
    
    
    
    //inter_dipole_expectation_value()
    
    inter_dipole_expectation_value( s  );
    
    
    
    fprintf( cfile, " %.5d    %e    %e    %e    %.16e     %.16e  \n",
            *ktime, *t, real( efield ), imag( efield )
            ,interband_dipole[0]
            ,interband_dipole[1]
            );
    
    
    
    fflush( cfile );
    
    
    
}





void observables::gnuplot_coherence_output( FILE *coh_out, int skip0=1, int skip1=1, int skip2=1 )
{


    for ( l = 0; l < g->N[2]/skip2; l++ )
        for ( j = 0; j < g->N[1]/skip1; j++ )
        {
            fprintf( coh_out, "\n" );
            for ( i = 0; i < g->N[0]/skip0; i++ )
            {
                
                
                it = i*skip0;
                jt = j*skip1;
                lt = l*skip2;
                
            
                fprintf( coh_out, "%e %e %.16e %.16e\n"
                        ,g->k[0][it]
                        ,g->k[1][jt]
                        ,real( coherence_pi[ g->index( &it, &jt, &lt ) ] )
                        ,imag( coherence_pi[ g->index( &it, &jt, &lt ) ] )
                    );
        
            }
        }
    
}


void observables::coherence_output( FILE *coh_out, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
  for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            it = i*skip0;
            jt = j*skip1;
            lt = l*skip2;
            
            
            fprintf( coh_out, "%e %e\n"
                    ,real( coherence_pi[ g->index( &it, &jt, &lt ) ] )
                    ,imag( coherence_pi[ g->index( &it, &jt, &lt ) ] )
                    );
            
        }
    }
    
    
}



void observables::cb_occup_output( FILE *coh_out, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
    for ( l = 0; l < g->N[2]/skip2; l++ )
        for ( j = 0; j < g->N[1]/skip1; j++ )
        {
            for ( i = 0; i < g->N[0]/skip0; i++ )
            {
            
            
                it = i*skip0;
                jt = j*skip1;
                lt = l*skip2;
            
            
                fprintf( coh_out, "%.16e \n"
                    ,real( occup_nc[ g->index( &it, &jt, &lt ) ] )
                    );
            
            }
        }
    
}


void observables::gnuplot_cb_occup_output( FILE *coh_out, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
 for ( l = 0; l < g->N[2]/skip2; l++ )
    for ( j = 0; j < g->N[1]/skip1; j++ )
    {
        fprintf( coh_out, "\n" );
        
        for ( i = 0; i < g->N[0]/skip0; i++ )
        {
            
            
            it = i*skip0;
            jt = j*skip1;
            lt = l*skip2;
            
            
            fprintf( coh_out, "%e %e %.16e %.16e\n"
                    ,g->k[0][it]
                    ,g->k[1][jt]
                    ,real( occup_nc[ g->index( &it, &jt, &lt ) ] )
                    ,imag( occup_source[ g->index( &it, &jt, &lt ) ] )
                    );
            
        }
    }
    
}



#endif /* observables_h */
