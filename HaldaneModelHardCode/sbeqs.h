//
//  blocheqs.h
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
//

#ifndef sbeqs_h
#define sbeqs_h
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <iomanip>
//#include "mpi.h"
//#include "mkl.h"



#include "constant.h"
#include "timegrid.h"
#include "timeobject.h"
#include "laser.h"
#include "momaxis.h"
#include "operators.h"
//#include "solidstructure.h"
//#include "observables.h"


using namespace std;

class sbeqs
{
    
public:

    
    /* INPUTs
     
     Derivative operator
     
     */
    momaxis *g;
    //laser *afield;



    int Nmaxt, n;
    double dt;
    double halfdt;

    complex ibeta;
    complex hbeta;
    
    complex compdt;
    
    
    
    /*
     
        Unitary Occupation propagator
     
     */
    
    
    /*
     
        Initials conditions for coherence, pi,
        valence population nc,
        and
        conduction population
     
     
     */
    
    
    sbeqs( momaxis *_g, laser *_ef );
    ~sbeqs(  );
    
    
    void coherence_time_mom_evol_SBEs( const solidstructure *crystal
                                            ,observables *qobs
                                            ,operators *acts
                                        ,complex efield2
                                        ,complex efield4
                                      );
    
    
    
    void occupation_time_mom_evol_SBEs( const solidstructure *crystal
                                             ,observables *qobs
                                             ,operators *acts
                                        ,complex efield2
                                        ,complex efield4
                                             );
    
    
    
};

sbeqs::sbeqs(  momaxis *_g, laser *_ef  )
{
    
    g       = new momaxis( _g-> N, _g -> kmaxs );
    
    //afield  = new laser( _ef->Npulses );
    
    Nmaxt       = _ef->Nmaxt;
    
    dt          = _ef->dt;
    compdt      = dt;
    halfdt      = dt;
    

    ibeta       = -I*dt;
    hbeta       = 1.;

    
}

sbeqs::~sbeqs()
{
    
    
    delete g;
    
    
}

/**********************************************************
 
 Time and Momentum-Space Propagations of the Coherence
 
 *****/

//InHomogeneous solver for the coherence of the SBEs
//memcpy( qobs->coherence_pi0, qobs->coherence_pi, ( gaxis.Ntotal )*sizeof( complex ) );
void sbeqs::coherence_time_mom_evol_SBEs( const solidstructure *crystal
                                  ,observables *qobs
                                  ,operators *acts
                                  ,complex efield2
                                  ,complex efield4
                                  )
{
    
    
    
    //Evaluating Rabbi Frequency Omega = d_cv * E(t)
    qobs->RabbiOmega( crystal, efield2 );
    
    memcpy( qobs->inhomo_pi, qobs->rabbi_omega, sizeof(complex)*g->Ntotal );
    
    qobs->RabbiOmega( crystal, efield4 );

    
    
    //Homogeneous solver for the coher.
    acts->forward_prop( qobs->coherence_pi
                       ,crystal
                       ,efield4
                       ,dt
                       ,hbeta
                       );
    
    
    
    for ( n = 0; n < g->Ntotal; n++ )
            qobs->coherence_pi[n] =  qobs->coherence_pi[n] + ( qobs->rabbi_omega[n] + qobs->inhomo_pi[n] )*ibeta;
    
    
    
    
    acts->backward_prop(     qobs->homo_pi //homogeneous solution
                            ,qobs->coherence_pi
                            ,crystal
                            ,efield4
                            ,dt
                            );//*/
    
    
}

/*  // forward prop of Omega at the inhomogenous part of the SBEs
 acts->forward_prop( qobs->rabbi_omega //rabbi-feq = dcv*E(t) and return der. times ibeta
 ,crystal
 ,efield4
 ,halfdt
 ,ibeta
 );
 
 
 acts->backward_prop(     qobs->inhomo_pi //Inhomogenoeus solution
 ,qobs->rabbi_omega
 ,crystal
 ,efield4
 ,halfdt
 );
 //*/

//for ( n = 0; n < g->Ntotal; n++ )
//    qobs->coherence_pi[n] = qobs->homo_pi[n]


//Solution for coherence propagation
/****************************************************************/
//acts->solution( qobs->coherence_pi, qobs->inhomo_pi, qobs->homo_pi );


/****************************************************************
 
 Time and Momentum-Space Propagations of the Occupations
 
 *****/


void sbeqs::occupation_time_mom_evol_SBEs( const solidstructure *crystal
                                        ,observables *qobs
                                        ,operators *acts
                                        ,complex efield2
                                        ,complex efield4
                                        )
{
    
    
    


    /******************************************************
     Inhomogenues solver for occupations
    Computing rabby freq. times the coherence at t_ktime.
    */
    qobs->RabbiOmegaOccup( qobs->coherence_pi0, crystal, efield2);
    
    
    
    
    /****************************************
     
     Occupation propagation, cond. band.
     computing coherence pi at t_ktime + dt/2.
     
     ***/

    acts->linear_interpolation( qobs->coherence_pi0, qobs->coherence_pi );
    
    
    
    
    
    /******************************************************
     Inhomogenues solver for occupations
     Computing rabby freq. times the coherence at t_ktime + dt/2.
     */
    
    
    memcpy( qobs->inhomo_nc, qobs->occup_source, sizeof(complex)*g->Ntotal );
    
    
    qobs->RabbiOmegaOccup( qobs->coherence_pi0, crystal, efield4 );
    

    
    /******************************************************
 
     Homogeneous solver for the occupations of the SBEs
 
    */

    //Homogenoeus solver for the occup.
    acts->forward_prop_occup(
                             qobs->occup_nc
                             ,efield4
                             ,dt
                             ,hbeta
                             );//*/
    
    
    for ( n = 0; n < g->Ntotal; n++ )
        qobs->occup_nc[n] =  qobs->occup_nc[n] + ( qobs->inhomo_nc[n] + qobs->occup_source[n] )*dt;
    
    
   
    acts->backward_prop_occup( qobs->homo_nc
                              ,qobs->occup_nc
                              ,efield4
                              ,dt
                              );//*/
    
    
    //Solution for occup. propagation
    //acts->solution( qobs->occup_nc, qobs->inhomo_nc, qobs->homo_nc );


}


/******************************************************
 acts->forward_prop_occup(
 qobs->occup_source //input-source and derivative computed + potential op
 ,efield4
 ,halfdt
 ,compdt
 );
 
 acts->backward_prop_occup( qobs->inhomo_nc
 ,qobs->occup_source
 ,efield4
 ,halfdt
 );
 
 //*/




#endif /* semiconductor bloch eqs_h */
