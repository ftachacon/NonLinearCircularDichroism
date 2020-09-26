//
//  operators.h
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/20/19.
//

#ifndef operators_h
#define operators_h

#include <stdlib.h>
#include <string>
//#include "mkl.h"
#include <complex>
#define complex complex<double>


#include "constant.h"
#include "momaxis.h"
#include "solidstructure.h"
#include "laser.h"
#include "observables.h"


using namespace std;


class operators
{
 
//private:

    
    
public:
    /*
 
        INPUTS PARAMETERS AND OBJECTS
 
    */
    
    momaxis *g;
    
    
    int i,j,l,n;
    int iint, iend;
    int jint, jend;
    
    
    /*
     
        PROCESS PARAMETERS AND OBJECTS
     
     */
    
    complex **a, **b, **c; //vectors for derivative on x and y directions
    complex **fder, *gam[Ngrad];
    
    
    
    //double T2;
    double forward_dt, backward_dt;
    
    
    
    complex unitary_pot0;
    
    

    complex fint, fend;
    complex alpha, ralpha,lalpha;
    
    
    

    
    
    
    
    operators( momaxis *_g );
    ~operators(  );
    
    
    void set_memory( );
    
    
    
    
    /*
     
     Derivative operator
     
     */
    
    
    void xdiff( const  complex up_y, const  complex down_y, complex *res );
    complex xdiff_v2( const  complex up_y, const  complex down_y );
    void ydiff( const  complex up_y, const  complex down_y, complex *res );
    
    
    void x_full_der_op( complex alp,   complex *f  );
    void y_full_der_op( complex alp,   complex *f  );
    
    
    
    void x_der_evaluation(complex alp, const  complex *f  );
    void y_der_evaluation(complex alp, const  complex *f  );
    
    
    

 
    
    
    
    
    
    void pot_unitary_propagator(  double egap
                                 ,double xchigap
                                 ,double ychigap
                                 ,double T2
                                 ,complex efield
                                 ,double dt
                                );

    
    
    
    //Construction of pseudo kinetic operator
    void xkinetic_unitary_propagator( const  complex *f
                                      ,complex efield
                                      ,double dt
                                     );
    
    
    

    
    void forward_prop(   complex *f
                               ,const solidstructure *s
                               ,complex efield
                               ,double dt
                               ,complex beta
                             );
    
    
    
    
    
    void backward_prop(   complex *u
                                ,complex *ra
                                ,const solidstructure *s
                                ,complex efield
                                ,double dt
                                );

    
    
    
    
    void forward_prop_occup( complex *f
                                    ,complex efield
                                    ,double dt
                                    ,complex beta
                                    );
    
    
    
    
    void backward_prop_occup( complex *u
                                ,complex *ra
                                ,complex efield
                                ,double dt
                                );
    
    
    
    
    void solution( complex *fsol
                  ,complex *ug
                  ,complex *homog );
    
    
    
    void copy_density_matrix(complex *dest, const complex *source);
    
    
    void nrerror( char error_text[] );
    
    void tridag( complex *a
                , complex *b
                , complex *c
                ,complex *r
                ,complex *u
                , complex *gam1
                , int *Nq );
    
    void linear_interpolation( complex *inf, const complex *f );
    
};


//Initalization of operators and construsctor
operators::operators( momaxis *_g )
{
 
    g = new momaxis( _g-> N, _g-> kmaxs );
    
    
    set_memory( );
    
}


operators::~operators()
{

    
    for ( i = 0; i < Ndim-1; i++ )
    {
        
        free( a[i] );
        free( b[i] );
        free( c[i] );
        
        
        free( fder[i] );
        free( gam[i] );
        
    }
    
    
    
    free( a );
    free( b );
    free( c );
    free( fder );
    
    delete g;

    
}


void  operators::set_memory( )
{
    
    
    
    a   = ( complex** )malloc( ( Ndim-1 )*sizeof( complex* ) );
    b   = ( complex** )malloc( ( Ndim-1 )*sizeof( complex* ) );
    c   = ( complex** )malloc( ( Ndim-1 )*sizeof( complex* ) );
    
    fder = ( complex** )malloc( ( Ngrad ) *sizeof( complex* ) );
    
    
    for ( i = 0; i < Ngrad; i++ )
    {
        
        
        a[i] = ( complex* )malloc( g->N[i] * sizeof( complex ) );
        b[i] = ( complex* )malloc( g->N[i] * sizeof( complex ) );
        c[i] = ( complex* )malloc( g->N[i] * sizeof( complex ) );
        
        
        memset( a[i],     0,    sizeof( complex )*(g->N[i]) );
        memset( b[i],     0,    sizeof( complex )*(g->N[i]) );
        memset( c[i],     0,    sizeof( complex )*(g->N[i]) );
        
        
        gam[i]   = ( complex* )malloc( g->N[i] * sizeof( complex ) );
        memset( gam[i],     0,    sizeof( complex )*(g->N[i]) );
        
        
        fder[i] = ( complex* )malloc( g->N[i] * sizeof( complex ) );
        memset( fder[i],     0,    sizeof( complex )*(g->N[i]) );
        
        
    }
  
    
}



void operators::xdiff( const complex up_y, const complex down_y, complex *res )
{
    
    *res = (up_y - down_y)/2./( g->dk[0] );
    
}


complex operators::xdiff_v2( const complex up_y, const complex down_y)
{
    
    return (up_y - down_y)/2./( g->dk[0] );
    
}



void operators::ydiff( complex up_y, complex down_y, complex *res )
{
    
    *res = (up_y - down_y)/2./( g->dk[1] );
    
}


/*derivative operator x-direction*/
void operators::x_full_der_op( complex alp, complex *f )
{
    
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j= 0; j < g->N[1]; j++  )
        {
            
            
            x_der_evaluation( alp, f );
            
            iint = 0;
            memcpy(  &f[ g->index( &(iint), &j, &l ) ]
                    ,&fder[0][iint]
                    , ( g->N[0] ) * sizeof( complex )
                    );
            
            
        }
    

}




void operators::x_der_evaluation( complex alp, const complex *f )
{
    
    iint = 0;
    iend = 1;
    
    
    xdiff( f[ g->index( &(iend), &j, &l ) ]
          ,f[ g->index( &(iint), &j, &l ) ]
          ,&fder[0][iint]
          );
    
    
    fder[0][iint]*= alp;
    
    
    
    for ( i = 1; i < g->N[0]-1; i++ )
    {
        
        iint = i-1;
        iend = i+1;
        
        
        xdiff(  f[ g->index( &(iend), &j, &l ) ]
               ,f[ g->index( &(iint), &j, &l ) ]
               ,&fder[0][i]
              );
        
        fder[0][i]*= alp;
        
    }
    
    
    
    iint = g->N[0]-2;
    iend = g->N[0]-1;
    
    
    
    xdiff( f[ g->index( &(iend), &j, &l ) ]
          ,f[ g->index( &(iint), &j, &l ) ]
          ,&fder[0][iend]
          );
    
    
    fder[0][iend]*= alp;
        
}



/*derivative operator y-direction*/
void operators::y_full_der_op( complex alp, complex *f )
{
    
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( i= 0; i < g->N[0]; i++  )
            {
                
                y_der_evaluation( alp, f );
                
                for ( j = 0; j < g->N[1]; j++ )
                    f[ g->index( &i, &j, &l ) ] = fder[1][j];
                
            }
    
}



void operators::y_der_evaluation( complex alp, const complex *f )
{
    
    
    jint = 0;
    jend = 1;
    
    
    ydiff( f[ g->index( &i, &(jend), &l ) ]
          ,f[ g->index( &i, &(jint), &l ) ]
          ,&fder[1][jint]
          );
    
    
    fder[1][jint]*= alp;
    
    
    for ( j = 1; j < g->N[1]-1; j++ )
    {
        
        jint = j-1;
        jend = j+1;
        
        ydiff(  f[ g->index( &i, &(jend), &l ) ]
              ,f[ g->index( &i, &(jint), &l ) ]
              ,&fder[1][j]
              );
        
        fder[1][j]*= alp;
        
    }
    
    
    jint = g->N[1]-2;
    jend = g->N[1]-1;
    
    
    ydiff(  f[ g->index( &i, &(jend), &l ) ]
          ,f[ g->index( &i, &(jint), &l ) ]
          ,&fder[1][jend]
          );
    
    fder[1][jend]*= alp;
    
    jint = 0;
    
}




void operators::xkinetic_unitary_propagator(  const  complex *f, complex efield, double dt )
{
    
    /************************Kinetic********************
     
        kinetic operator UnitaryOperator ()
     
        xfder =  - I*dt/4.*( I*Efied(t) der )*f
     
        */

    
    alpha = ( -I*( dt )/2. )*( I*real( efield ) );
    
    //cout << "\n alpha =" << alpha ;
    //cout << "; j = " << j << " and l = " << l << endl;
    //x_der_evaluation( alpha, f );
    
    
    iint = 0;
    iend = 1;
    
    //cout << "\n f[0] " << f[ g->index( &(iend), &j, &l )];
    
    fder[0][iint] = xdiff_v2( f[ g->index( &(iend), &j, &l ) ]
                             ,f[ g->index( &(iint), &j, &l ) ]
                             )*alpha;
    
    
    for ( i = 1; i < g->N[0]-1; i++ )
    {
        
        iint = i-1;
        iend = i+1;
        
        
        fder[0][i] = xdiff_v2(  f[ g->index( &(iend), &j, &l ) ]
                                ,f[ g->index( &(iint), &j, &l ) ]
                                )*alpha;
        
    }
    
    
    
    iint = g->N[0]-2;
    iend = g->N[0]-1;
    
    
    
    fder[0][iend] = xdiff_v2( f[ g->index( &(iend), &j, &l ) ]
                             ,f[ g->index( &(iint), &j, &l ) ]
                             )*alpha;
    

    
}



void operators::pot_unitary_propagator( double egap
                                       ,double xchigap
                                       ,double ychigap
                                       ,double T2
                                       ,complex efield
                                       ,double dt )
{
    
    /*********************************
     Unitary Pot Operator(f)
     
        fr = Identity - I*dt/2.*Up
     
     //*/
    
    /*unitary_pot0 = 1. - I*( *dt )/2.*(  *egap
                                            + *xchigap * real( efield )
                                            + *ychigap * imag( efield )
                                            - I/(*T2)
                                            );
    //*/
    
    unitary_pot0 = 1. -  I*( dt )/2.*(  egap
                                      + xchigap * real( efield )
                                      + ychigap * imag( efield )
                                     - I/T2
                                      );
    
    
}


/*void operators::full_pseudokinetic_OP( complex *efield, complex *f )
{
 
}*/





void operators::forward_prop( complex *f
                         ,const solidstructure *s
                         ,complex efield
                         ,double dt
                         ,complex beta )
{
    
    /****************************
     Full UnitaryOperator
     
     
        fr = I*dt*(Identity - I*dt/2.( I*Efied(t) der + Up) ) f
     
     
     */
    
    forward_dt  = dt;
    
    //cout << "\n\nFrom Operators dt =" << forward_dt << endl ;
    
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++  )
        {
            
            
            xkinetic_unitary_propagator( f, efield, forward_dt );
            
            
            for ( i = 0; i < g->N[0]; i++  )
            {
                
                
                    pot_unitary_propagator( s->energy_gap[ g->index(&i,&j,&l) ]
                                           ,s->chi_gap[0][ g->index(&i,&j,&l) ]
                                           ,s->chi_gap[1][ g->index(&i,&j,&l) ]
                                           ,s->T2
                                           ,efield
                                           ,forward_dt
                                        );
                

                f[g->index(&i,&j,&l)] = ( fder[0][i] + unitary_pot0*f[g->index(&i,&j,&l)] )*beta;

                
                
            }
            
            
        }
    
}




void operators::backward_prop( complex *u
                                     ,complex *ra
                                     , const solidstructure *s
                                     ,complex efield
                                     ,double dt )
{
    
    /************************************************************************

     Full solution to the inhomogenous and homogeous part of SBEs
     
        [Identity + I*dt/2.( I*Efied(t) der + Up) ] u = ra*alpha
     
     */
    
    backward_dt  = -dt;
    lalpha      = ( -I*backward_dt/g->dk[0]/4. ) *( I*real( efield ) );
    
    //alpha = ( -I*( dt )/2. )*( I*real( efield ) );
    
    
    for ( l = 0; l < g->N[2]; l++ )
      for ( j = 0; j < g->N[1]; j++  )
       {
          
          
          i = 0;
          
          
          pot_unitary_propagator(  s->energy_gap[ g->index(&i,&j,&l) ]
                                 ,s->chi_gap[0][ g->index(&i,&j,&l) ]
                                 ,s->chi_gap[1][ g->index(&i,&j,&l) ]
                                 ,s->T2
                                 ,efield
                                 ,backward_dt
                                 );
          
          
          
          a[0][i]    = - lalpha;
          b[0][i]    = - lalpha + unitary_pot0;
          c[0][i]    =   lalpha;
          
          
          for ( i = 1; i < g->N[0]-1; i++  )
          {
            
              
              
              pot_unitary_propagator( s->energy_gap[ g->index(&i,&j,&l) ]
                                     ,s->chi_gap[0][ g->index(&i,&j,&l) ]
                                     ,s->chi_gap[1][ g->index(&i,&j,&l) ]
                                     ,s->T2
                                     ,efield
                                     ,backward_dt
                                   );
            
              
              a[0][i]    = - lalpha;
              b[0][i]    =   unitary_pot0;
              c[0][i]    =   lalpha;
              
            
          }
        
          
          i = g->N[0]-1;
          
          
          pot_unitary_propagator( s->energy_gap[ g->index(&i,&j,&l) ]
                                 ,s->chi_gap[0][ g->index(&i,&j,&l) ]
                                 ,s->chi_gap[1][ g->index(&i,&j,&l) ]
                                 ,s->T2
                                 ,efield
                                 ,backward_dt
                                 );
          
          
          a[0][i]    = - lalpha;
          b[0][i]    =   lalpha +  unitary_pot0;
          c[0][i]    =   lalpha;
          
        
          i = 0;
           
          tridag( a[0]
                 ,b[0]
                 ,c[0]
                 ,&ra[  g->index(&i,&j,&l) ]
                 ,&u[ g->index(&i,&j,&l) ]
                 ,gam[0]
                 ,&g->N[0]
                 );
           
           memcpy( &ra[  g->index(&i,&j,&l) ], &u[ g->index(&i,&j,&l) ], sizeof(complex)*(g->N[0]) );
           
     }
    
}




void operators::solution( complex *fsol, complex *ug, complex *homog )
{

    
    for ( n = 0; n < g->Ntotal; n++ )
        fsol[n] = ug[ n ] + homog[ n ];
    
    
}




/*
 
 
 
 Operators for evolution of the occupations
 
 
 
 */



void operators::forward_prop_occup( complex *f
                                    ,complex efield
                                    ,double dt
                                    ,complex beta)
{
    
    /**********************************
     
     Full UnitaryOperator
        fr = dt*[ Identity - I*dt/2.( I*Efied(t) der ) ] f
     
     
     */
    
    forward_dt  = dt;

    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j= 0; j < g->N[1]; j++  )
        {
        
        
            
            xkinetic_unitary_propagator( f, efield, forward_dt );
        
        
            for ( i = 0; i < g->N[0]; i++  )
                f[g->index(&i,&j,&l)] = ( f[g->index(&i,&j,&l)] + fder[0][i] )*beta;
            
            
        
        }
    
}


void operators::copy_density_matrix(complex *dest, const complex *source)
{
    for ( n = 0; n < g->Ntotal; n++ )
        dest[n] = source[n];
}



void operators::backward_prop_occup(  complex *u
                                     ,complex *ra
                                     ,complex efield
                                     ,double dt
                                           )
{
    
    /*
     
     Full solution to
     {
     
        [Identity + I*dt/2.( I*Efied(t) der ) ] u = ra*alpha
     
     }
     
     */
    
    backward_dt  = -dt;
    lalpha =  ( -I*backward_dt/4.  )*( I*real( efield )/g->dk[0] );
    
    
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j= 0; j < g->N[1]; j++  )
        {
        
            i = 0;
        
        
            a[0][i]    = - lalpha;
            b[0][i]    = 1. - lalpha;
            c[0][i]    =  lalpha;
        
        
            for ( i = 1; i < g->N[0]-1; i++  )
            {
                
                a[0][i]    = - lalpha;
                b[0][i]    =   1.;
                c[0][i]    =  lalpha;
                
            }
        
        
            i = g->N[0]-1;
            
        
            a[0][i]    = - lalpha;
            b[0][i]    =   1. + lalpha ;
            c[0][i]    =  lalpha;
        
        
            i = 0;
        
            tridag( a[0]
                   ,b[0]
                   ,c[0]
                   ,&ra[  g->index(&i,&j,&l) ]
                   ,&u[ g->index(&i,&j,&l) ]
                   ,gam[0]
                   ,&g->N[0]
                   );
        
            
           memcpy( &ra[  g->index(&i,&j,&l) ], &u[ g->index(&i,&j,&l) ], sizeof(complex)*(g->N[0]) );
            
            
        
        }
    
}






//Tridiag solver
void operators::tridag( complex *a
                       ,complex *b
                       ,complex *c
                       ,complex *r
                       ,complex *u
                       ,complex *gam1
                       ,int *Nq)
{
    
    //a down diagonal
    //b main diagonal
    //c upper diagonal
    
    complex bet=b[0]; // Declare and define auxiliar array
    
    if(b[0]==0.0) nrerror ("error 1 en tridag");
    
    u[0] = r[0]/bet;
    
    
    for( n = 1; n < *Nq; n++ )
    {
        
        gam1[n] = c[n-1]/bet;
        
        bet     = b[n] - a[n]*gam1[n];
        
        if (bet==0.0) nrerror("error 2 en tridag");
        
        u[n]  = ( r[n] - a[n]*u[n-1] )/bet;
    }
    
    
    for( n = (*Nq-2); n >=0; n--)
        u[n]-= gam1[n+1]*u[n+1];
    
    
}




void operators::nrerror (char error_text[])
{
    
    
    fprintf(stderr, "Runtime error..\n");
    
    fprintf(stderr, "%s\n", error_text);
    
    fprintf(stderr,"... now exiting system..\n");
    
    
}




void operators::linear_interpolation( complex *inf, const complex *f )
{
    
    
    for ( l = 0; l < g->N[2]; l++ )
        for ( j = 0; j < g->N[1]; j++  )
            for ( i = 0; i < g->N[0]; i++  )
            {
                
                inf[ g->index( &i, &j, &l ) ] = (f[ g->index( &i, &j, &l ) ] + inf[ g->index( &i, &j, &l ) ])/2.;
            }
    
    
}


#endif /* operators_h */
