//
//  momaxis.h
//  
//
//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
//

#ifndef momaxis_h
#define momaxis_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <complex>
#include "constant.h"

#define complex complex<double>
using namespace std;


class momaxis
{
    
public:
    
    double **k, tkx, tky;                             //kx, ky, and kz mesh
    
    double *dk;                             //grid or mesh steps
    double *weight;
    
    double *kmaxs, lattice_a0;              //Number of points
    
    int *N;
    
    long long int Ntotal;
    
    double dV;
    double Kpoint1[Ngrad];
    double Kpoint2[Ngrad];
    double DeltaKp[Ngrad];
    string imethod;
    
    momaxis( int *_N
            ,double *_kmaxs, double const *_a0 );                // Creator/Constructor Object ...
    
    ~momaxis( );                              // Destructor ...
    
    
    void basic_mem();                        // Setting memory ...
    
    void set_grid_mem( );                    // Setting memory ...
    
    void steps_sizes( int *_N
                     ,double *_kmaxs );            // Momentum calculation
    
    
    void box( int dir_index );                // Basic axis or box ...
    void set_brillouin_zone_grid( );         // Creating axis ...
    
    
    void mom_outputs( FILE *output
                     , int skip0
                     , int skip1
                     , int skip2 );          //Momentum output
    
    
    long int index( int const *i
                    ,int const *j
                    ,int const *l );
    
    void Trapz();
    void Simpson();
    void integral_method( const string NMethod );
    void checker( int *_N );

    void print_info();                       // Printing information of momentum grid
    
};




//Constructor
momaxis::momaxis( int *_N, double *_kmaxs, double const *_a0 )
{
    
    
    lattice_a0  = *_a0;
    steps_sizes( _N, _kmaxs );
    set_brillouin_zone_grid( );
    imethod = "Trapz";
    
}


void momaxis::integral_method( const string NMethod )
{
    
    imethod = NMethod;
    for (int i = 0; i<Ntotal;i++)
        weight[i]=1.;
    
    if (imethod == "Simpson")
    {
        Simpson( );
    }
    if (imethod == "Trapz")
    {
        Trapz( );
    }
    
    
}


//Destructor
momaxis::~momaxis()
{
    
    free(   dk   );
    free( kmaxs  );
    free(   N    );
    free( weight );
    
    for (int i = 0; i < Ndim; i++)
    {
        
        free( k[i] );
        
    }
    
    free( k );
}







void momaxis::checker( int *_N )
{
    if (_N[0] <= 0)
    {
        cout << "\n\n/****************OJO***************/\nNx should be larger or eq. to one, Nx >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
    
    if (_N[1] <= 0)
    {
        cout << "\n\n/***************OJO******************/\nNy should be larger or eq. to one, Ny >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
    if (_N[2] <= 0)
    {
        cout << "\n\n/*****************OJO***************/\nNz should be larger or eq. to one, Nz >= 1" << endl<< endl<< endl;
        exit (EXIT_FAILURE);
    }
    
}


//Creating basic memory
void momaxis::basic_mem()
{
    
    
    dk      = ( double* ) malloc( Ndim*sizeof( double ) );
    kmaxs   = ( double* ) malloc( Ndim*sizeof( double ) );
    N       = (    int* ) malloc( Ndim*sizeof(    int ) );
    
    
    memset( dk, 1, Ndim*sizeof( double ) );
    memset( kmaxs, 0, Ndim*sizeof( double ) );
    memset( N, 1, Ndim*sizeof( int ) );
    
    k   = ( double** ) malloc( Ndim*sizeof( double*) );
    
    
}




//Setting axis or box memory
void momaxis::set_grid_mem( )
{
    
    
    for (int i=0; i<Ndim; i++)
    {
        
        k[i] = ( double* )malloc( N[i] *sizeof( double ) );
        memset( k[i], 0, N[i]*sizeof( double ) );
        
    }
    
    weight = ( double* )malloc( Ntotal * sizeof( double ) );
    memset( weight, 1, Ntotal * sizeof( double )  );
    

    
}



//Momentum steps and sizes
void  momaxis::steps_sizes(int *_N, double *_kmax )
{
    
    /*******************************************************
     
     Conditional on kmaxs
     
     N[1]=N[2] has to be one to reduce the dimensionality
     
     of the problem to 1D,
     
     N[2]= 1 to reduce the dimentionality of the problem to 2D
     
     then, it follows:
     
     *********/
    
    checker( _N );
    
    basic_mem( );
    
    
    for ( int i=0; i<Ndim; i++ )
    {
        
        kmaxs[i]    = _kmax[i];
        N[i]        = _N[i];
        
    }
    
    
    //Building condition for symmetric grid
    for ( int i=0; i<Ndim; i++ )
    {
        
        //if ( N[i]%2==0 )
            dk[i] = 2.*kmaxs[i] / double( N[i] -1 );
        //else
          //  dk[i] = 2.*kmaxs[i] / double( N[i] );
        
    }
    
    
    //1D condition
    if( N[1]==1 && N[2]==1 )
    {
        
        dk[1] = 1.;
        dk[2] = 1.;
        
    }
    
    
    //2D condition
    if( N[2]==1 )
        dk[2] = 1.;
    
    
    Ntotal  = N[0]*N[1]*N[2];
    dV      = dk[0]*dk[1]*dk[2];
    
    set_grid_mem( );
    
    
}




//Box or grid  along dir-index
void momaxis::box( int dir_index )
{
    
    double kmin;
    
    if ( N[dir_index] > 1 )
    {
        
    
        kmin = kmaxs[dir_index]*(-1.) ; //Shifting whole BZ
        //if ( dir_index==1 )
        //    kmin = dk[dir_index]/2.;//kmaxs[dir_index]*(-1.) ;
        
        for ( int i = 0; i < N[dir_index]; i++ )
            k[dir_index][i] = kmin + dk[dir_index]*i;
        
        
    }
    else
    {
        
        
        k[dir_index][0] = 0.;
        
        
    }
    
    //k[1][0] = -0.1   ;
    //k[1][0]=0.2;
}


//index calculation
long int momaxis::index( int const *i, int const *j, int const *l )
{
    
    return N[1]*N[0]*(*l) + N[0]*(*j) + (*i);
    
}


void momaxis::Trapz()
{
    
    int ia0 =0;
    int ie0 =N[0]-1;
    
    int ja0=0;
    int je0=N[1]-1;
    
    int le=0;
    
    weight[ index( &ia0, &ja0, &le) ] = 1./4.;
    weight[ index( &ie0, &ja0, &le) ] = 1./4.;

    weight[ index( &ia0, &je0, &le) ] = 1./4.;
    weight[ index( &ie0, &je0, &le) ] = 1./4.;
    
    
    for ( int i=1; i<ie0; i++)
    {
        weight[ index( &i, &ja0, &le) ] = 2./4.;
        weight[ index( &i, &je0, &le) ] = 2./4.;
    }
    
    
    for ( int j=1; j<je0; j++ )
    {
        weight[ index( &ia0, &j, &le) ] = 2./4.;
        weight[ index( &ie0, &j, &le) ] = 2./4.;
    }
    
    
    for (int j=1; j<je0; j++)
        for (int i=1; i<ie0; i++)
            weight[ index( &i, &j, &le) ] = 4./4.;

}


void momaxis::Simpson()
{
    
    int i_aux, j_aux;
    int ia0 = 0;
    int ie0 = N[0]-1;
    
    int ja0 = 0;
    int je0 = N[1]-1;
    
    int le=0;
    
    weight[ index( &ia0, &ja0, &le) ] = 1./9.;
    weight[ index( &ie0, &ja0, &le) ] = 1./9.;
    
    weight[ index( &ia0, &je0, &le) ] = 1./9.;
    weight[ index( &ie0, &je0, &le) ] = 1./9.;
    
    
    for ( int i = 1; i < int( (N[0]+1)/2.); i++ )
    {
        

        i_aux=2*i-1;
        
        weight[ index( &i_aux, &ja0, &le) ] = 4./9.;
        weight[ index( &i_aux, &je0, &le) ] = 4./9.;
        
    }
    
    
    for ( int i=1; i < int( N[0]/2.); i++ )
    {
        i_aux=2*i;
        weight[ index( &i_aux, &ja0, &le) ] = 2./9.;
        weight[ index( &i_aux, &je0, &le) ] = 2./9.;
    }
    
    
    for ( int j = 1; j < int( (N[1]+1) /2. ); j++ )
    {
        
        j_aux = 2*j-1;
        weight[ index( &ia0, &j_aux, &le) ] = 4./9.;
        weight[ index( &ie0, &j_aux, &le) ] = 4./9.;
        
    }
    
    for ( int j = 1; j < int( N[1]/2. ); j++ )
    {
        
        j_aux = 2*j;
        weight[ index( &ia0, &j_aux, &le) ] = 2./9.;
        weight[ index( &ie0, &j_aux, &le) ] = 2./9.;
        
    }

    
    
    
    for ( int j=1; j < int( (N[1]+1)/2. ); j++ )
        for ( int i=1; i < int( (N[0]+1)/2. ); i++ )
        {
            
            i_aux=2*i-1;
            j_aux=2*j-1;
            
            weight[ index( &i_aux, &j_aux, &le) ] = 16./9.;
            
        }
    
    for ( int j=1; j < int( N[1]/2. ); j++ )
        for ( int i=1; i < int( (N[0]+1)/2. ); i++ )
        {
            j_aux=2*j;
            i_aux=2*i-1;
            weight[ index( &i_aux, &j_aux, &le) ] = 8./9.;
        }
    
    for ( int j=1; j < int( (N[1]+1)/2. ); j++ )
        for ( int i=1; i < int( N[0]/2. ); i++ )
        {

            i_aux=2*i;
            j_aux=2*j-1;
            weight[ index( &i_aux, &j_aux, &le) ] = 8./9.;
        }
    

    
    
    for ( int j=1; j < int( N[1]/2. ); j++ )
        for ( int i=1; i < int( N[0]/2. ); i++ )
        {
            i_aux=2*i;
            j_aux=2*j;
            weight[ index( &i_aux, &j_aux, &le) ] = 4./9.;
        }
    
    
}




//First Brillouine Zone
void momaxis::set_brillouin_zone_grid()
{

    
    //x-direction
    box( 0 );
    
    
    //y-direction
    box( 1 );
    
    
    //z-direction
    box( 2 );
    
    
    
    Kpoint1[0] = -4.*pi/lattice_a0/sqrt(3.)/3.;
    Kpoint1[1] = 0.;

    Kpoint2[0] = 2.*pi/lattice_a0/sqrt(3.)/3.;
    Kpoint2[1] = 2.*pi/lattice_a0/3.;
    
    DeltaKp[0] = Kpoint2[0] - Kpoint1[0];
    DeltaKp[1] = Kpoint2[1] - Kpoint1[1];
    
}



//Outputs of grid
void momaxis::mom_outputs( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
{
    
    
    fprintf( output,"%d \n", skip0  );
    fprintf( output,"%d \n", N[0]  );
    fprintf( output,"%e \n\n", dk[0]  );
    
    
    
    for ( int i=0; i<N[0]/skip0; i++ )
        fprintf( output, "%e \n", k[0][ i*skip0 ] );
    
    
    fprintf( output,"\n%d \n", skip1  );
    fprintf( output,"%d \n", N[1]  );
    fprintf( output,"%e \n", dk[1]  );
    fprintf(output,"\n");
    
    
    for ( int i=0; i<N[1]/skip1; i++ )
        fprintf( output, "%e \n", k[1][ i*skip1 ] );
    
    
    fprintf( output,"\n%d \n", skip2  );
    fprintf( output,"%d \n", N[2]  );
    fprintf( output,"%e \n", dk[2]  );
    fprintf(output,"\n");
    
    
    for ( int i=0; i<N[2]/skip2; i++ )
        fprintf( output, "%e \n", k[2][ i*skip2 ] );
    
    
    fflush( output );
    
    
}

void momaxis::print_info()
{
    cout << "\n==========================================\n";
    cout << "First  K' point = ("<< Kpoint1[0] << " , " << Kpoint1[1] << " ) a.u.\n";
    cout << "Second K' point = ("<< Kpoint2[0] << " , " << Kpoint2[1] << " ) a.u.\n";
    cout << "Delta-K' point = (" << DeltaKp[0] << " , " << DeltaKp[1] << " ) a.u.\n========================\n\n";
}

#endif /* momaxis_h */
