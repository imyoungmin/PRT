//
// Created by Im YoungMin on 2019-02-08.
//

#include "PRT.h"

using namespace prt;

/**
 * Constructor.
 * @param nSamples Number of samples.
 * @param nBands Number of spherical-harmonics bands.
 */
PRT::PRT( size_t nSamples, unsigned int nBands )
{
	_N_SAMPLES = static_cast<size_t>( pow( ceil( sqrt( nSamples ) ), 2.0 ) );		// Make sure the number of samples is a square number.
	_N_BANDS = nBands;
	_generateSamples();
}

/**
 * Generate samples uniformly distributed on the unit sphere.
 */
void PRT::_generateSamples()
{
	std::random_device rd;					// Request random data from OS.
	std::mt19937 generator( rd() );
	uniform_real_distribution<double> uniform( 0, 1 );
	const auto N = static_cast<size_t>( sqrt( _N_SAMPLES ) );
	for( int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			double x = ( i + uniform( generator ) ) / static_cast<double>( N );		// Uniformly distributed coordinates in the unit square
			double y = ( j + uniform( generator ) ) / static_cast<double>( N );		// with the (0,0) located at the lower left corner.

			double theta = 2.0 * acos( sqrt( 1.0 - x ) );							// Spherical coordinates: \theta, \phi \in [0, 2\pi)
			double phi = 2.0 * M_PI * y;

			// Generate the _N_BANDS^2 spherical-harmonics functions for this sample.
			vector<double> sh;
			for( unsigned l = 0; l < _N_BANDS; l++ )
				for( int m = -l; m <= l; m++ )
					sh.emplace_back( _y_lm( l, m, theta, phi ) );					// The ith function evaluation is at i = l(l+1) + m.

			_samples.emplace_back( Sample( theta, phi, sh ) );						// Add new sample.
		}
	}
}

/**
 * Spherical harmonics function.
 * @param l Band index.
 * @param m Offset index within a band.
 * @param theta Unit sphere theta angle (longitude).
 * @param phi Unit sphere phi angle (latitude).
 * @return Function value.
 */
double PRT::_y_lm( unsigned int l, int m, double theta, double phi )
{
	// We have three cases: m > 0, m < 0, and m = 0.
	if( m > 0 )
		return M_SQRT2 * _K_lm( l, m ) * cos( m * phi ) * _P_lm( l, m, cos( theta ) );

	if( m < 0 )
		return M_SQRT2 * _K_lm( l, m ) * sin( -m * phi ) * _P_lm( l, -m, cos( theta ) );

	return _K_lm( l, 0 ) * _P_lm( l, 0, cos( theta ) );
}

/**
 * Associated Legendre polynomial.
 * https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 * @param l Band index.
 * @param m Offset index within a band.
 * @param x Input value to polynomial.
 * @return Associated Legendre polynomial evaluated at x.
 */
double PRT::_P_lm( unsigned int l, int m, double x )
{
	if( l == m )
		return pow( -1, l ) * _doubleFactorial( 2 * l - 1 ) * pow( 1 - x * x, static_cast<double>( l )/2.0 );

	if( l == m + 1 )
		return x * ( 2 * m + 1 ) * _P_lm( l - 1, m, x );

	return ( x * ( 2 * l - 1 ) * _P_lm( l - 1, m, x ) - ( l + m - 1) * _P_lm( l - 2, m, x ) ) / ( l - m );
}

/**
 * Spherical harmonics normalizing constant.
 * @param l Band index.
 * @param m Offset index within a band.
 * @return Normalizing constant.
 */
double PRT::_K_lm( unsigned int l, int m )
{
	if( m == 0 )
		return sqrt( ( 2.0 * l + 1.0 ) / ( 4.0 * M_PI ) );

	return sqrt( ( 2.0 * l + 1.0 ) / ( 4.0 * M_PI ) *
				static_cast<double>( _factorial( l - abs( m ) ) ) / static_cast<double>( _factorial( l + abs( m ) ) ) );
}

/**
 * Compute the factorial of a non-negative integer.
 * @param x Input value.
 * @return x!
 */
unsigned int PRT::_factorial( unsigned int x )
{
	if( x <= 1)
		return 1;

	unsigned int factorial = 1;
	for( unsigned int i = 2; i <= x; i++ )
		factorial *= i;
	return factorial;
}

/**
 * Double factorial of an integer.
 * https://en.wikipedia.org/wiki/Double_factorial
 * @param x Input value.
 * @return x!!
 */
unsigned int PRT::_doubleFactorial( int x )
{
	if( x <= 1 )
		return 1;
	return x * _doubleFactorial( x - 2);
}

/**
 * Sample constructor.
 * @param t Spherical theta coordinate.
 * @param p Spherical phi coordinate.
 * @param sh Spherical harmonics function evaluations.
 */
Sample::Sample( double theta, double phi, const vector<double>& sh )
{
	_theta = theta;			// Spherical coordinates.
	_phi = phi;
	_sh = vector<double>( sh );

	// Position in rectangular coordinates.
	_position = { sin( _theta ) * cos( _phi ), sin( _theta ) * sin( _phi ), cos( _theta ) };
}
