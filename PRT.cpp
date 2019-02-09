//
// Created by Im YoungMin on 2019-02-08.
//

#include "PRT.h"

using namespace prt;

/**
 * Constructor.
 * @param nSamples Number of samples.
 */
PRT::PRT( size_t nSamples )
{
	_N_SAMPLES = static_cast<size_t>( pow( ceil( sqrt( nSamples ) ), 2.0 ) );		// Make sure the number of samples is a square number.
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
			_samples.emplace_back( Sample( theta, phi ) );							// Add new sample.
		}
	}
}

/**
 * Sample constructor.
 * @param t Spherical theta coordinate.
 * @param p Spherical phi coordinate.
 */
Sample::Sample( double theta, double phi )
{
	_theta = theta;			// Spherical coordinates.
	_phi = phi;

	// Position in rectangular coordinates.
	_position = { sin( _theta ) * cos( _phi ), sin( _theta ) * sin( _phi ), cos( _theta ) };
}
