//
// Created by Im YoungMin on 2019-02-08.
//

#ifndef PRT_PRT_H
#define PRT_PRT_H

#include <armadillo>
#include <random>

using namespace arma;
using namespace std;

/**
 * Precomputed Radiance Transfer namespace.
 */
namespace prt
{
	/**
	 * Store information on samples on the unit sphere, corresponding to the infinitely far environment light.
	 */
	class Sample
	{
	private:
		vec3 _position;					// Position in rectangular coordinates [x,y,z].
		double _theta, _phi;			// Spherical coordinates: \theta, \phi -- \rho = 1.
		vector<double> _sh;				// Spherical harmonics function evaluations: _N_BANDS^2.

	public:
		explicit Sample( double theta, double phi, const vector<double>& sh );
	};

	/**
	 * Implementing the Precomputed Radiance Transfer process.
	 */
	class PRT
	{
	private:
		unsigned int _N_BANDS;			// Number of Spherical-Harmonics bands.
		size_t _N_SAMPLES;				// Samples on unit sphere (must be a square number).
		vector<Sample> _samples;

		void _generateSamples();
		double _y_lm( unsigned int l, int m, double theta, double phi );
		double _K_lm( unsigned int l, int m );
		unsigned int _factorial( unsigned int x );
		unsigned int _doubleFactorial( int x );
		double _P_lm( unsigned int l, int m, double x );

	public:
		explicit PRT( size_t nSamples, unsigned int nBands = 4 );
	};
}


#endif //PRT_PRT_H
