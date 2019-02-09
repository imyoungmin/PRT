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
	public:
		Sample( double theta, double phi );
	};

	/**
	 * Implementing the Precomputed Radiance Transfer process.
	 */
	class PRT
	{
	private:
		size_t _N_SAMPLES;				// Samples on unit sphere (must be a square number).
		vector<Sample> _samples;

		void _generateSamples();

	public:
		explicit PRT( size_t nSamples );
	};
}


#endif //PRT_PRT_H
