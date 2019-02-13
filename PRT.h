//
// Created by Im YoungMin on 2019-02-08.
//

#ifndef PRT_PRT_H
#define PRT_PRT_H

#include <armadillo>
#include <random>
#include "stb_image.h"
#include "Configuration.h"
#include "Transformations.h"

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
		unsigned int _N_BANDS;				// Number of Spherical-Harmonics bands.
		size_t _N_SAMPLES;					// Samples on unit sphere (must be a square number).
		vector<Sample> _samples;

		unsigned char* _cubeMapFaces[6];	// Cube map faces.
		int _cubeMapFaceWidth = -1;
		int _cubeMapFaceNrChannels = -1;
		vec3 _cubeMapFacesNormals[6] = {
				 Tx::X_AXIS,				// Right.
				-Tx::X_AXIS,				// Left.
				 Tx::Y_AXIS,				// Top.
				-Tx::Y_AXIS,				// Bottom.
				 Tx::Z_AXIS,				// Front.
				-Tx::Z_AXIS					// Back.
		};

		void _generateSamples();
		void _generateCubeMap( const vector<string>& facesFileNames );
		void _getPixel( unsigned int x, unsigned int y, unsigned int face, unsigned char* output ) const;
		double _y_lm( int l, int m, double theta, double phi );
		double _K_lm( int l, int m );
		int _factorial( int x );
		int _doubleFactorial( int x );
		double _P_lm( int l, int m, double x );

	public:
		PRT();
		void init( size_t nSamples, const vector<string>& facesFileNames, unsigned int nBands = 4 );
		~PRT();
		const unsigned char* getCubeMapFaceData( int faceIndex ) const;
		int getCubeMapFaceWidth() const;
		int getCubeMapFaceNrChannels() const;
		void queryCubeMap( const vec3& query, unsigned char* output ) const;
	};
}


#endif //PRT_PRT_H
