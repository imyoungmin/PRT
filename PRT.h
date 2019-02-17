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
#include "PRTConstants.h"
#include "PRTObject3D.h"

using namespace arma;
using namespace std;
using namespace std::chrono;

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
		vec3 _position;						// Position in rectangular coordinates [x,y,z].
		double _theta, _phi;				// Spherical coordinates: \theta, \phi -- \rho = 1.
		vector<double> _sh;					// Spherical harmonics function evaluations: _N_BANDS^2.

	public:
		explicit Sample( double theta, double phi, const vector<double>& sh );
		double getSHValueAt( int index ) const;
		const vec3& getPosition() const;
	};

	/**
	 * Implementing the Precomputed Radiance Transfer process.
	 */
	class PRT
	{
	private:
		size_t _N_SAMPLES;					// Samples on unit sphere (must be a square number).
		vector<Sample> _samples;

		unsigned char* _cubeMapFaces[6];	// Cube map faces.
		int _cubeMapFaceWidth = -1;
		int _cubeMapFaceNrChannels = -1;
		vec3 _cubeMapFacesNormals[6] = {
				{  1,  0,  0 },				// Right.
				{ -1,  0,  0 },				// Left.
				{  0,  1,  0 },				// Top.
				{  0, -1,  0 },				// Bottom.
				{  0,  0,  1 },				// Front.
				{  0,  0, -1 }				// Back.
		};

		vec3* _lightCoefficients;			// Projection of sampled lighting into the _N_BANDS^2 spherical harmonics functions.  Notice: RGB channels.
		float* _lightCoefficientsArray;		// Flat version of light coefficients, given in RGBRGBRGB...RGB

		vector<unique_ptr<Object3D>> _objects;			// List of PRT 3D objects to be shade.

		GLuint _renderingProgram;			// OpenGL rendering program.

		void _generateSamples();
		void _generateCubeMap( const vector<string>& facesFileNames );
		void _getPixel( unsigned int x, unsigned int y, unsigned int face, unsigned char* output ) const;
		void _queryCubeMap( const vec3& query, unsigned char* output ) const;
		void _projectLighting();
		void _unshadowedDiffuseTransferProjection();
		void _shadowedDiffuseTransferProjection();
		int _visibility( const vec3& p, const unsigned int s, const Triangle* trianglePtr );
		double _y_lm( int l, int m, double theta, double phi );
		double _K_lm( int l, int m );
		int _factorial( int x );
		int _doubleFactorial( int x );
		double _P_lm( int l, int m, double x );

	public:
		PRT();
		void init( size_t nSamples, const vector<string>& facesFileNames, GLuint program );
		~PRT();
		const unsigned char* getCubeMapFaceData( int faceIndex ) const;
		int getCubeMapFaceWidth() const;
		int getCubeMapFaceNrChannels() const;
		void addObject( const char* name, const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color );
		void precomputeRadianceTransfer();
		void renderObjects( const mat44& Projection, const mat44& Camera, const mat44& Model );
	};
}


#endif //PRT_PRT_H
