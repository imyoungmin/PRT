//
// Created by Im YoungMin on 2019-02-14.
//

#ifndef PRT_PRTOBJECT3D_H
#define PRT_PRTOBJECT3D_H

#include <OpenGL/gl3.h>
#include <armadillo>
#include "Transformations.h"

using namespace std;
using namespace arma;

namespace prt
{
	/**
	 * Organizing vertices and normals into triangles by keeping references to original vectors.
	 */
	class Triangle
	{
	private:
		const vec3& _v0Ref, _v1Ref, _v2Ref;		// Vertices.
		const vec3& _n0Ref, _n1Ref, _n2Ref;		// Normal vectors.

	public:
		explicit Triangle( const vec3& ver0, const vec3& ver1, const vec3& ver2, const vec3& nor0, const vec3& nor1, const vec3& nor2 );
	};


	/**
	 * Especialized 3D geometry/object for PRT.
	 */
	class Object3D
	{
	private:
		vector<Triangle> _triangles;			// Geometry data: vertices, normals, and triangles with references.
		vector<vec3> _vertices;
		vector<vec3> _normals;
		vector<vec3*> _shCoefficients;			// Spherical harmonics projection coefficients for each color channel.
		vec3 _color;							// Object uniform color: currently supported only RGB (no transparency).
		unsigned int _nrCoefficients;			// Number of spherical harmonics coefficients per vertex per channel.
		GLsizei _verticesCount;					// Number of vertices.
		GLuint _bufferID;						// Buffer ID given by OpenGL.

		GLsizei _getData( vector<float>& outVs, vector<float>& outNs ) const;

	public:
		explicit Object3D( const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color, unsigned int nrCoefficients );
		GLuint getBufferID() const;
		GLsizei getVerticesCount() const;
		void resetSHCoefficients();
		void scaleSHCoefficients( double s );
		void accumulateSHCoefficients( unsigned int vIndex, unsigned int shIndex, const vec3& value );
		const vec3& getVertexPositionAt( unsigned int index ) const;
		const vec3& getVertexNormalAt( unsigned int index ) const;
		const vec3& getColor() const;
		~Object3D();
	};
}


#endif //PRT_PRTOBJECT3D_H
