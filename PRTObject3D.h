//
// Created by Im YoungMin on 2019-02-14.
//

#ifndef PRT_PRTOBJECT3D_H
#define PRT_PRTOBJECT3D_H

#include <OpenGL/gl3.h>
#include <armadillo>
#include "Transformations.h"
#include "PRTConstants.h"

using namespace std;
using namespace arma;

namespace prt
{
	class Object3D;								// Forward declaration.

	/**
	 * Organizing vertices and normals into triangles by keeping references to original vectors.
	 */
	class Triangle
	{
	private:
		const Object3D* _objectPtr;				// Object this triangle belongs to.
		const unsigned int _v0;					// Vertex indices.
		const unsigned int _v1;
		const unsigned int _v2;
		vec3 _normal;							// Triangle's supporting plane normal vector, assuming vertices are given in CCW direction.
		double _d;								// Given the equation of the supporting plane: n \dot x = d

	public:
		explicit Triangle( const Object3D* oPtr, unsigned int ver0, unsigned int ver1, unsigned int ver2 );
		bool rayIntersection( const vec3& p, const vec3& d ) const;
	};


	/**
	 * Especialized 3D geometry/object for PRT.
	 */
	class Object3D
	{
	private:
		string _name;							// Object name.
		vector<Triangle*> _triangles;			// Geometry data: vertices, normals, and triangles with references.
		vector<vec3> _vertices;
		vector<Triangle*> _verticesTriangles;	// Each vertex needs to know which triangle it belongs to to avoid self intersections.
		vector<vec3> _normals;
		vector<vec3*> _shCoefficients;			// Spherical harmonics projection coefficients for each color channel.
		vec3 _color;							// Object uniform color: currently supported only RGB (no transparency).
		GLsizei _verticesCount;					// Number of vertices.
		GLuint _bufferID;						// Buffer ID given by OpenGL.
		GLuint _tboID;							// We use TBOs to store the spherical harmonics coefficients for each vertex.
		GLuint _tboTextureID;					// Associated texture for accessing TBO.

		GLsizei _getData( vector<float>& outVs, vector<float>& outNs ) const;

	public:
		explicit Object3D( const char* name, const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color );
		GLuint getBufferID() const;
		GLuint getTBOID() const;
		GLuint getTBOTextureID() const;
		GLsizei getVerticesCount() const;
		void resetSHCoefficients();
		void scaleSHCoefficients( double s );
		void accumulateSHCoefficients( unsigned int vIndex, unsigned int shIndex, const vec3& value );
		void loadSHCoefficientsIntoTexture();
		const vec3& getVertexPositionAt( unsigned int index ) const;
		const Triangle* getVertexTrianglePtrAt( unsigned int index ) const;
		const vec3& getVertexNormalAt( unsigned int index ) const;
		const vec3& getColor() const;
		const string& getName() const;
		bool rayIntersection( const vec3& p, const vec3& d, const Triangle* trianglePtr );
		void deallocateGeometries();
		~Object3D();
	};
}


#endif //PRT_PRTOBJECT3D_H
