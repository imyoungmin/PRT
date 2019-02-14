//
// Created by Im YoungMin on 2019-02-14.
//

#include "PRTObject3D.h"

using namespace prt;

/**
 * Constructor.
 * @param vertices Object vertices.
 * @param normals Object normal vectors.
 * @param T Transformation matrix that takes object coordinates into world coordinates.
 * @param color Object RGB color.
 * @param nrCoefficients Number of spherical harmonics projection coefficients.
 * @return
 */
Object3D::Object3D( const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color, unsigned int nrCoefficients )
{
	// Transform vertices and normals to world coordinates and register triangles with references.
	mat33 R = Tx::getInvTransModelView( T, false );		// Don't assume uniform scaling: to transform the normals.

	int i = 0;
	while( i < vertices.size() )
	{
		for( int j = 0; j < 3; j++, i++ )				// Three vertices per triangle.
		{
			vec4 hTV = T * vec4( { vertices[i][0], vertices[i][1], vertices[i][2], 1.0 } );
			_vertices.emplace_back( hTV.head( 3 ) );
			_normals.emplace_back( R * normals[i] );

			// Allocate space for spherical harmonics coefficients for each vertex (for the three channels RGB).
			_shCoefficients.push_back( new vec3[nrCoefficients] );
		}

		// Register triangle.
		_triangles.emplace_back( _vertices[i-3], _vertices[i-2], _vertices[i-1], _normals[i-3], _normals[i-2], _normals[i-1] );
	}

	// Create OpenGL objects and variables.
	glGenBuffers( 1, &(_bufferID) );
	glBindBuffer( GL_ARRAY_BUFFER, _bufferID );
	vector<float> vertexPositions;
	vector<float> normalComponents;
	_verticesCount = _getData( vertexPositions, normalComponents );

	// Allocate space for vertices, normals, and spherical harmonics projection coefficients.
	const size_t size3D = sizeof(float) * vertexPositions.size();							// Size of 3D arrays in bytes.
	const size_t sizeSH = sizeof(float) * _verticesCount * nrCoefficients * 3;				// Make space for 3 channels per coefficient: RGB-RGB-RGB-...
	glBufferData( GL_ARRAY_BUFFER, 2 * size3D + sizeSH, nullptr, GL_DYNAMIC_DRAW );			// Dynamic as we'll copy sh coefficients later.
	glBufferSubData( GL_ARRAY_BUFFER, 0, size3D, vertexPositions.data() );					// Copy positions.
	glBufferSubData( GL_ARRAY_BUFFER, size3D, size3D, normalComponents.data() );			// Copy normals.
	// TODO: Copy sh coefficients into buffer.
}

/**
 * Retrieve the OpenGL buffer ID.
 * @return VBO.
 */
GLuint Object3D::getBufferID() const
{
	return _bufferID;
}

/**
 * Retrieve number of vertices.
 * @return Vertex count.
 */
GLsizei Object3D::getVerticesCount() const
{
	return _verticesCount;
}

/**
 * Collect the vertex, uv, and normal coordinates into linear vectors of scalars.
 * @param outVs Output vertices data.
 * @param outNs Output normals data.
 * @return Number of processed vertices (= normals).
 */
GLsizei Object3D::_getData( vector<float>& outVs, vector<float>& outNs ) const
{
	auto N = static_cast<GLsizei>( _vertices.size() );

	for( int i = 0; i < N; i++ )
	{
		// Vertices.
		outVs.push_back( _vertices[i][0] );			// X-coordinate.
		outVs.push_back( _vertices[i][1] );			// Y-coordinate.
		outVs.push_back( _vertices[i][2] );			// Z-coordinate.

		// Normals.
		outNs.push_back( _normals[i][0] );			// X-coordinate.
		outNs.push_back( _normals[i][1] );			// Y-coordinate.
		outNs.push_back( _normals[i][2] );			// Z-coordinate.
	}

	return N;
}

/**
 * Destructor.
 */
Object3D::~Object3D()
{
	// Deallocate arrays of coefficients for RGB for each vertex.
	for( const vec3* coefficients : _shCoefficients )
		delete [] coefficients;
}


//////////////////////////////////////////////////// Triangle class ////////////////////////////////////////////////////

/**
 * Constructor.
 * @param ver0 First vertex (in CCW order).
 * @param ver1 Second vertex.
 * @param ver2 Third vertex.
 * @param nor0 First normal.
 * @param nor1 Second normal.
 * @param nor2 Third normal.
 */
Triangle::Triangle( const vec3& ver0, const vec3& ver1, const vec3& ver2, const vec3& nor0, const vec3& nor1, const vec3& nor2 ):
	_v0Ref( ver0 ), _v1Ref( ver1 ), _v2Ref( ver2 ), _n0Ref( nor0 ), _n1Ref( nor1 ), _n2Ref( nor2 )
{}

