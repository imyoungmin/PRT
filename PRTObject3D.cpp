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
	_nrCoefficients = nrCoefficients;
	_color = color;

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
			_shCoefficients.push_back( new vec3[_nrCoefficients] );
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

	// Allocate space for vertices.
	const size_t size3D = sizeof(float) * vertexPositions.size();							// Size of 3D arrays in bytes.
	glBufferData( GL_ARRAY_BUFFER, size3D, vertexPositions.data(), GL_STATIC_DRAW );		// Don't need to send the normals.
}

/**
 * Copy spherical harmonics coefficients into a texture.
 * This function should be called whenever we have projected the transfer function for this object.
 * We use a texture buffer object because we can't send an array of floats for each vertex to the vertex shader.
 * How to: https://gist.github.com/roxlu/5090067.
 */
void Object3D::loadSHCoefficientsIntoTexture()
{
	glGenBuffers( 1, &_tboID );								// Create a texture buffer object
	glBindBuffer( GL_TEXTURE_BUFFER, _tboID );
	glGenTextures( 1, &_tboTextureID );

	vector<float> coeffs;
	for( const vec3* coefficients : _shCoefficients )		// From vertex 0 -> (total - 1).
	{														// Create a row per vertex in texture.
		for( int i = 0; i < _nrCoefficients; i++ )
		{
			coeffs.push_back( 0.1/*coefficients[i][0]*/ );			// Load three channels: R.
			coeffs.push_back( 0.05/*coefficients[i][1]*/ );			// G.
			coeffs.push_back( 0/*coefficients[i][2]*/ );			// B.
		}
	}

	// Copy data into texture buffer object: it's now available in texture for rendering.
	glBufferData( GL_TEXTURE_BUFFER, sizeof(float) * coeffs.size(), coeffs.data(), GL_STATIC_DRAW );
	glBindBuffer( GL_TEXTURE_BUFFER, 0 );

	// We can now free the coefficients' memory from the application.
	_deallocateGeometries();
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
 * Retrieve the OpenGL TBO for spherical harmonics coefficients.
 * @return TBO ID.
 */
GLuint Object3D::getTBOID() const
{
	return _tboID;
}

/**
 * Retrieve the OpenGL texture associated with the TBO for spherical harmonics coefficients.
 * @return Texture TBO ID.
 */
GLuint Object3D::getTBOTextureID() const
{
	return _tboTextureID;
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
 * Reset spherical harmonics coefficients for each vertex and for each channel.
 */
void Object3D::resetSHCoefficients()
{
	for( int i = 0; i < _verticesCount; i++ )	// For each vertex in the object, initialize the corresponding coefficients.
	{
		for( int j = 0; j < _nrCoefficients; j++ )
			_shCoefficients[i][j] = { 0, 0, 0 };
	}
}

/**
 * Scale all of the spherical harmonics cofficients in this object by a scalar.
 * @param s Scale factor.
 */
void Object3D::scaleSHCoefficients( double s )
{
	for( int i = 0; i < _verticesCount; i++ )
	{
		for( int j = 0; j < _nrCoefficients; j++ )
			_shCoefficients[i][j] *= s;
	}
}

/**
 * Accumulate spherical harmonics coefficients at a given vertex and spherical harmonics index.
 * @param vIndex Vertex index.
 * @param shIndex Spherical harmonics index.
 * @param value Operand to add.
 */
void Object3D::accumulateSHCoefficients( unsigned int vIndex, unsigned int shIndex, const vec3& value )
{
	_shCoefficients[vIndex][shIndex] += value;
}

/**
 * Retrieve a vertex position.
 * @param index Vertex index.
 * @return 3D position.
 */
const vec3& Object3D::getVertexPositionAt( unsigned int index ) const
{
	return _vertices[ index ];
}

/**
 * Retrieve a vertex normal.
 * @param index Normal index.
 * @return 3D normal vector.
 */
const vec3& Object3D::getVertexNormalAt( unsigned int index ) const
{
	return _normals[ index ];
}

/**
 * Retrieve object's color.
 * @return RGB vector.
 */
const vec3& Object3D::getColor() const
{
	return _color;
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
 * Free memory for unnecessary geometries.
 */
void Object3D::_deallocateGeometries()
{
	for( int i = 0; i < _verticesCount; i++ )
	{
		if( _shCoefficients[i] )
		{
			delete [] _shCoefficients[i];			// Deallocate arrays of coefficients for RGB for each vertex.
			_shCoefficients[i] = nullptr;
		}
	}
}

/**
 * Destructor.
 */
Object3D::~Object3D()
{
	_deallocateGeometries();
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

