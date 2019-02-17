//
// Created by Im YoungMin on 2019-02-14.
//

#include "PRTObject3D.h"

using namespace prt;

/**
 * Constructor.
 * @param name Object name.
 * @param vertices Object vertices.
 * @param normals Object normal vectors.
 * @param T Transformation matrix that takes object coordinates into world coordinates.
 * @param color Object RGB color.
 * @return
 */
Object3D::Object3D( const char* name, const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color )
{
	_name = string( name );
	_color = color;

	// Transform vertices and normals to world coordinates and register triangles with references.
	mat33 R = Tx::getInvTransModelView( T, false );		// Don't assume uniform scaling: to transform the normals.

	unsigned int i = 0;
	while( i < vertices.size() )
	{
		for( int j = 0; j < 3; j++, i++ )				// Three vertices per triangle.
		{
			vec4 hTV = T * vec4( { vertices[i][0], vertices[i][1], vertices[i][2], 1.0 } );
			_vertices.emplace_back( hTV.head( 3 ) );
			_normals.emplace_back( normalise( R * normals[i] ) );

			// Allocate space for spherical harmonics coefficients for each vertex (for the three channels RGB).
			_shCoefficients.push_back( new vec3[N_BANDS * N_BANDS] );
		}

		// Register triangle and reference this new triangle to each vertex (circular link).
		auto triangle = new Triangle( this, i-3, i-2, i-1 );
		_triangles.emplace_back( triangle );
		_verticesTriangles.push_back( triangle );		// Past three vertices belong to same triangle.
		_verticesTriangles.push_back( triangle );
		_verticesTriangles.push_back( triangle );
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
		for( int i = 0; i < N_BANDS * N_BANDS; i++ )
		{
			coeffs.push_back( coefficients[i][0] );			// Load three channels: R.
			coeffs.push_back( coefficients[i][1] );			// G.
			coeffs.push_back( coefficients[i][2] );			// B.
		}
	}

	// Verify maximum size of texture buffer object.
	int limit;
	glGetIntegerv( GL_MAX_TEXTURE_BUFFER_SIZE, &limit );
	size_t coefficientsSize = sizeof(float) * coeffs.size();
	if( coefficientsSize > limit )
	{
		cerr << "[PRTObject3D][" << _name << "] Cannot write to TBO to store spherical harmonics coefficients.  Max size exceeded: " << coefficientsSize << endl;
		exit( EXIT_FAILURE );
	}

	// Copy data into texture buffer object: it's now available in texture for rendering.
	glBufferData( GL_TEXTURE_BUFFER, coefficientsSize, coeffs.data(), GL_STATIC_DRAW );
	glBindBuffer( GL_TEXTURE_BUFFER, 0 );
	cout << "[PRTObject3D][" << _name << "] Successfully wrote " << coeffs.size() << " spherical harmonics coefficients to TBO!" << endl;
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
		for( int j = 0; j < N_BANDS * N_BANDS; j++ )
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
		for( int j = 0; j < N_BANDS * N_BANDS; j++ )
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
 * @return Reference to 3D position.
 */
const vec3& Object3D::getVertexPositionAt( unsigned int index ) const
{
	return _vertices[ index ];
}

/**
 * Retrieve a vertex normal.
 * @param index Normal index.
 * @return Reference to 3D normal vector.
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
 * Retrieve object's name.
 * @return Object's name.
 */
const string& Object3D::getName() const
{
	return _name;
}

/**
 * Retrieve a vertex pointer to the triangle it belongs to.
 * @param index Vertex triangle index.
 * @return Pointer to vertex's triangle.
 */
const Triangle* Object3D::getVertexTrianglePtrAt( unsigned int index ) const
{
	return _verticesTriangles[ index ];
}

/**
 * Check whether a point (i.e. a vertex on some object) intersects this object.
 * @param p Ray origin (i.e. a vertex position).
 * @param d Ray direction (i.e. in the direction of a lighting sample).
 * @param trianglePtr Pointer to triangle that p belongs to.
 * @return True if ray intersects this object, false otherwise.
 */
bool Object3D::rayIntersection( const vec3& p, const vec3& d, const Triangle* trianglePtr )
{
	for( const Triangle* triangle : _triangles )
	{
		if( triangle == trianglePtr )				// Skip checking triangle that the input vertex position belongs to.
			continue;

		if( triangle->rayIntersection( p, d ) )
			return true;
	}

	return false;
}

/**
 * Collect the vertex, uv, and normal coordinates into linear vectors of scalars.
 * @param outVs Output vertices data.
 * @param outNs Output normals data.
 * @return Number of processed vertices (= normals).
 */
GLsizei Object3D::_getData( vector<float>& outVs, vector<float>& outNs ) const
{
	// Clear contents of vector containers.
	outVs.clear();
	outNs.clear();

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
void Object3D::deallocateGeometries()
{
	for( int i = 0; i < _verticesCount; i++ )
	{
		if( _shCoefficients[i] )
		{
			delete [] _shCoefficients[i];			// Deallocate arrays of coefficients for RGB for each vertex.
			_shCoefficients[i] = nullptr;
		}
	}

	for( int i = 0; i < _triangles.size(); i++ )
	{
		if( _triangles[i] )
		{
			delete _triangles[i];
			_triangles[i] = nullptr;
		}
	}
}

/**
 * Destructor.
 */
Object3D::~Object3D()
{
	deallocateGeometries();
	cout << "[PRT][" << _name << "] destroyed!" << endl;
}


//////////////////////////////////////////////////// Triangle class ////////////////////////////////////////////////////

/**
 * Constructor.
 * @param oPtr Pointer to object that owns the triangle.
 * @param ver0 First vertex (in CCW order).
 * @param ver1 Second vertex.
 * @param ver2 Third vertex.
 */
Triangle::Triangle( const Object3D* oPtr, unsigned int ver0, unsigned int ver1, unsigned int ver2 ):
	_objectPtr( oPtr), _v0( ver0 ), _v1( ver1 ), _v2( ver2 )
{
	// Get the vertices from bound object.
	const vec3& v0 = _objectPtr->getVertexPositionAt( ver0 );
	const vec3& v1 = _objectPtr->getVertexPositionAt( ver1 );
	const vec3& v2 = _objectPtr->getVertexPositionAt( ver2 );

	_normal = cross( v1 - v0, v2 - v0 );	// Calculate supporting plane normal's vector.
	_normal = normalise( _normal );
	_d = dot( _normal, v0 );				// Parameter d in plane equation n \dot x = d.
}

/**
 * Check if a ray r(t) = p + td intersects this triangle.
 * @param p Ray origin.
 * @param d Ray direction.
 * @return True if ray intersects triangle, false otherwise.
 */
bool Triangle::rayIntersection( const vec3& p, const vec3& d ) const
{
	double nDotD = dot( _normal, d );
	double eps = 0.00001;									// Some numerical tolerance.
	if( -eps < nDotD && nDotD < eps )						// Ray is (almost) parallel to plane? Dot product == 0.
		return false;

	double t = ( _d - dot( _normal, p ) ) / nDotD;			// Parameter in ray equation that yields intersectino with triangle's plane.
	if( t < eps )											// Check if p already on the plane (t == 0), or triangle is behind origin of ray (t < 0).
		return false;										// This avoids having a point on the triangle marked as if it intersect it.

	// Basically we require t to be extrictly positive to consider a possible intersection.
	vec3 q = p + t * d;

	// If reached this point, check if q is inside or outside triangle.
	// Get the vertices from bound object.
	const vec3& v0 = _objectPtr->getVertexPositionAt( _v0 );
	const vec3& v1 = _objectPtr->getVertexPositionAt( _v1 );
	const vec3& v2 = _objectPtr->getVertexPositionAt( _v2 );

	return ( dot( cross( v1 - v0, q - v0 ), _normal ) >= -eps &&
			 dot( cross( v2 - v1, q - v1 ), _normal ) >= -eps &&
			 dot( cross( v0 - v2, q - v2 ), _normal ) >= -eps );
}

