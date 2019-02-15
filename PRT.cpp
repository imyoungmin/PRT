//
// Created by Im YoungMin on 2019-02-08.
//

#include "PRT.h"

using namespace prt;

/**
 * Default constructor.
 */
PRT::PRT() = default;

/**
 * Initialize precomputed radiance transfer object.
 * @param nSamples Number of samples.
 * @param facesFileNames File names of cube map faces.
 * @param program Rendering program.
 * @param nBands Number of spherical-harmonics bands.
 */
void PRT::init( size_t nSamples, const vector<string>& facesFileNames, GLuint program, unsigned int nBands )
{
	_N_SAMPLES = static_cast<size_t>( pow( ceil( sqrt( nSamples ) ), 2.0 ) );		// Make sure the number of samples is a square number.
	_N_BANDS = nBands;
	_generateSamples();
	_generateCubeMap( facesFileNames );
	_renderingProgram = program;

	_lightCoefficients = new vec3[ _N_BANDS * _N_BANDS ];							// Allocate memory for the _N_BANDS^2 lighting projection RGB coefficients.
	_lightCoefficientsArray = new float[ _N_BANDS * _N_BANDS * ELEMENTS_PER_VERTEX ];
}

/**
 * Read and store the 6 faces of the environment light as a cube map.
 * @param facesFileNames
 */
void PRT::_generateCubeMap( const vector<string>& facesFileNames )
{
	if( facesFileNames.size() != 6 )
	{
		cerr << "A cube map expects 6 faces, only " << facesFileNames.size() << " were provided" << endl;
		exit( EXIT_FAILURE );
	}

	// Read cube map 6 faces.
	int width = -1, height = -1;
	int nrChannels = -1;
	int faceIndex = 0;
	for( const string& fileName : facesFileNames )
	{
		string faceFullFileName = conf::SKYBOXES_FOLDER + fileName;
		// Don't flip vertical axis here: cubemaps expect image coordinates to start top left and +y to grow downwards.
		unsigned char *data = stbi_load( faceFullFileName.c_str(), &width, &height, &nrChannels, 0 );
		if( data )
		{
			if( width != height )						// Check for square cube map faces.
			{
				cerr << "Cube map face " << fileName << " is not square: " << width << "x" << height << endl;
				exit( EXIT_FAILURE );
			}

			if( _cubeMapFaceWidth < 0 )					// First image?
			{
				_cubeMapFaceWidth = width;
				_cubeMapFaceNrChannels = nrChannels;
			}

			if( width != _cubeMapFaceWidth || nrChannels != _cubeMapFaceNrChannels )		// Check for compatible face sizes.
			{
				cerr << "Cube map face " << fileName << " has an incompatible size/number of channels" << endl;
				exit( EXIT_FAILURE );
			}

			_cubeMapFaces[faceIndex] = data;			// Keep data for later projection of "light" onto the sample spherical harmonics functions.
			faceIndex++;
		}
		else
		{
			cerr << "Cube map face failed to load at path: " << faceFullFileName << endl;
			exit( EXIT_FAILURE );
		}
	}
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

			// Generate the _N_BANDS^2 spherical-harmonics functions for this sample.
			vector<double> sh;
			for( int l = 0; l < _N_BANDS; l++ )
			{
				for( int m = -l; m <= l; m++ )
					sh.emplace_back( _y_lm( l, m, theta, phi ) );					// The ith function evaluation is at i = l(l+1) + m.
			}

			_samples.emplace_back( theta, phi, sh );								// Add new sample.
		}
	}
}

/**
 * Query pixel value from the environment cube map.
 * @param query Direction.
 * @param output Pixel values for each channel registered for the cube map's faces.
 */
void PRT::_queryCubeMap( const vec3& query, unsigned char* output ) const
{
	// First, choose which face to query based on the query direction.
	vec3 d = normalise( query );
	unsigned int chosenFace = 6;
	double chosenFaceCosine = -2;			// Select the face whose normal vector has the smallest angle with the query direction.
	for( unsigned int i = 0; i < 6; i++ )
	{
		double cosine = dot( _cubeMapFacesNormals[i], d );
		if( cosine > chosenFaceCosine )
		{
			chosenFace = i;
			chosenFaceCosine = cosine;
		}
	}

	// Now that we have the face to query, find intersection of query direction (as a ray) with the chosen face plane.
	vec3 p;
	double t;
	double xf, yf;							// The x- and y-coordinates on the chose face at the intersection point, in [0,1].
	switch( chosenFace )
	{
		case 0:								// Case RIGHT: plane equation is x = +1.
			t = 1.0 / d[0];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[2] + 1.0 ) / 2.0;
			yf = ( p[1] - 1.0 ) / -2.0;
			xf = 1.0 - xf;					// Adjustment.
			break;
		case 1:								// Case LEFT: plane equation is x = -1.
			t = -1.0 / d[0];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[2] - 1.0 ) / -2.0;
			yf = ( p[1] - 1.0 ) / -2.0;
			xf = 1.0 - xf;					// Adjustment.
			break;
		case 2:								// Case TOP: plane equation is y = +1.
			t = 1.0 / d[1];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[0] + 1.0 ) / 2.0;
			yf = ( p[2] - 1.0 ) / -2.0;
			yf = 1.0 - yf;					// Adjustment.
			break;
		case 3:								// Case BOTTOM: plane equation is y = -1.
			t = -1.0 / d[1];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[0] + 1.0 ) / 2.0;
			yf = ( p[2] + 1.0 ) / 2.0;
			yf = 1.0 - yf;					// Adjustment.
			break;
		case 4:								// Case FRONT: plane equation is z = +1.
			t = 1.0 / d[2];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[0] - 1.0 ) / -2.0;
			yf = ( p[1] - 1.0 ) / -2.0;
			xf = 1.0 - xf;					// Adjustment.
			break;
		case 5:								// Case BACK: plane equation is z = -1.
			t = -1.0 / d[2];
			p = t * d;						// Stretch query direction to intersect face plane.
			xf = ( p[0] + 1.0 ) / 2.0;
			yf = ( p[1] - 1.0 ) / -2.0;
			xf = 1.0 - xf;					// Adjustment.
			break;
		default:
			cerr << "Attempting to query a nonexistent face in the cube map" << endl;
			exit( EXIT_FAILURE );
	}

	// Discrete pixel coordinates.
	auto xd = static_cast<unsigned>( min( floor( _cubeMapFaceWidth * xf ), _cubeMapFaceWidth - 1.0 ) );
	auto yd = static_cast<unsigned>( min( floor( _cubeMapFaceWidth * yf ), _cubeMapFaceWidth - 1.0 ) );
	return _getPixel( xd, yd, chosenFace, output );
}

/**
 * Read pixel values from a cube map face image.
 * @param x Column index (left - right).
 * @param y Row index (top - bottom).
 * @param face Which face to query.
 * @param output Result pixel values (triplet, etc), depending on number of channels, in integer format.
 */
void PRT::_getPixel( unsigned int x, unsigned int y, unsigned int face, unsigned char* output ) const
{
	face = ( face > 5 )? 5 : face;										// Check boundaries.
	x = ( x > _cubeMapFaceWidth - 1 )? _cubeMapFaceWidth - 1 : x;
	y = ( y > _cubeMapFaceWidth - 1 )? _cubeMapFaceWidth - 1 : y;

	int pixelPosition = ( y * _cubeMapFaceWidth + x ) * _cubeMapFaceNrChannels;
	for( int i = 0; i < _cubeMapFaceNrChannels; i ++ )					// Fill in pixel value for each channel.
		output[i] = _cubeMapFaces[face][pixelPosition + i];				// Caller must allocate space in output array.
}

/**
 * Project the sampled environment lighting onto the _N_BANDS^2 spherical harmonics functions.
 */
void PRT::_projectLighting()
{
	// Initialize lighting coefficients.
	for( int i = 0; i < _N_BANDS * _N_BANDS; i++ )
		_lightCoefficients[i] = { 0, 0, 0 };

	// Projection: accumulation process.
	unsigned char uLight[_cubeMapFaceNrChannels];
	for( const Sample& sample : _samples )					// For each sample, compute the light projection coefficients.
	{
		_queryCubeMap( sample.getPosition(), uLight );		// Light color at the sample position (received from cube map).
		vec3 light = { uLight[0] / 255.0, uLight[1] / 255.0, uLight[2] / 255.0 };
		for( int i = 0; i < _N_BANDS * _N_BANDS; i++ )
			_lightCoefficients[i] += light * sample.getSHValueAt( i );
	}

	// Projection: final scaling.
	// c_k \approx \frac{4\pi}{n} \sum_{j=1}^{n}( f(\omega_j) y_k(\omega_j) ).
	for( int i = 0; i < _N_BANDS * _N_BANDS; i++ )
		_lightCoefficients[i] *= 4.0 * M_PI / _samples.size();

	// Flatten (copy) coefficients into a linear array to be sent to shaders.
	int l = 0;
	for( int i = 0; i < _N_BANDS * _N_BANDS; i++ )
	{
		_lightCoefficientsArray[l] = 1/*_lightCoefficients[i][0]*/; l++;		// R.
		_lightCoefficientsArray[l] = 1/*_lightCoefficients[i][1]*/; l++;		// G.
		_lightCoefficientsArray[l] = 1/*_lightCoefficients[i][2]*/; l++;		// B.
	}
}

/**
 * Project the objects' vertices onto the _N_BANDS^2 spherical harmonics functions.
 * This is for the unshadow diffuse transfer function.
 */
void PRT::_unshadowedDiffuseTransferProjection()
{
	for( Object3D& object : _objects )
	{
		object.resetSHCoefficients();
		const vec3& color = object.getColor();

		// Diffuse transfer projection at each vertex.
		for( unsigned int i = 0; i < object.getVerticesCount(); i++ )
		{
			const vec3& n = object.getVertexNormalAt( i );			// Vertex normal.
			for( const Sample& sample : _samples )					// For each (\theta, \phi) direction.
			{
				double h = max( 0.0, dot( n, sample.getPosition() ) );		// Geometric function.
				if( h > 0 )
				{
					for( unsigned int j = 0; j < _N_BANDS * _N_BANDS; j++ )
						object.accumulateSHCoefficients( i, j, sample.getSHValueAt( j ) * h * color );
				}
			}
		}

		// Projection: final scaling.
		// c_k \approx \frac{4\pi}{n} \sum_{j=1}^{n}( f(\omega_j) y_k(\omega_j) ).
		object.scaleSHCoefficients( 4.0 * M_PI / _samples.size() );

		// Load computed coefficients into a texture.
		object.loadSHCoefficientsIntoTexture();
	}
}

/**
 * Spherical harmonics function.
 * @param l Band index.
 * @param m Offset index within a band.
 * @param theta Unit sphere theta angle (longitude).
 * @param phi Unit sphere phi angle (latitude).
 * @return Function value.
 */
double PRT::_y_lm( int l, int m, double theta, double phi )
{
	// We have three cases: m > 0, m < 0, and m = 0.
	if( m > 0 )
		return M_SQRT2 * _K_lm( l, m ) * cos( m * phi ) * _P_lm( l, m, cos( theta ) );

	if( m < 0 )
		return M_SQRT2 * _K_lm( l, m ) * sin( -m * phi ) * _P_lm( l, -m, cos( theta ) );

	return _K_lm( l, 0 ) * _P_lm( l, 0, cos( theta ) );
}

/**
 * Associated Legendre polynomial.
 * https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 * @param l Band index.
 * @param m Offset index within a band.
 * @param x Input value to polynomial.
 * @return Associated Legendre polynomial evaluated at x.
 */
double PRT::_P_lm( int l, int m, double x )
{
	if( l == m )
		return pow( -1, l ) * _doubleFactorial( 2 * l - 1 ) * pow( 1 - x * x, static_cast<double>( l )/2.0 );

	if( l == m + 1 )
		return x * ( 2 * m + 1 ) * _P_lm( l - 1, m, x );

	return ( x * ( 2 * l - 1 ) * _P_lm( l - 1, m, x ) - ( l + m - 1) * _P_lm( l - 2, m, x ) ) / ( l - m );
}

/**
 * Spherical harmonics normalizing constant.
 * @param l Band index.
 * @param m Offset index within a band.
 * @return Normalizing constant.
 */
double PRT::_K_lm( int l, int m )
{
	if( m == 0 )
		return sqrt( ( 2.0 * l + 1.0 ) / ( 4.0 * M_PI ) );

	return sqrt( ( 2.0 * l + 1.0 ) / ( 4.0 * M_PI ) *
				static_cast<double>( _factorial( l - abs( m ) ) ) / static_cast<double>( _factorial( l + abs( m ) ) ) );
}

/**
 * Compute the factorial of a non-negative integer.
 * @param x Input value.
 * @return x!
 */
int PRT::_factorial( int x )
{
	if( x <= 1)
		return 1;

	int factorial = 1;
	for( int i = 2; i <= x; i++ )
		factorial *= i;
	return factorial;
}

/**
 * Double factorial of an integer.
 * https://en.wikipedia.org/wiki/Double_factorial
 * @param x Input value.
 * @return x!!
 */
int PRT::_doubleFactorial( int x )
{
	if( x <= 1 )
		return 1;
	return x * _doubleFactorial( x - 2);
}

/**
 * Destructor.
 */
PRT::~PRT()
{
	// Deallocate environment lighting image data (i.e. cubemap faces).
	for( auto data : _cubeMapFaces )
	{
		if( data != nullptr )
			stbi_image_free( data );
	}

	// Deallocate lighting projection coefficients.
	delete [] _lightCoefficients;
	delete [] _lightCoefficientsArray;
}


/**
 * Retrieve a cube map's face data.
 * @param faceIndex Index in [0,6].
 * @return Image data.
 */
const unsigned char* PRT::getCubeMapFaceData( int faceIndex ) const
{
	return _cubeMapFaces[faceIndex];
}

/**
 * Get the width (or side length of the square) of any of the cube map's faces.
 * @return Square side length.
 */
int PRT::getCubeMapFaceWidth() const
{
	return _cubeMapFaceWidth;
}

/**
 * Get the number of channels in the cube map's face images.
 * @return Number of channels, e.g. 3 for RGB, 4 for RGBA.
 */
int PRT::getCubeMapFaceNrChannels() const
{
	return _cubeMapFaceNrChannels;
}

/**
 * Add a shading 3D object to the PRT scene.
 * @param vertices List of 3D vertices in object coordinates.
 * @param normals List of normal vectors in object coordinates.
 * @param T Transformation matrix that takes object coordinates to world coordinates.
 * @param color RGB object color (currently transparency is non supported).
 */
void PRT::addObject( const vector<vec3>& vertices, const vector<vec3>& normals, const mat44& T, const vec3& color )
{
	// New object, with _N_BANDS^2 spherical harmonics projection coefficients per vertex.
	// Note the use of emplace_back to avoid the extra copy of push_back.
	_objects.emplace_back( vertices, normals, T, color, _N_BANDS * _N_BANDS );
}

/**
 * Execute the precomputation process.
 */
void PRT::precomputeRadianceTransfer()
{
	cout << "[PRT] Now projecting sampled lighting... ";
	_projectLighting();
	cout << "Done!" << endl;

	cout << "[PRT] Now projecting unshadowed diffuse transfer function... ";
	_unshadowedDiffuseTransferProjection();
	cout << "Done!" << endl;
}

/**
 * Render scene objects.
 * @param Projection The 4x4 projection matrix.
 * @param Camera The 4x4 view matrix.
 * @param Model The 4x4 (interaction) model matrix.
 */
void PRT::renderObjects( const mat44& Projection, const mat44& Camera, const mat44& Model )
{
	// Send constant/uniform information to shaders.
	int model_location = glGetUniformLocation( _renderingProgram, "Model" );						// Transformation matrices.
	int view_location = glGetUniformLocation( _renderingProgram, "View");
	int proj_location = glGetUniformLocation( _renderingProgram, "Projection" );

	int lightCoeffs_location = glGetUniformLocation( _renderingProgram, "lightSHCoefficients" );	// And light spherical harmonics coefficients.

	if( model_location != -1 && view_location != -1 && proj_location != -1 && lightCoeffs_location != -1 )
	{
		float model_matrix[ELEMENTS_PER_MATRIX];
		Tx::toOpenGLMatrix( model_matrix, Model );
		glUniformMatrix4fv( model_location, 1, GL_FALSE, model_matrix );		// Send model matrix.

		float view_matrix[ELEMENTS_PER_MATRIX];
		Tx::toOpenGLMatrix( view_matrix, Camera );
		glUniformMatrix4fv( view_location, 1, GL_FALSE, view_matrix );			// Send view matrix.

		float proj_matrix[ELEMENTS_PER_MATRIX];
		Tx::toOpenGLMatrix( proj_matrix, Projection );
		glUniformMatrix4fv( proj_location, 1, GL_FALSE, proj_matrix );			// Send projection matrix.

		glUniform3fv( lightCoeffs_location, _N_BANDS * _N_BANDS, _lightCoefficientsArray );
	}

	for( const Object3D& o : _objects )
	{
		glBindBuffer( GL_ARRAY_BUFFER, o.getBufferID() );

		// Set up our vertex and spherical harmonics coefficients attributes and uniforms.
		GLint position_location = glGetAttribLocation( _renderingProgram, "position" );
		GLint shCoefficients_location = glGetUniformLocation( _renderingProgram, "shCoefficients" );
		if( position_location != -1  && shCoefficients_location != -1 )
		{
			glEnableVertexAttribArray( static_cast<GLuint>( position_location ) );	// Send vertices.
			glVertexAttribPointer( static_cast<GLuint>( position_location ), ELEMENTS_PER_VERTEX, GL_FLOAT, GL_FALSE, 0, nullptr );

			glActiveTexture( GL_TEXTURE0 );											// Send the spherical harmonics coefficients.
			glBindTexture( GL_TEXTURE_BUFFER, o.getTBOTextureID() );
			glTexBuffer( GL_TEXTURE_BUFFER, GL_R32F, o.getTBOID() );
			glUniform1i( shCoefficients_location, 0 );

			// Draw triangles.
			glDrawArrays( GL_TRIANGLES, 0, o.getVerticesCount() );

			// Disable attribute arrays.
			glDisableVertexAttribArray( static_cast<GLuint>( position_location ) );
		}
	}
}

/////////////////////////////////////////////////// Sample class ///////////////////////////////////////////////////////

/**
 * Sample constructor.
 * @param t Spherical theta coordinate.
 * @param p Spherical phi coordinate.
 * @param sh Spherical harmonics function evaluations.
 */
Sample::Sample( double theta, double phi, const vector<double>& sh )
{
	_theta = theta;			// Spherical coordinates.
	_phi = phi;
	_sh = vector<double>( sh );

	// Position in rectangular coordinates.
	_position = { sin( _theta ) * cos( _phi ), sin( _theta ) * sin( _phi ), cos( _theta ) };
}

/**
 * Get the spherical harmonics value at a given index.
 * @param index Querying index.
 * @return Spherical harmonics functional value.
 */
double Sample::getSHValueAt( int index ) const
{
	return _sh[index];
}

/**
 * Get sample position.
 * @return 3D position.
 */
const vec3& Sample::getPosition() const
{
	return _position;
}