#version 410 core

const int N_BANDS = 4;									// Number of spherical harmonics bands.

in vec3 position;										// Position in world coordinates.

uniform samplerBuffer shCoefficients;					// Spherical harnonics coefficients for current object.
uniform vec3 lightSHCoefficients[N_BANDS * N_BANDS];	// Spherical harmonics coefficients for lighting.

uniform mat4 Model;										// Model transform for rotation and zoom.
uniform mat4 View;										// View matrix takes points from world into camera coordinates.
uniform mat4 Projection;

out vec3 vertexColor;

void main( void )
{
	vec4 p = Model * vec4( position, 1.0 );				// Vertex in true world coordinates.
	gl_Position = Projection * View * p;

	// Color is given by the dot product of the spherical harmonics coefficients of the vertex and the lighting.
	vertexColor = vec3( 0, 0, 0 );
	int NN = N_BANDS * N_BANDS;
	for( int i = 0; i < NN; i++ )
	{
		int offset = ( gl_VertexID * NN + i ) * 3;				// Where vertex coefficients start for current index.
		float r = texelFetch( shCoefficients, offset + 0 ).r;	// Extract RGB components of the vertex SH coefficients.
		float g = texelFetch( shCoefficients, offset + 1 ).r;
		float b = texelFetch( shCoefficients, offset + 2 ).r;

		vertexColor += lightSHCoefficients[i] * vec3( r, g, b );
	}
}
