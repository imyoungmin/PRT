#version 410 core

const int N_BANDS = 3;									// Number of spherical harmonics bands.

in vec3 position;										// Position in world coordinates.
in vec3 normal;											// Vertex normal.
in vec3 shCoefficients[N_BANDS * N_BANDS];				// Spherical harmonics coefficients for this vertex (given in RGB).

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
	vertexColor = vec3( 0 );
	for( int i = 0; i < N_BANDS * N_BANDS; i++ )
		vertexColor += shCoefficients[i] * lightSHCoefficients[i];
}
