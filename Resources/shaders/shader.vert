#version 410 core

in vec3 position;
in vec3 normal;
in vec2 texCoords;

uniform mat4 Model;										// Model transform takes points from model into world coordinates.
uniform mat4 View;										// View matrix takes points from world into camera coordinates.
uniform mat3 InvTransModelView;							// Inverse-transposed 3x3 principal submatrix of ModelView matrix.
uniform mat4 Projection;
uniform float pointSize;
uniform bool useBlinnPhong;

out vec3 vPosition;										// Position in view (camera) coordinates.
out vec3 vNormal;										// Normal vector in view coordinates.
out vec2 oTexCoords;									// Interpolate texture coordinates into fragment shader.

void main( void )
{
	vec4 p = Model * vec4( position, 1.0 );				// Vertex in world coordinates.
	gl_Position = Projection * View * p;

	if( useBlinnPhong )
	{
		vPosition = (View * p).xyz;						// Send vertex and normal to fragment shader in camera coodinates.
		vNormal = InvTransModelView * normal;
	}

	gl_PointSize = pointSize;
	oTexCoords = texCoords;
}
