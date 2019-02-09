#version 410 core

in vec3 position;

uniform mat4 Model;						// Model transform takes points from model into world coordinates.
uniform mat4 View;						// View matrix takes points from world into camera coordinates.
uniform mat4 Projection;

out vec3 texCoords;						// Interpolate texture coordinates into fragment shader.

void main( void )
{
	vec4 p = Model * vec4( position, 1.0 );				// Vertex in world coordinates.
	vec4 pos = Projection * View * p;
	gl_Position = pos.xyww;								// Trick shader to always give a depth of 1 to cube fragments.
	texCoords = position;								// Vertex becomes a direction for the skybox cubemap.
}