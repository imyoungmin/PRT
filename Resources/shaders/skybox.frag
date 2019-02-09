#version 410 core

uniform samplerCube skybox;								// Skybox texture.

in vec3 texCoords;

out vec4 color;

/**
 * Main function.
 */
void main( void )
{
	color = texture( skybox, texCoords );
}
