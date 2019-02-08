#version 410 core

uniform vec4 lightPosition;								// In camera coordinates.
uniform vec3 lightColor;								// Only RGB.

uniform vec4 ambient, diffuse, specular;				// The [r,g,b,a] ambient, diffuse, and specular material properties, respectively.
uniform float shininess;
uniform bool useBlinnPhong;
uniform bool useTexture;
uniform bool drawPoint;

uniform sampler2D objectTexture;						// 3D object texture.

in vec3 vPosition;										// Position in view (camera) coordinates.
in vec3 vNormal;										// Normal vector in view coordinates.
in vec2 oTexCoords;

out vec4 color;

/**
 * Apply color given a selected light and shadow map.
 * @param lightColor RGB color of light source.
 * @param lightPosition 3D coordinates of light source with respect to the camera.
 * @param N Normalized normal vector to current fragment (if using Blinn-Phong shading) in camera coordinates.
 * @param E Normalized view direction (if using Blinn-Phong shading) in camera coordinates.
 * @return Fragment color (minus ambient component).
 */
vec3 shade( vec3 lightColor, vec3 lightPosition, vec3 N, vec3 E )
{
	vec3 diffuseColor = diffuse.rgb,
		 specularColor = specular.rgb;
	float shadow;
	
	if( useBlinnPhong )
	{
		vec3 L = normalize( lightPosition - vPosition );
		
		vec3 H = normalize( L + E );
		float incidence = dot( N, L );
		
		// Diffuse component.
		float cDiff = max( incidence, 0.0 );
		diffuseColor = cDiff * ( (useTexture)? texture( objectTexture, oTexCoords ).rgb * diffuseColor : diffuseColor );
		
		// Specular component.
		if( incidence > 0 && shininess > 0.0 )		// Negative shininess turns off specular component.
		{
			float cSpec = pow( max( dot( N, H ), 0.0 ), shininess );
			specularColor = cSpec * specularColor;
		}
		else
			specularColor = vec3( 0.0, 0.0, 0.0 );
	}
	else
		specularColor = vec3( 0.0, 0.0, 0.0 );
	
	// Fragment color with respect to this light (excluding ambient component).
	return ( diffuseColor + specularColor ) * lightColor;
}

/**
 * Main function.
 */
void main( void )
{
	vec3 ambientColor = ambient.rgb;		// Ambient component is constant across lights.
    float alpha = ambient.a;
	vec3 N, E;								// Unit-length normal and eye direction (only necessary for shading with Blinn-Phong reflectance model).
	
	if( useBlinnPhong )
	{
		N = normalize( vNormal );
		E = normalize( -vPosition );
	}
	
    // Final fragment color is the sum of light contributions.
    vec3 totalColor = ambientColor + shade( lightColor, lightPosition.xyz, N, E );
	
    if( drawPoint )
    {
        if( dot( gl_PointCoord - 0.5, gl_PointCoord - 0.5 ) > 0.25 )		// For rounded points.
        	discard;
        else
            color = vec4( totalColor, alpha );
    }
    else
        color = vec4( totalColor, alpha );
}
