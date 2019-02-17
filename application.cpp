/**
 * OpenGL main application.
 */

#include <iostream>
#include <armadillo>
#include <OpenGL/gl3.h>
#include <string>
#include "GLFW/glfw3.h"
#include "ArcBall/Ball.h"
#include "OpenGL.h"
#include "Transformations.h"
#include "PRT.h"

using namespace std;
using namespace std::chrono;
using namespace arma;

// Perspective projection matrix.
mat44 Proj;

// Text scaling.
float gTextScaleX;
float gTextScaleY;

// Camera controls globals.
vec3 gPointOfInterest;
vec3 gEye;
vec3 gUp;

bool gLocked;							// Track if mouse button is pressed down.
bool gUsingArrowKey;					// Track if we are using the arrow keys for rotating scene.
bool gRotatingLights;					// Enable/disable rotating lights about the scene.
float gZoom;							// Camera zoom.
const float ZOOM_IN = 1.015;
const float ZOOM_OUT = 0.985;
BallData* gArcBall;						// Arc ball.

// Framebuffer size metrics.
int fbWidth;
int fbHeight;
float gRetinaRatio;						// How many screen dots exist per OpenGL pixel.

OpenGL ogl;								// Initialize application OpenGL.

prt::PRT gPRT;							// Precomputed radiance transfer object.

// Lights.
vector<Light> gLights;					// Light source objects.
int gLightsCount = 0;

// Frame rate variables and functions.
static const int NUM_FPS_SAMPLES = 64;
float gFpsSamples[NUM_FPS_SAMPLES];
unsigned char gCurrentSample = 0;		// Should start storing from gCurrentSample >= 1.

/**
 * Calculate the number of frames per second using a window.
 * @param dt Amount of seconds for current frame.
 * @return Frames per second.
 */
float calculateFPS( float dt )
{
	gCurrentSample++;
	gCurrentSample = static_cast<unsigned char>( max( 1, static_cast<int>( gCurrentSample ) ) );
	if( dt <= 0 )
		cout << "error" << endl;
	gFpsSamples[(gCurrentSample - 1) % NUM_FPS_SAMPLES] = 1.0f / dt;
	float fps = 0;
	int i = 0;
	for( i = 0; i < min( NUM_FPS_SAMPLES, static_cast<int>( gCurrentSample ) ); i++ )
		fps += gFpsSamples[i];
	fps /= i;
	return fps;
}

/**
 * Reset rotation and zoom.
 */
void resetArcBall()
{
	Ball_Init( gArcBall );
	Ball_Place( gArcBall, qOne, 0.75 );
}

/**
 * GLFW error callback.
 * @param error Error code.
 * @param description Error description.
 */
void errorCallback( int error, const char* description )
{
	cerr << error << ": " << description << endl;
}

/**
 * Rotate scene in x or y direction.
 * @param x Rotation amount in x direction, usually in the range [-1,1].
 * @param y Rotation amount in y direction, usually in the range [-1,1].
 */
void rotateWithArrowKey( float x, float y )
{
	if( gLocked )		// Do not rotate scene with arrow key if it's currently rotating with mouse.
		return;

	gUsingArrowKey = true;						// Start blocking rotation with mouse button.
	HVect arcballCoordsStart, arcballCoordsEnd;

	arcballCoordsStart.x = 0.0;					// Start rotation step.
	arcballCoordsStart.y = 0.0;
	Ball_Mouse( gArcBall, arcballCoordsStart );
	Ball_BeginDrag( gArcBall );

	arcballCoordsEnd.x = x;						// End rotation step.
	arcballCoordsEnd.y = y;
	Ball_Mouse( gArcBall, arcballCoordsEnd );
	Ball_Update( gArcBall );
	Ball_EndDrag( gArcBall );
	gUsingArrowKey = false;						// Exiting conflicting block for rotating with arrow keys.
}

/**
 * GLFW keypress callback.
 * @param window GLFW window.
 * @param key Which key was pressed.
 * @param scancode Unique code for key.
 * @param action Key action: pressed, etc.
 * @param mods Modifier bits: shift, ctrl, alt, super.
 */
void keyCallback( GLFWwindow* window, int key, int scancode, int action, int mods )
{
	if( action != GLFW_PRESS && action != GLFW_REPEAT )
		return;
	
	const float rotationStep = 0.0025;
	
	switch( key )
	{
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose( window, GL_TRUE );
			break;
		case GLFW_KEY_LEFT:
			rotateWithArrowKey( -rotationStep, 0 );
			break;
		case GLFW_KEY_RIGHT:
			rotateWithArrowKey( +rotationStep, 0 );
			break;
		case GLFW_KEY_UP:
			rotateWithArrowKey( 0, +rotationStep );
			break;
		case GLFW_KEY_DOWN:
			rotateWithArrowKey( 0, -rotationStep );
			break;
		case GLFW_KEY_R:
			resetArcBall();
			gZoom = 1.0;
			break;
		case GLFW_KEY_L:
			gRotatingLights = !gRotatingLights;
			break;
		default: return;
	}
}

/**
 * GLFW mouse button callback.
 * @param window GLFW window.
 * @param button Which mouse button has been actioned.
 * @param action Mouse action: pressed or released.
 * @param mode Modifier bits: shift, ctrl, alt, supper.
 */
void mouseButtonCallback( GLFWwindow* window, int button, int action, int mode )
{
	if( button != GLFW_MOUSE_BUTTON_LEFT )		// Ignore mouse button other than left one.
		return;

	if( gUsingArrowKey )						// Wait for arrow keys to stop being used as rotation mechanism.
		return;

	if( action == GLFW_PRESS )
	{
		//glfwSetInputMode( window, GLFW_CURSOR, GLFW_CURSOR_DISABLED );
		HVect arcballCoords;
		double x, y;
		int w, h;
		glfwGetWindowSize( window, &w, &h );
		glfwGetCursorPos( window, &x, &y );
		arcballCoords.x = static_cast<float>( 2.0*x/static_cast<float>(w) - 1.0 );
		arcballCoords.y = static_cast<float>( -2.0*y/static_cast<float>(h) + 1.0 );
		Ball_Mouse( gArcBall, arcballCoords );
		Ball_BeginDrag( gArcBall );
		gLocked = true;							// Determines whether the mouse movement is used for rotating the object.
	}
	else
	{
		//glfwSetInputMode( window, GLFW_CURSOR, GLFW_CURSOR_NORMAL );
		Ball_EndDrag( gArcBall );
		gLocked = false;						// Go back to normal.
	}
}

/**
 * GLFW mouse motion callback.
 * @param window GLFW window.
 * @param x Mouse x-position.
 * @param y Mouse y-position.
 */
void mousePositionCallback( GLFWwindow* window, double x, double y )
{
	if( glfwGetMouseButton( window, GLFW_MOUSE_BUTTON_LEFT ) == GLFW_PRESS && gLocked )
	{
		HVect arcballCoords;
		int w, h;
		glfwGetWindowSize( window, &w, &h );
		arcballCoords.x = static_cast<float>( 2.0*x/static_cast<float>(w) - 1.0 );
		arcballCoords.y = static_cast<float>( -2.0*y/static_cast<float>(h) + 1.0 );
		Ball_Mouse( gArcBall, arcballCoords );
		Ball_Update( gArcBall );
	}
}

/**
 * GLFW mouse scroll callback.
 * @param window GLFW window.
 * @param xOffset X offset.
 * @param yOffset Y offset.
 */
void mouseScrollCallback( GLFWwindow* window, double xOffset, double yOffset )
{
	gZoom *= (yOffset > 0)? ZOOM_IN: ZOOM_OUT;
	gZoom = max( 0.5f, min( gZoom, 2.5f ) );
}

/**
 * GLFW frame buffer resize callback.
 * @param window GLFW window.
 * @param w New frame buffer width.
 * @param h New frame buffer height.
 */
void resizeCallback( GLFWwindow* window, int w, int h )
{
	fbWidth = w;		// w and h are width and height of framebuffer, not window.
	fbHeight = h;

	//Proj = Tx::frustrum( -0.5, 0.5, -0.5, 0.5, 1.0, 100 );

	// Projection used for 3D.
	double ratio = static_cast<double>(w)/static_cast<double>(h);
	Proj = Tx::perspective( M_PI/3.0, ratio, 0.01, 1000.0 );

	// Projection used for text rendering.
	int windowW, windowH;
	glfwGetWindowSize( window, &windowW, &windowH );
	gTextScaleX = 1.0f / windowW;
	gTextScaleY = 1.0f / windowH;
}

/**
 * Load skybox cubemap texture.
 * Texture targets are as follows:
 * GL_TEXTURE_CUBE_MAP_POSITIVE_X	Right
 * GL_TEXTURE_CUBE_MAP_NEGATIVE_X	Left
 * GL_TEXTURE_CUBE_MAP_POSITIVE_Y	Top
 * GL_TEXTURE_CUBE_MAP_NEGATIVE_Y	Bottom
 * GL_TEXTURE_CUBE_MAP_POSITIVE_Z	Front
 * GL_TEXTURE_CUBE_MAP_NEGATIVE_Z	Back
 * @return Texture ID.
 */
GLuint loadCubemap()
{
	GLuint textureID;
	glGenTextures( 1, &textureID );
	glBindTexture( GL_TEXTURE_CUBE_MAP, textureID );

	int width = gPRT.getCubeMapFaceWidth(), nrChannels = gPRT.getCubeMapFaceNrChannels();
	for( int i = 0; i < 6; i++ )
	{
		// Read data proloaded in PRT object.
		const unsigned char *data = gPRT.getCubeMapFaceData( i );
		if( data )
			glTexImage2D( static_cast<GLenum>(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i), 0, GL_RGB, width, width, 0, (nrChannels == 4)? GL_RGBA: GL_RGB, GL_UNSIGNED_BYTE, data );
		else
		{
			cerr << "Cubemap texture failed to load for face " << i << endl;
			exit( EXIT_FAILURE );
		}
	}
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE );	// Third texture dimension.

	return textureID;
}

/**
 * Application main function.
 * @param argc Number of input arguments.
 * @param argv Input arguments.
 * @return Exit code.
 */
int main( int argc, const char * argv[] )
{
	srand( static_cast<unsigned>( time( nullptr ) ) );
	
	gPointOfInterest = { 0, 1, 0 };		// Camera controls globals.
	gEye = { 6, 5, 8 };
	gUp = Tx::Y_AXIS;
	
	gLocked = false;					// Track if mouse button is pressed down.
	gRotatingLights = false;			// Start with still lights.
	gUsingArrowKey = false;				// Track pressing action of arrow keys.
	gZoom = 1.0;						// Camera zoom.
	
	GLFWwindow* window;
	glfwSetErrorCallback( errorCallback );
	
	if( !glfwInit() )
		exit( EXIT_FAILURE );
	
	// Indicate GLFW which version will be used and the OpenGL core only.
	glfwWindowHint( GLFW_SAMPLES, 4 );
	glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 4 );
	glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 1 );
	glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
	glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE );
	
	cout << glfwGetVersionString() << endl;

	// Create window object (with screen-dependent size metrics).
	int windowWidth = 1080;
	int windowHeight = 720;
	window = glfwCreateWindow( windowWidth, windowHeight, "Real-Time Rendering", nullptr, nullptr );

	if( !window )
	{
		glfwTerminate();
		exit( EXIT_FAILURE );
	}
	
	glfwMakeContextCurrent( window );
	glfwSwapInterval( 1 );
	
	// Hook up callbacks.
	glfwSetFramebufferSizeCallback( window, resizeCallback );
	glfwSetKeyCallback( window, keyCallback );
	glfwSetMouseButtonCallback( window, mouseButtonCallback );
	glfwSetCursorPosCallback( window, mousePositionCallback );
	glfwSetScrollCallback( window, mouseScrollCallback );

	// Initialize projection matrices and viewport.
	glfwGetFramebufferSize( window, &fbWidth, &fbHeight );
	gRetinaRatio = static_cast<float>( fbWidth ) / windowWidth;
	cout << "Retina pixel ratio: " << gRetinaRatio << endl;
	resizeCallback( window, fbWidth, fbHeight );

	// Check maximum texture buffer size.
	int limit;
	glGetIntegerv( GL_MAX_TEXTURE_BUFFER_SIZE, &limit );
	cout << "Maximum texture buffer size: " << limit << endl;
	
	gArcBall = new BallData;						// Initialize arcball.
	resetArcBall();
	
	///////////////////////////////////// Intialize OpenGL and rendering shaders ///////////////////////////////////////
	
	ogl.init();
	Shaders shaders;

	// Initialize shaders for skybox program.
	cout << "Initializing skybox shaders... ";
	GLuint skyboxProgram = shaders.compile( conf::SHADERS_FOLDER + "skybox.vert", conf::SHADERS_FOLDER + "skybox.frag" );
	cout << "Done!" << endl;

	// Intialize shaders for precomputed radiance transfer.
	cout << "Initializing PRT shaders... ";
	GLuint prtProgram = shaders.compile( conf::SHADERS_FOLDER + "prt.vert", conf::SHADERS_FOLDER + "prt.frag" );
	cout << "Done!" << endl;

	///////////////////////////////// Initialize precomputed radiance transfer object //////////////////////////////////

	vector<string> faces = {
			"skybox1/right.tga",
			"skybox1/left.tga",
			"skybox1/top.tga",
			"skybox1/bottom.tga",
			"skybox1/front.tga",
			"skybox1/back.tga"
	};
	gPRT.init( 10*10, faces, prtProgram );

	// Loading objects' original data to be used in PRT.
	vector<vec3> vertices, normals;
	vector<vec2> uvs;
	Object3D::loadOBJ( "cube.obj", vertices, uvs, normals );
	gPRT.addObject( "CubeTop", vertices, normals, Tx::translate( 0.0, 0.25, 0.0 ) * Tx::scale( 2.0, 0.05, 2.0 ), { 1.0, 1.0, 1.0 } );
	gPRT.addObject( "CubeBottom", vertices, normals, Tx::translate( 0.0, 0.05, 0.0 ) * Tx::scale( 3.0, 0.05, 3.0 ), { 1.0, 1.0, 1.0 } );
//	Object3D::loadOBJ( "bust.obj", vertices, uvs, normals );
//	gPRT.addObject( "Bust", vertices, normals, Tx::scale( 0.75 ) * Tx::rotate( -M_PI/2.0, Tx::X_AXIS ) * Tx::rotate( M_PI, Tx::Z_AXIS ), { 1.0, 1.0, 1.0 } );
	Object3D::loadOBJ( "mask.obj", vertices, uvs, normals );
	gPRT.addObject( "Mask", vertices, normals, Tx::translate( 0.0, 0.3, 0.0 ) * Tx::scale( 0.75 ), { 1.0, 1.0, 1.0 } );
//	Object3D::loadOBJ( "deer.obj", vertices, uvs, normals );
//	gPRT.addObject( "Deer", vertices, normals, Tx::translate( 0.0, 0.2, 0.0 ) * Tx::scale( 0.75 ) * Tx::rotate( M_PI/2.0, Tx::Y_AXIS ), { 1.0, 1.0, 1.0 } );

	gPRT.precomputeRadianceTransfer();

	//////////////////////////////////////////////// Create lights /////////////////////////////////////////////////////
	
	gLightsCount = 1;
	const double lRadius = sqrt( 11 * 11 * 2 );
	const double theta = 2.0 * M_PI / gLightsCount;
	const float phi = static_cast<float>( rand() ) / static_cast<float>( RAND_MAX / M_PI_4 );
	const float lHeight = 15;
	const float lRGB[3] = { 0.8, 0.8, 0.8 };
	gLights.push_back( Light({ lRadius * sin( theta + phi ), lHeight, lRadius * cos( theta + phi ) }, { lRGB[0], lRGB[1], lRGB[2] }, eye( 4, 4 ), 0 ) );

	//////////////////////////////////////////////// Create skyboxes ///////////////////////////////////////////////////

	GLuint skyboxTexture = loadCubemap();
	GLint skybox_sampler_location = glGetUniformLocation( skyboxProgram, "skybox" );
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double currentTime = 0.0;
	const double timeStep = 0.01;
	const float textColor[] = { 0.0, 0.8, 1.0, 1.0 };
	char text[128];
	
	glEnable( GL_DEPTH_TEST );
	glDepthFunc( GL_LEQUAL );
	glFrontFace( GL_CCW );

	// Frame rate variables.
	long gNewTicks = duration_cast<milliseconds>( system_clock::now().time_since_epoch() ).count();
	long gOldTicks = gNewTicks;
	float transcurredTimePerFrame;
	string FPS = "FPS: ";

	ogl.setUsingUniformScaling( true );							// Important! We'll be using uniform scaling in the following scene rendering.

	// Rendering loop.
	while( !glfwWindowShouldClose( window ) )
	{
		glClearColor( 0.0f, 0.0f, 0.01f, 1.0f );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		HMatrix abr;
		Ball_Value( gArcBall, abr );
		mat44 ArcBall = {
			{ abr[0][0], abr[0][1], abr[0][2], abr[0][3] },
			{ abr[1][0], abr[1][1], abr[1][2], abr[1][3] },
			{ abr[2][0], abr[2][1], abr[2][2], abr[2][3] },
			{ abr[3][0], abr[3][1], abr[3][2], abr[3][3] } };
		mat44 Model = ArcBall.t() * Tx::scale( gZoom );
		mat44 Camera = Tx::lookAt( gEye, gPointOfInterest, gUp );

		glViewport( 0, 0, fbWidth, fbHeight );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		
		///////////////////////////////////////// Define new lights' positions /////////////////////////////////////////
		
		if( gRotatingLights )								// Check if rotating lights is enabled (with key 'L').
		{
			for( int i = 0; i < gLightsCount; i++ )
				gLights[i].rotateBy( static_cast<float>( 0.01 * M_PI ) );
		}

		//////////////////////////////////////////////// Render scene //////////////////////////////////////////////////

		glEnable( GL_CULL_FACE );
		glUseProgram( prtProgram );
		gPRT.renderObjects( Proj, Camera, Model );

		/////////////////////////////////////////////// Render skybox //////////////////////////////////////////////////

		glDepthMask( GL_FALSE );											// Disable depth writing.  This way it's always drawn in background.
		glCullFace( GL_FRONT );

		ogl.useProgram( skyboxProgram );
		glActiveTexture( GL_TEXTURE0 );										// Enable texture unit.
		glBindTexture( GL_TEXTURE_CUBE_MAP, skyboxTexture );
		glUniform1i( skybox_sampler_location, 0 );
		mat44 CameraPrime = eye( 4, 4 );									// Obtain the principal 3x3 submatrix of the View transform to remove the translation.
		CameraPrime.submat( 0, 0, 2, 2 ) = Camera.submat( 0, 0, 2, 2 );
		ogl.drawCube( Proj, CameraPrime, Model );							// Cube will be colored with cube texture map.

		glCullFace( GL_BACK );
		glDepthMask( GL_TRUE );

		/////////////////////////////////////////////// Rendering text /////////////////////////////////////////////////

		glUseProgram( ogl.getGlyphsProgram() );				// Switch to text rendering.  The text rendering is the only program created within the OpenGL class.

		glEnable( GL_BLEND );
		glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
		glDisable( GL_CULL_FACE );

		gNewTicks = duration_cast<milliseconds>( system_clock::now().time_since_epoch() ).count();
		transcurredTimePerFrame = (gNewTicks - gOldTicks) / 1000.0f;
		sprintf( text, "FPS: %.2f", ( ( transcurredTimePerFrame <= 0 )? -1 : calculateFPS( transcurredTimePerFrame ) ) );
		gOldTicks = gNewTicks;

		ogl.renderText( text, ogl.atlas48, -1 + 10 * gTextScaleX, 1 - 30 * gTextScaleY, static_cast<float>( gTextScaleX * 0.6 ),
						static_cast<float>( gTextScaleY * 0.6 ), textColor );

		glDisable( GL_BLEND );

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		glfwSwapBuffers( window );
		glfwPollEvents();
		
		currentTime += timeStep;
	}
	
	glfwDestroyWindow( window );
	glfwTerminate();
	
	// Delete OpenGL programs.
	glDeleteProgram( skyboxProgram );
	glDeleteProgram( prtProgram );

	return 0;
}

