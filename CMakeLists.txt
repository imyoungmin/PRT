cmake_minimum_required(VERSION 3.13)
project(RTRendering)

set(CMAKE_CXX_STANDARD 14)

add_executable(RTRendering application.cpp
        ArcBall/Ball.h ArcBall/Ball.cpp ArcBall/BallAux.h ArcBall/BallAux.cpp ArcBall/BallMath.h ArcBall/BallMath.cpp
        OpenGL.h OpenGL.cpp
        OpenGLGeometry.h OpenGLGeometry.cpp
        Shaders.h Shaders.cpp
        Atlas.h Atlas.cpp
        Configuration.h
        Object3D.h Object3D.cpp
        Transformations.h Transformations.cpp
		Light.h Light.cpp
		stb_image.h stb_image.cpp)

target_link_libraries(RTRendering
        "-framework OpenGL"
        "armadillo"
        "freetype"
        "glfw")

target_include_directories(RTRendering PUBLIC "/usr/local/include/"
        "/usr/local/include/freetype2/")
