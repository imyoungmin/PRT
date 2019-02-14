cmake_minimum_required(VERSION 3.13)

link_directories(/usr/local/lib)

project(PRT)

set(CMAKE_CXX_STANDARD 14)

add_executable(PRT application.cpp
        ArcBall/Ball.h ArcBall/Ball.cpp ArcBall/BallAux.h ArcBall/BallAux.cpp ArcBall/BallMath.h ArcBall/BallMath.cpp
        OpenGL.h OpenGL.cpp
        OpenGLGeometry.h OpenGLGeometry.cpp
        Shaders.h Shaders.cpp
        Atlas.h Atlas.cpp
        Configuration.h
        Object3D.h Object3D.cpp
        Transformations.h Transformations.cpp
		Light.h Light.cpp
		stb_image.h stb_image.cpp
		PRT.h PRT.cpp PRTObject3D.cpp PRTObject3D.h)

target_link_libraries(PRT
        "-framework OpenGL"
        "armadillo"
        "freetype"
        "glfw")

target_include_directories(PRT PUBLIC "/usr/local/include/"
        "/usr/local/include/freetype2/")
