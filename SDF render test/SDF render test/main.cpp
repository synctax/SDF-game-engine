#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <iostream>

//GLM
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// Other includes
#include "Shader.h"
#include "Camera.h"
#include "sdfTree.h"
#include "SOIL2/SOIL2.h"

#define PI (3.14159265)

using namespace glm;

bool keys[256];
Camera playerCam(vec3(0,0,0), vec3(0,0,1));
float playerSpeed = 10;
float playerLookSpeed = 1.5;
double lastTime = glfwGetTime();
int nbFrames = 0;
float deltaT = 16.66;



// Window dimensions
//const GLuint WIDTH = 500, HEIGHT = 500;
const GLuint WIDTH = 320, HEIGHT = 800;
//const GLuint WIDTH = 1280, HEIGHT = 800;

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods){
    if (action == GLFW_PRESS){
        keys[key] = true;
    }else if (action == GLFW_RELEASE){
        keys[key] = false;
    }
}

void handleKeys(Shader* shader){
    mat4 rotationMat;
    vec3 p;
    float s = playerSpeed*(deltaT/1000);
    float ts = playerLookSpeed*(deltaT/1000);
    unsigned char key;
    for (int i = 0; i < 256; i++){
        if (keys[i]){
            key = i;
            switch((char)key+32){
                case 'f':
                    shader->update();
                    break;
                    
                case 'w':
                    playerCam.moveForward(s);
                    break;
                    
                case 's':
                    playerCam.moveForward(-s);
                    break;
                    
                case 'a':
                    playerCam.moveRight(-s);
                    break;
                    
                case 'd':
                    playerCam.moveRight(s);
                    break;
                    
                case 'q':
                    playerCam.moveUp(s);
                    break;
                    
                case 'e':
                    playerCam.moveUp(-s);
                    break;
                    
                case 'j':
                    playerCam.lookRight(-ts);
                    break;
                    
                case 'l':
                    playerCam.lookRight(ts);
                    break;
                    
                case 'i':
                    playerCam.lookUp(ts);
                    break;
                    
                case 'k':
                    playerCam.lookUp(-ts);
                    break;
            }
        }
    }
}

struct sdfNode{
    int data[5];
};

void updateDeltaT(GLFWwindow* window){
    double currentTime = glfwGetTime();
    nbFrames++;
    // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
    
    if ( currentTime - lastTime >= 1.0 ){ // If last prinf() was more than 1 sec ago
        // printf and reset timer
        char buffer[64];
        deltaT = 1000.0/double(nbFrames);
        snprintf(buffer, sizeof buffer, "%f",deltaT);
        glfwSetWindowTitle(window, buffer);
        nbFrames = 0;
        lastTime += 1.0;
    }
}



// The MAIN function, from here we start the application and run the game loop
int main( )
{
    // Init GLFW
    glfwInit( );
    
    // Set all the required options for GLFW
    glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 3 );
    glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 3 );
    glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
    glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE );
    
    // Create a GLFWwindow object that we can use for GLFW's functions
    GLFWwindow *window = glfwCreateWindow( WIDTH, HEIGHT, "LearnOpenGL", nullptr, nullptr );
    
    int screenWidth, screenHeight;
    glfwGetFramebufferSize( window, &screenWidth, &screenHeight );
    
    if ( nullptr == window )
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate( );
        
        return EXIT_FAILURE;
    }
    
    glfwMakeContextCurrent( window );
    
    // Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    glewExperimental = GL_TRUE;
    // Initialize GLEW to setup the OpenGL Function pointers
    if ( GLEW_OK != glewInit( ) )
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Define the viewport dimensions
    glViewport( 0, 0, screenWidth, screenHeight );
    
    // Build and compile our shader program
    Shader ourShader( "/Users/benw/Documents/Projects/SDF game Engine/SDF-game-engine/SDF render test/SDF render test/core.vs", "/Users/benw/Documents/Projects/SDF game Engine/SDF-game-engine/SDF render test/SDF render test/core.frag" );
    
    // Set up vertex data (and buffer(s)) and attribute pointers
    GLfloat vertices[] =
    {
        // Positions          // Colors           // Texture Coords
        1.0f,  1.0f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f, // Top Right
        1.0f, -1.0f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f, // Bottom Right
        -1.0f, -1.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f, // Bottom Left
        -1.0f,  1.0f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f  // Top Left
    };
    GLuint indices[] =
    {  // Note that we start from 0!
        0, 1, 3, // First Triangle
        1, 2, 3  // Second Triangle
    };
    GLuint VBO, VAO, EBO;
    glGenVertexArrays( 1, &VAO );
    glGenBuffers( 1, &VBO );
    glGenBuffers( 1, &EBO );
    
    glBindVertexArray( VAO );
    
    glBindBuffer( GL_ARRAY_BUFFER, VBO );
    glBufferData( GL_ARRAY_BUFFER, sizeof( vertices ), vertices, GL_STATIC_DRAW );
    
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, EBO );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof( indices ), indices, GL_STATIC_DRAW );
    
    // Position attribute
    glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof( GLfloat ), ( GLvoid * ) 0 );
    glEnableVertexAttribArray(0);
    // Color attribute
    glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof( GLfloat ), ( GLvoid * )( 3 * sizeof( GLfloat ) ) );
    glEnableVertexAttribArray(1);
    // Texture Coordinate attribute
    glVertexAttribPointer( 2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof( GLfloat ), ( GLvoid * )( 6 * sizeof( GLfloat ) ) );
    glEnableVertexAttribArray( 2 );
    
    glBindVertexArray( 0 ); // Unbind VAO
    
    
    glfwSetKeyCallback(window, keyCallback);
    
    // Game loop
    while ( !glfwWindowShouldClose( window ) )
    {
        updateDeltaT(window);
        glfwPollEvents();
        handleKeys(&ourShader);
        
        // Render
        // Clear the colorbuffer
        glClearColor( 0.2f, 0.3f, 0.3f, 1.0f );
        glClear( GL_COLOR_BUFFER_BIT );
        
        // Draw the triangle
        ourShader.Use( );
        
        // Draw container
        glBindVertexArray( VAO );
        glDrawElements( GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0 );
        glBindVertexArray( 0 );
        
        
        GLint posLoc = glGetUniformLocation(ourShader.Program, "playerPosition");
        glProgramUniform3f(ourShader.Program, posLoc, playerCam.position.x, playerCam.position.y, playerCam.position.z);
        GLint lookLoc = glGetUniformLocation(ourShader.Program, "playerView");
        glProgramUniform3f(ourShader.Program, lookLoc, playerCam.direction.x, playerCam.direction.y, playerCam.direction.z);
        
        // Swap the screen buffers
        glfwSwapBuffers( window );
        
    }
    
    // Properly de-allocate all resources once they've outlived their purpose
    glDeleteVertexArrays( 1, &VAO );
    glDeleteBuffers( 1, &VBO );
    glDeleteBuffers( 1, &EBO );
    
    // Terminate GLFW, clearing any resources allocated by GLFW.
    glfwTerminate( );
    
    return EXIT_SUCCESS;
}







