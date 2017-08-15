//
//  Camera.h
//  SDF renderer with shader
//
//  Created by Ben Wyngaard on 8/1/17.
//  Copyright © 2017 Ben Wyngaard. All rights reserved.
//

#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp>

#ifndef Camera_h
#define Camera_h

using namespace glm;

class Camera{
public:
    vec3 position, direction;
    Camera(vec3 p, vec3 d){
        this->position = p;
        this->direction = d;
    }
    
    void moveForward(float speed){
        this->position += this->direction*speed;
    }
    
    void moveRight(float speed){
        vec3 p = cross(this->direction, vec3(0,1,0));
        this->position -= p*speed;
    }
    
    void lookRight(float speed){
        vec3 p = normalize(cross(vec3(0,1,0),this->direction));
        this->direction += p*speed;
        normalize(this->direction);
    }
    
    void lookUp(float speed){
        vec3 p = normalize(cross(vec3(0,1,0),this->direction));
        p = normalize(cross(p, this->direction));
        this->direction -= p*speed;
        normalize(this->direction);
    }
    
    void moveUp(float speed){
        this->position.y += speed;
    }
};


#endif /* Camera_h */
