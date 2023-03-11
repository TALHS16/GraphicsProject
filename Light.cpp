#include "Light.h"


Light::Light(glm::vec3 position):position(position){}

glm::vec3& Light::getPosition(){
    return this->position;
}

glm::vec3& Light::getAmbient(){
    return this->ambient;
}

glm::vec3& Light::getDiffuse(){
    return this->diffuse;
}

glm::vec3& Light::getSpecular(){
    return this->specular;
}

