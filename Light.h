#pragma once
#include <glm/glm.hpp>
#include <vector>

class Light:
{
public:
    Light(glm::vec3 position);
    glm::vec3& getPosition();
    glm::vec3& getAmbient();
    glm::vec3& getDiffuse();
    glm::vec3& getSpecular();
    
private:
    glm::vec3 position;
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};
