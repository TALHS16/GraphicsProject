#pragma once
#include <glm/glm.hpp>
#include <vector>

class Material
{
public:
    Material(glm::vec3 ambient = {0.0f,0.0f,0.0f},glm::vec3 diffuse = {0.0f,0.0f,0.0f}, glm::vec3 specular = {0.0f,0.0f,0.0f});
    glm::vec3& get_ambient();
    glm::vec3& get_diffuse();
    glm::vec3& get_specular();
    float* get_shine_factor();
protected:
	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;
    float shine_factor;

};
