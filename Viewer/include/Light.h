#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "Material.h"

class Light
{
public:
    Light(glm::vec3 position = {0.0f,0.0f,0.0f});
    glm::vec3& get_position();
    glm::vec3 get_norm_position();
    Material materials;
private:
	glm::vec3 position;
};
