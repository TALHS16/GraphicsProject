#include "Light.h"
#include "Renderer.h"

Light::Light(glm::vec3 position){
    this->position = position;
}
glm::vec3& Light::get_position()
{
    return this->position;
}

glm::vec3 Light::get_norm_position()
{
    return Norm3Copy(this->position);
}
