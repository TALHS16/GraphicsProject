#include "Material.h"


Material::Material(glm::vec3 ambient,glm::vec3 diffuse, glm::vec3 specular):ambient(ambient), diffuse(diffuse), specular(specular)
{
    
}
glm::vec3& Material::get_ambient(){
    return this->ambient;
}
glm::vec3& Material::get_diffuse(){
    return this->diffuse;
}
glm::vec3& Material::get_specular(){
    return this->specular;
}

float* Material::get_shine_factor()
{
    return &this->shine_factor;
}
