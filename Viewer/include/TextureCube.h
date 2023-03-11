#pragma once
#include <string>
#include <vector>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
class TextureCube
{
public:
	TextureCube();
	int loadCubemap();
private:
	std::vector<std::string> faces;

};

