#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Utils.h"

glm::vec3 Utils::Vec3fFromStream(std::istream& issLine)
{
	float x, y, z;
	issLine >> x >> std::ws >> y >> std::ws >> z;
	return glm::vec3(x, y, z);
}

glm::vec2 Utils::Vec2fFromStream(std::istream& issLine)
{
	float x, y;
	issLine >> x >> std::ws >> y;
	return glm::vec2(x, y);
}

std::shared_ptr<MeshModel> Utils::LoadMeshModel(const std::string& filePath)
{
	std::vector<Face> faces;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> textureCoords;
	std::ifstream ifile(filePath.c_str());

	// while not end of file
	while (!ifile.eof())
	{
		// get line
		std::string curLine;
		std::getline(ifile, curLine);

		// read the type of the line
		std::istringstream issLine(curLine);
		std::string lineType;

		issLine >> std::ws >> lineType;

		// based on the type parse data
		if (lineType == "v")
		{
			vertices.push_back(Utils::Vec3fFromStream(issLine));
		}
		else if (lineType == "vn")
		{
			normals.push_back(Utils::Vec3fFromStream(issLine));
		}
		else if (lineType == "vt")
		{
			textureCoords.push_back(Utils::Vec2fFromStream(issLine));
		}
		else if (lineType == "f")
		{
			faces.push_back(Face(issLine));
		}
		else if (lineType == "#" || lineType == "")
		{
			// comment / empty line
		}
		else
		{
			std::cout << "Found unknown line Type \"" << lineType << "\"";
		}
	}
	std::shared_ptr<MeshModel> new_meshModel = std::make_shared<MeshModel>(faces, vertices, CalNormals(vertices, faces), Utils::GetFileName(filePath), textureCoords);
	new_meshModel->set_model_name(GetFileName(filePath));
	return new_meshModel;
}

std::vector<glm::vec3> Utils::CalNormals(std::vector<glm::vec3> vertices, std::vector<Face> faces)
{
	std::vector<int> faces_count(vertices.size());
	std::vector<glm::vec3> normals(vertices.size());
	for (int i = 0; i < faces_count.size(); i++)
	{
		faces_count[i] = 0;
	}
	for (int i = 0; i < faces.size(); i++) {
		int point0_index = faces.at(i).GetVertexIndex(0) - 1;
		int point1_index = faces.at(i).GetVertexIndex(1) - 1;
		int point2_index = faces.at(i).GetVertexIndex(2) - 1;

		glm::vec3 point0 = vertices.at(point0_index);
		glm::vec3 point1 = vertices.at(point1_index);
		glm::vec3 point2 = vertices.at(point2_index);

		glm::vec3 normal = glm::cross( point0 - point1, point2 - point1);

		normal = glm::normalize(-normal);

		normals.at(point0_index) += normal;
		normals.at(point1_index) += normal;
		normals.at(point2_index) += normal;
		
		faces_count.at(point0_index) += 1;
		faces_count.at(point1_index) += 1;
		faces_count.at(point2_index) += 1;
	}

	for (int i = 0; i < normals.size(); i++)
	{
		normals[i] /= faces_count[i];
		normals[i] = glm::normalize(normals[i]);
	}

	return normals;
}

std::string Utils::GetFileName(const std::string& filePath)
{
	if (filePath.empty()) {
		return {};
	}

	auto len = filePath.length();
	auto index = filePath.find_last_of("/\\");

	if (index == std::string::npos) {
		return filePath;
	}

	if (index + 1 >= len) {

		len--;
		index = filePath.substr(0, len).find_last_of("/\\");

		if (len == 0) {
			return filePath;
		}

		if (index == 0) {
			return filePath.substr(1, len - 1);
		}

		if (index == std::string::npos) {
			return filePath.substr(0, len);
		}

		return filePath.substr(index + 1, len - index - 1);
	}

	return filePath.substr(index + 1, len - index);
}
