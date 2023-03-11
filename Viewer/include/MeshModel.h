#pragma once
#include <glm/glm.hpp>
#include <string>
#include "Face.h"
#include "Transform.h"
#include "Camera.h"

class MeshModel:public Transform
{
public:
	MeshModel(std::vector<Face> faces, std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, const std::string& model_name);
	virtual ~MeshModel();
	void move(glm::vec3 moveCordinate);
	void scale(float scale_factor);
	const Face& GetFace(int index) const;
	int GetFacesCount() const;
	int GetNormalsCount() const;
	const std::string& GetModelName() const;
	void Print_Verrices() const;
	void Print_Faces() const;
	glm::vec2 get_vertex2(int index, int first_cordinate, int second_cordinate);
    glm::vec2 get_normal2(int index, int first_cordinate, int second_cordinate);
	glm::vec3 get_vertex3_after_transform(int index);
	glm::vec3 get_normal3_after_transform(int index)const;
	void set_model_name(const std::string& model_name2);
	glm::vec3& get_model_color();
	void set_model_color(const glm::vec3& model_color);
	void NormolizeModel();
    bool& get_show_axes();
    bool& get_show_faces_normals();
    bool& get_show_verteces_normals();
    bool& get_show_object();
    bool& get_show_bound_triangle();
    bool& get_fill_triangle();
    glm::mat4& get_final_transform_matrix();
	glm::mat4& get_local_transform_matrix();
	glm::mat4& get_world_transform_matrix();
    void set_final_transform_matrix(glm::mat4 final);
    Camera* activeCamera;
    glm::vec2 world_transform_point(glm::vec3 point);
    glm::vec2 local_transform_point(glm::vec3 point);
    glm::mat4 inv_camera_final_mat;
    float up, down, left, right, zMax, zMin;
    glm::vec2 get_normal_of_face(int face_index);
    glm::vec2 full_tranform(glm::vec3 point);
    glm::vec3 full_tranform3(glm::vec3 point);
    void refresh_transform_matrix();
	std::vector<glm::vec3> faces_normals;
    glm::vec3 get_vertex3(int index);
    std::vector<glm::vec3> get_model_fill_color();
    

private:
    std::vector<Face> faces;
    std::vector<glm::vec3> model_fill_color;
    
    inline glm::vec4 get_vertex4(int index)const{ return glm::vec4(this->vertices[index], 1.0f);};
    inline glm::vec4 get_normal4(int index)const{ return glm::vec4(this->normals[index], 1.0f);};
    
    bool show_axes;
    bool show_faces_normals;
    bool show_verteces_normals;
    bool show_object;
    bool show_bound_Triangle;
    bool fill_triangle;
	std::vector<glm::vec3> vertices;
	
	
	std::vector<glm::vec3> normals;
	std::string model_name;
	glm::vec3 color_model;
};
