#include "MeshModel.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <glm/gtx/string_cast.hpp>
#include "Renderer.h"

MeshModel::MeshModel(std::vector<Face> faces, std::vector<glm::vec3> vertices, std::vector<glm::vec3> normals, const std::string& model_name) :
	faces(faces),
	vertices(vertices),
	normals(normals)
{
	NormolizeModel();
	this->transform_matrix = glm::mat4(1.0f);
  this->activeCamera=NULL;
    show_axes = 0;
    show_faces_normals = 0;
    show_verteces_normals = 0;
    show_object = 1;
    show_bound_Triangle = 0;
    fill_triangle = 0;
    this->up=-1;
    this->right=-1;
    this->left = Camera::screen_hight*Camera::screen_width;
    this->down = Camera::screen_hight*Camera::screen_width;
    this->zMax = -1;
    this->zMin = Camera::screen_hight*Camera::screen_width;
    cout << faces.size() << endl;
    for (int i = 0; i < faces.size(); i++)
    {
        model_fill_color.push_back({ (float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX });
    }
    
}

void MeshModel::Print_Verrices() const
{
	std::cout<< "MeshModel Name: "<<this->GetModelName()<<std::endl;
	std::cout<< "Vertices - "<<this->vertices.size()<<"-"<<std::endl;
	for (glm::vec3 v : this->vertices){
		printf("%.5f\t%.5f\t%.5f\n", v[0], v[1],v[2]);
	}
		
}


void MeshModel::Print_Faces() const
{
	std::cout<< "Faces - "<<this->faces.size()<<"-"<<std::endl;
	for( int i=0;i<this->faces.size();i++){
		printf("%d: { ", i);
		for(int j=0;j<3;j++){
				Face f = this->GetFace(i);
				printf("[v: %d  vt: %d  vn: %d]", f.GetVertexIndex(j), f.GetTextureIndex(j), f.GetNormalIndex(j));
			}
		printf(" } \n");
	}
	
		
}

glm::vec3 MeshModel::get_vertex3(int index){
    return this->vertices[index];
}

std::vector<glm::vec3> MeshModel::get_model_fill_color()
{
    return model_fill_color;
}


glm::vec2 MeshModel::get_normal_of_face(int face_index){
    glm::vec4 vertex4 = glm::vec4(this->faces_normals[face_index], 1.0f);
    vertex4 = this->activeCamera->final_matrix() * this->rotate_word_matrix * this->rotate_local_matrix * vertex4;
    
    vertex4/=vertex4[3];
    
    return  glm::vec2(vertex4[0], vertex4[1]);
}

glm::vec3& MeshModel::get_model_color()
{
	return this->color_model;
}

void MeshModel::set_model_name(const std::string& model_name2)
{
	this->model_name = model_name2;
}

void MeshModel::set_model_color(const glm::vec3& color_model)
{
	this->color_model = color_model;
}


bool& MeshModel::get_show_axes(){
    return this->show_axes;
}
bool& MeshModel::get_show_faces_normals(){
    return this->show_faces_normals;
}
bool& MeshModel::get_show_object(){
    return this->show_object;
}
bool& MeshModel::get_show_verteces_normals(){
    return this->show_verteces_normals;
}

bool& MeshModel::get_show_bound_triangle(){
    return this->show_bound_Triangle;
}

bool& MeshModel::get_fill_triangle()
{
    return this->fill_triangle;
}



glm::vec2 MeshModel::get_vertex2(int index, int first_cordinate, int second_cordinate)
{
    glm::vec4 vertex4_local = this->get_vertex4(index);
    vertex4_local = this->activeCamera->final_matrix() *  this->get_local_transform_matrix() * vertex4_local;
    
    vertex4_local/=vertex4_local[3];
    
    if(vertex4_local[1]>up)
        this->up=vertex4_local[1];
    if(vertex4_local[1]<down)
        this->down=vertex4_local[1];
    if(vertex4_local[0]>right)
        this->right=vertex4_local[0];
    if(vertex4_local[0]<left)
        this->left=vertex4_local[0];
    
    glm::vec4 vertex4 = this->get_vertex4(index);
    vertex4 = this->activeCamera->final_matrix() * this->get_world_transform_matrix() * this->get_local_transform_matrix() * vertex4;
    
    vertex4/=vertex4[3];
    
    return  glm::vec2(vertex4[first_cordinate], vertex4[second_cordinate]);
}


glm::vec2 MeshModel::get_normal2(int index, int first_cordinate, int second_cordinate)
{
    glm::vec4 vertex4 = this->get_normal4(index);
    vertex4 = this->activeCamera->final_matrix() * this->rotate_word_matrix * this->rotate_local_matrix * vertex4;
    
    vertex4/=vertex4[3];
    
    return  glm::vec2(vertex4[first_cordinate], vertex4[second_cordinate]);
}

const Face& MeshModel::GetFace(int index) const
{
	return faces[index];
}

int MeshModel::GetFacesCount() const
{
	return faces.size();
}
int MeshModel::GetNormalsCount() const
{
	return normals.size();
}

const std::string& MeshModel::GetModelName() const
{
	return model_name;
}

glm::vec3 MeshModel::get_vertex3_after_transform(int index)
{
    glm::vec4 vertex4_local = this->get_vertex4(index);
    vertex4_local = this->activeCamera->final_matrix() *  this->get_local_transform_matrix() * vertex4_local;
    
    vertex4_local/=vertex4_local[3];
    
    if(vertex4_local[1]>up)
        this->up=vertex4_local[1];
    if(vertex4_local[1]<down)
        this->down=vertex4_local[1];
    if(vertex4_local[0]>right)
        this->right=vertex4_local[0];
    if(vertex4_local[0]<left)
        this->left=vertex4_local[0];
    if(vertex4_local[2]>zMax)
        this->zMax=vertex4_local[2];
    if(vertex4_local[2]<zMin)
        this->zMin=vertex4_local[2];
    
    glm::vec4 vertex4 = this->get_vertex4(index);
    vertex4 = this->activeCamera->final_matrix() * this->get_world_transform_matrix() * this->get_local_transform_matrix() * vertex4;
    
    vertex4/=vertex4[3];
    return vertex4;
}


void MeshModel::refresh_transform_matrix(){
    
    
    Transform::refresh_transform_matrix();
    
}

void MeshModel::set_final_transform_matrix(glm::mat4 final)
{
	final_transform_matrix = final;
}

glm::mat4& MeshModel::get_final_transform_matrix()
{
	return final_transform_matrix;
}

void MeshModel::scale(float scale_factor)
{
	for (int i = 0; i < this->vertices.size(); i++) {
		this->vertices[i] *= scale_factor;
	}
}

MeshModel::~MeshModel()
{
}

void MeshModel::NormolizeModel() {
	float minX, maxX, minY, maxY, maxZ, minZ;
	minX = vertices[0].x;
	maxX = vertices[0].x;
	minY = vertices[0].y;
	maxY = vertices[0].y;
	minZ = vertices[0].z;
	maxZ = vertices[0].z;
	for (int i = 1; i < vertices.size(); i++) {
		minX = std::min(minX, vertices[i].x);
		minY = std::min(minY, vertices[i].y);
		minZ = std::min(minZ, vertices[i].z);
		maxX = std::max(maxX, vertices[i].x);
		maxY = std::max(maxY, vertices[i].y);
		maxZ = std::max(maxZ, vertices[i].z);
	}
	float size = std::max(maxX - minX, maxY - minY);
	size = std::max(size, maxZ - minZ);
	float scale_factor = 10.f / size;
	glm::mat4 scaleMat = glm::scale(glm::mat4(1.0f), glm::vec3(scale_factor, scale_factor, scale_factor));
	for (int i = 0; i < vertices.size(); i++) {
		glm::vec4 v(vertices[i], 1.f);
		v = scaleMat * v;
		vertices[i] = v;
	}
}

glm::vec2 MeshModel::world_transform_point(glm::vec3 point){
    glm::vec4 point4 = glm::vec4(point, 1.0f);
    point4 = this->activeCamera->final_matrix() * this->world_transform_matrix * point4;
    
    point4/=point4[3];
    
    return glm::vec2(point4[0], point4[1]);
}

glm::vec2 MeshModel::local_transform_point(glm::vec3 point){
    glm::vec4 point4 = glm::vec4(point, 1.0f);
    point4 = this->activeCamera->final_matrix() * this->transform_matrix * point4;
    
    point4/=point4[3];
    
    return glm::vec2(point4[0], point4[1]);
}

glm::vec2 MeshModel::full_tranform(glm::vec3 point){
    glm::vec3 point3 = this->full_tranform3(point);
    
    return glm::vec2(point3[0], point3[1]);
}

glm::vec3 MeshModel::full_tranform3(glm::vec3 point){
    glm::vec4 point4 = glm::vec4(point, 1.0f);
    point4 = this->activeCamera->final_matrix() * this->world_transform_matrix * this->transform_matrix * point4;
    
    point4/=point4[3];
    
    return glm::vec3(point4[0], point4[1], point4[2]);
}


glm::mat4& MeshModel::get_local_transform_matrix(){
    return this->transform_matrix;
}
glm::mat4& MeshModel::get_world_transform_matrix(){
    return this->world_transform_matrix;
}


