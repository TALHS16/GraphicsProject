#pragma once
#include <glm/glm.hpp>
#include <vector>

class Transform
{
public:
	Transform();
	~Transform() {};
	void init_trans();
	void SetMove_factor(glm::vec3 move_factor);
	void SetScale_factor(float scale_factor);
	void SetRotate_factor(glm::vec3 rotate_factor);
	void SetRotate_world_factor(glm::vec3 rotate_wold_factor);
	float& get_scale_factor();
	glm::vec3& get_move_factor();
	glm::vec3& get_world_move_factor();
	glm::vec3& get_rotate_factor();
	glm::vec3& get_rotate_world_factor();
    virtual void refresh_transform_matrix();
	void set_final_transform_matrix(glm::mat4 final);
	glm::mat4& get_final_transform_matrix();
	glm::mat4& get_local_transform_matrix();
	glm::mat4& get_world_transform_matrix();


protected:
	float scale_factor;
	glm::vec3 move_factor;
	glm::vec3 rotate_factor;
	glm::vec3 rotate_world_factor;
	glm::vec3 move_world_factor;
	glm::mat4 transform_matrix;
	glm::mat4 world_transform_matrix;
	glm::mat4 final_transform_matrix;
	glm::mat4 rotate_local_matrix;
	glm::mat4 rotate_word_matrix;
};

