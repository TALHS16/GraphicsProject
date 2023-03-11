#include "Transform.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <glm/gtx/string_cast.hpp>

Transform::Transform() {
	init_trans();
}
void Transform::init_trans()
{
	this->transform_matrix = glm::mat4(1.0f);
	scale_factor = 1.f;
	move_factor = { 0,0,0 };
	rotate_factor = { 0,0,0 };
	rotate_world_factor = { 0,0,0 };
	move_world_factor = { 0,0,0 };
	world_transform_matrix = glm::mat4(1.0f);
	final_transform_matrix = glm::mat4(1.0f);
}
void Transform::SetMove_factor(glm::vec3 move_factor)
{
	move_factor = move_factor;
}

void Transform::SetScale_factor(float scale_factor)
{
	scale_factor = scale_factor;
}

void Transform::SetRotate_factor(glm::vec3 rotate_factor)
{
	rotate_factor = rotate_factor;
}

void Transform::SetRotate_world_factor(glm::vec3 rotate_world_factor)
{
	rotate_world_factor = rotate_world_factor;
}

float& Transform::get_scale_factor()
{
	return this->scale_factor;
}
glm::vec3& Transform::get_move_factor()
{
	return this->move_factor;
}
glm::vec3& Transform::get_world_move_factor()
{
	return this->move_world_factor;
}
glm::vec3& Transform::get_rotate_factor()
{
	return this->rotate_factor;
}

glm::vec3& Transform::get_rotate_world_factor()
{
	return this->rotate_world_factor;
}

void Transform::refresh_transform_matrix()
{
	this->transform_matrix = glm::mat4(1.0f);
	//scale
	glm::mat4 scale_matrix = glm::scale(glm::mat4(1.0f), glm::vec3(scale_factor, scale_factor, scale_factor));

	this->transform_matrix = scale_matrix * this->transform_matrix;


	//move
	glm::mat4 move_matrix = glm::mat4(1.0f);
	move_matrix[0][3] = this->move_factor[0];
	move_matrix[1][3] = this->move_factor[1];
	move_matrix[2][3] = this->move_factor[2];
	move_matrix = glm::transpose(move_matrix);
	glm::mat4 invert_move_matrix = glm::inverse(move_matrix);
	this->transform_matrix = move_matrix * this->transform_matrix;

	this->rotate_local_matrix = glm::mat4(1.0f);

	//rotate by shear
	if (0 < this->rotate_factor[0] && this->rotate_factor[0] < 360) {
		float radian_rotate = (360 - this->rotate_factor[0]) * (M_PI / 180);
		glm::mat4 rotate_matrix = glm::mat4(1.0f);
		float a;
		float b = sin(radian_rotate);
		float c;
		if (0 < this->rotate_factor[0] && this->rotate_factor[0] < 180) {
			a = -tan(radian_rotate / 2);
			c = a;
			rotate_matrix[1][1] = 1 + a * b;
			rotate_matrix[1][2] = a + a * b * c + c;
			rotate_matrix[2][1] = b;
			rotate_matrix[2][2] = b * c + 1;
		}
		if (180 <= this->rotate_factor[0] && this->rotate_factor[0] < 360) {
			a = tan((radian_rotate - M_PI) / 2);
			c = a;
			rotate_matrix[1][1] = -1 - a * b;
			rotate_matrix[1][2] = a + a * b * c + c;
			rotate_matrix[2][1] = b;
			rotate_matrix[2][2] = -b * c - 1;
		}
		rotate_matrix = glm::transpose(rotate_matrix);
		this->transform_matrix = move_matrix * rotate_matrix * invert_move_matrix * this->transform_matrix;
		this->rotate_local_matrix = rotate_matrix * this->rotate_local_matrix;
	}
	//rotate y
	if (0 < this->rotate_factor[1] && this->rotate_factor[1] < 360) {

		float radian_rotate_y = (360 - this->rotate_factor[1]) * (M_PI / 180);
		glm::mat4 rotate_y_matrix = glm::mat4(1.0f);

		if (0 < this->rotate_factor[1] && this->rotate_factor[1] < 180) {
			rotate_y_matrix[0][0] = cos(radian_rotate_y);
			rotate_y_matrix[2][0] = -sin(radian_rotate_y);
			rotate_y_matrix[0][2] = sin(radian_rotate_y);
			rotate_y_matrix[2][2] = cos(radian_rotate_y);
		}
		if (180 <= this->rotate_factor[1] && this->rotate_factor[1] < 360) {
			radian_rotate_y = radian_rotate_y - M_PI;
			rotate_y_matrix[0][0] = -cos(radian_rotate_y);
			rotate_y_matrix[2][0] = sin(radian_rotate_y);
			rotate_y_matrix[0][2] = -sin(radian_rotate_y);
			rotate_y_matrix[2][2] = -cos(radian_rotate_y);
		}
		rotate_y_matrix = glm::transpose(rotate_y_matrix);
		this->transform_matrix = move_matrix * rotate_y_matrix * invert_move_matrix * this->transform_matrix;
		this->rotate_local_matrix = rotate_y_matrix * this->rotate_local_matrix;

	}
	//rotate z
	if (0 < this->rotate_factor[2] && this->rotate_factor[2] < 360) {

		float radian_rotate_z = (360 - this->rotate_factor[2]) * (M_PI / 180);
		glm::mat4 rotate_z_matrix = glm::mat4(1.0f);

		if (0 < this->rotate_factor[2] && this->rotate_factor[2] < 180) {
			rotate_z_matrix[0][0] = cos(radian_rotate_z);
			rotate_z_matrix[0][1] = -sin(radian_rotate_z);
			rotate_z_matrix[1][0] = sin(radian_rotate_z);
			rotate_z_matrix[1][1] = cos(radian_rotate_z);
		}
		if (180 <= this->rotate_factor[2] && this->rotate_factor[2] < 360) {
			radian_rotate_z = radian_rotate_z - M_PI;
			rotate_z_matrix[0][0] = -cos(radian_rotate_z);
			rotate_z_matrix[0][1] = sin(radian_rotate_z);
			rotate_z_matrix[1][0] = -sin(radian_rotate_z);
			rotate_z_matrix[1][1] = -cos(radian_rotate_z);
		}
		rotate_z_matrix = glm::transpose(rotate_z_matrix);
		this->transform_matrix = move_matrix * rotate_z_matrix * invert_move_matrix * this->transform_matrix;

		this->rotate_local_matrix = rotate_z_matrix * this->rotate_local_matrix;

	}


	this->rotate_word_matrix = glm::mat4(1.0f);



	//world rotate X
	if (0 < this->rotate_world_factor[0] && this->rotate_world_factor[0] < 360) {
		float radian_rotate_world = (360 - this->rotate_world_factor[0]) * (M_PI / 180);
		glm::mat4 rotate_world_matrix = glm::mat4(1.0f);

		if (0 < this->rotate_world_factor[0] && this->rotate_world_factor[0] < 180) {
			rotate_world_matrix[1][1] = cos(radian_rotate_world);
			rotate_world_matrix[1][2] = -sin(radian_rotate_world);
			rotate_world_matrix[2][1] = sin(radian_rotate_world);
			rotate_world_matrix[2][2] = cos(radian_rotate_world);
		}
		if (180 <= this->rotate_world_factor[0] && this->rotate_world_factor[0] < 360) {
			radian_rotate_world = radian_rotate_world - M_PI;
			rotate_world_matrix[1][1] = -cos(radian_rotate_world);
			rotate_world_matrix[1][2] = sin(radian_rotate_world);
			rotate_world_matrix[2][1] = -sin(radian_rotate_world);
			rotate_world_matrix[2][2] = -cos(radian_rotate_world);
		}
		rotate_world_matrix = glm::transpose(rotate_world_matrix);
		this->rotate_word_matrix = rotate_world_matrix * this->rotate_word_matrix;

	}

	//world rotate y
	if (0 < this->rotate_world_factor[1] && this->rotate_world_factor[1] < 360) {
		float radian_rotate_world = (360 - this->rotate_world_factor[1]) * (M_PI / 180);
		glm::mat4 rotate_world_matrix = glm::mat4(1.0f);

		if (0 < this->rotate_world_factor[1] && this->rotate_world_factor[1] < 180) {
			rotate_world_matrix[0][0] = cos(radian_rotate_world);
			rotate_world_matrix[2][0] = -sin(radian_rotate_world);
			rotate_world_matrix[0][2] = sin(radian_rotate_world);
			rotate_world_matrix[2][2] = cos(radian_rotate_world);
		}
		if (180 <= this->rotate_world_factor[1] && this->rotate_world_factor[1] < 360) {
			radian_rotate_world = radian_rotate_world - M_PI;
			rotate_world_matrix[0][0] = -cos(radian_rotate_world);
			rotate_world_matrix[2][0] = sin(radian_rotate_world);
			rotate_world_matrix[0][2] = -sin(radian_rotate_world);
			rotate_world_matrix[2][2] = -cos(radian_rotate_world);
		}
		rotate_world_matrix = glm::transpose(rotate_world_matrix);
		this->rotate_word_matrix = rotate_world_matrix * this->rotate_word_matrix;

	}

	//world rotate z
	if (0 < this->rotate_world_factor[2] && this->rotate_world_factor[2] < 360) {
		float radian_rotate_world = (360 - this->rotate_world_factor[2]) * (M_PI / 180);
		glm::mat4 rotate_world_matrix = glm::mat4(1.0f);

		if (0 < this->rotate_world_factor[2] && this->rotate_world_factor[2] < 180) {
			rotate_world_matrix[0][0] = cos(radian_rotate_world);
			rotate_world_matrix[0][1] = -sin(radian_rotate_world);
			rotate_world_matrix[1][0] = sin(radian_rotate_world);
			rotate_world_matrix[1][1] = cos(radian_rotate_world);
		}
		if (180 <= this->rotate_world_factor[2] && this->rotate_world_factor[2] < 360) {
			radian_rotate_world = radian_rotate_world - M_PI;
			rotate_world_matrix[0][0] = -cos(radian_rotate_world);
			rotate_world_matrix[0][1] = sin(radian_rotate_world);
			rotate_world_matrix[1][0] = -sin(radian_rotate_world);
			rotate_world_matrix[1][1] = -cos(radian_rotate_world);
		}
		rotate_world_matrix = glm::transpose(rotate_world_matrix);
		this->rotate_word_matrix = rotate_world_matrix * this->rotate_word_matrix;

	}

	this->world_transform_matrix = this->rotate_word_matrix;

	//move
	move_matrix = glm::mat4(1.0f);
	move_matrix[0][3] = this->move_world_factor[0];
	move_matrix[1][3] = this->move_world_factor[1];
	move_matrix[2][3] = this->move_world_factor[2];
	move_matrix = glm::transpose(move_matrix);
	invert_move_matrix = glm::inverse(move_matrix);
	this->world_transform_matrix = move_matrix * this->world_transform_matrix;
}

void Transform::set_final_transform_matrix(glm::mat4 final)
{
	final_transform_matrix = final;
}

glm::mat4& Transform::get_final_transform_matrix()
{
	return final_transform_matrix;
}

glm::mat4& Transform::get_local_transform_matrix()
{
	return transform_matrix;
}

glm::mat4& Transform::get_world_transform_matrix()
{
	return world_transform_matrix;
}
