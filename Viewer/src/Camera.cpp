#include "Camera.h"
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/glm.hpp>


float Camera::screen_hight = 0;
float Camera::screen_width = 0;

Camera::Camera()
{
	init();
	

}
void Camera::init()
{
	init_trans();
	flag_camera_view = 0;
	view_transformation = glm::mat4(1.f);
	projection_transformation = glm::mat4(1.f);
	eye = glm::vec3(0.0f, 0.0f, 1000.0f);
	at = glm::vec3(0.0f, 0.0f, 0.0f);
	up = glm::vec3(0.0f, 1.0f, 0.0f);
	this->ortho_vec = { -1280 / 2, 1280 / 2, -720 / 2, 720 / 2, -1, 1 };
	this->pers_vec = { 45.f,1280.f / 720.f, 1.f ,1000.f };
	dolly_zoom = 0;
}

Camera::~Camera()
{
	
}

void Camera::SetProjection()
{
	if (flag_camera_view)
		projection_transformation = glm::perspective(glm::radians(this->pers_vec[0]), this->pers_vec[1], this->pers_vec[2], this->pers_vec[3]);//glm::frustum(this->pers_vec[0], this->pers_vec[1], this->pers_vec[2], this->pers_vec[3], this->pers_vec[4], this->pers_vec[5]);
	else
		projection_transformation = glm::ortho(this->ortho_vec[0], this->ortho_vec[1], this->ortho_vec[2], this->ortho_vec[3], this->ortho_vec[4], this->ortho_vec[5]);
	
}

void Camera::set_camera_view(bool current_item)
{
	flag_camera_view = current_item;
}

int Camera::get_camera_view()
{
	return flag_camera_view;
}

const glm::mat4x4& Camera::GetProjectionTransformation() 
{
	SetProjection();
	return projection_transformation;
}

glm::mat4 Camera::final_matrix()
{
	return GetProjectionTransformation() * GetViewTransformation() * glm::inverse(this->get_world_transform_matrix() * this->get_local_transform_matrix());
}

void Camera::SetCameraLookAt()
{
	view_transformation = glm::lookAt(eye, at, up);
}

const glm::mat4x4& Camera::GetViewTransformation()
{
	SetCameraLookAt();
	return view_transformation;
}

glm::vec3& Camera::get_eye()
{
	return eye;
}

glm::vec3& Camera::get_at()
{
	return at;
}

glm::vec3& Camera::get_up()
{
	return up;
}




void Camera::refresh()
{
	this->ortho_vec[0] = -screen_width / 2;
	this->ortho_vec[1] = screen_width / 2;
	this->ortho_vec[2] = -screen_hight / 2;
	this->ortho_vec[3] = screen_hight / 2;
}

std::vector<float>& Camera::get_ortho()
{
	return this->ortho_vec;
}

std::vector<float>& Camera::get_pers()
{
	return this->pers_vec;
}

float Camera::get_dolly_zoom() const
{
	return dolly_zoom;
}
void Camera::setDollyZoom(float x)
{
	dolly_zoom = x;
	pers_vec[0] = 45 + 0.5 * x;
	eye[2] = 1000 - 10 * x;
}