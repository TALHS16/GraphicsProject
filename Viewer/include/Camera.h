#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "Transform.h"

class Camera: public Transform
{
public:
	static float screen_width;
	static float screen_hight;
	Camera();
	virtual ~Camera();
	void init();
	void SetCameraLookAt();
	void SetProjection();
	const glm::mat4x4& GetProjectionTransformation();
	const glm::mat4x4& GetViewTransformation();
	std::vector<float>& get_ortho();
	std::vector<float>& get_pers();
	void set_camera_view(bool current_item);
	int get_camera_view();
	void refresh();
	glm::vec3& get_eye();
	glm::vec3& get_at();
	glm::vec3& get_up();
	float get_dolly_zoom() const;
	glm::mat4 final_matrix(); 
	void setDollyZoom(float x);
private:
	glm::vec3 eye;
	glm::vec3 at;
	glm::vec3 up;
	glm::mat4x4 view_transformation;
	glm::mat4x4 projection_transformation;
	std::vector<float> ortho_vec;
	std::vector<float> pers_vec;
	float dolly_zoom;
	bool flag_camera_view;
};
