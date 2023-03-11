#pragma once

#include <vector>
#include <memory>

#include "Camera.h"
#include "MeshModel.h"
#include "Light.h"
#include "Light.h"
#include "TextureCube.h"


using namespace std;

class Scene {
public:
	Scene();

	void AddModel(const shared_ptr<MeshModel>& mesh_model);
    void RemoveModelByIndex(int index);
	int GetModelCount() const;
	MeshModel& GetModel(int index) const;
	MeshModel& GetActiveModel() const;
	
	void AddCamera(const shared_ptr<Camera>& camera);
	int GetCameraCount() const;
	Camera& GetCamera(int index);
	Camera& GetActiveCamera() const;

	void SetActiveCameraIndex(int index);
	int GetActiveCameraIndex() const;

	void SetActiveModelIndex(int index);
	int GetActiveModelIndex() const;
    
    void AddLight(const shared_ptr<Light>& camera);
    int GetLightCount() const;
    Light& GetLight(int index);
    Light& GetActiveLight() const;
	vector<shared_ptr<Light>> getAllLight();
    int GetActiveLightIndex() const;
	void SetActiveLightIndex(int index);
	
private:
	vector<shared_ptr<MeshModel>> mesh_models;
	vector<shared_ptr<Camera>> cameras;
    vector<shared_ptr<Light>> lights;

	int active_camera_index;
	int active_model_index;
    int active_light_index;
};
