#include "Scene.h"
#include "MeshModel.h"
#include <string>

Scene::Scene() :
	active_camera_index(0),
	active_model_index(0),
    active_light_index(0)
{
    
}

void Scene::RemoveModelByIndex(int index)
{
    if(index<0 || index>this->GetModelCount())
        return;
    
    mesh_models.erase(mesh_models.begin()+index);
}

void Scene::AddModel(const std::shared_ptr<MeshModel>& mesh_model)
{
    if(this->GetCameraCount()>0)
        mesh_model->activeCamera = &this->GetActiveCamera();
	mesh_models.push_back(mesh_model);
}

int Scene::GetModelCount() const
{
	return mesh_models.size();
}

MeshModel& Scene::GetModel(int index) const
{
	return *mesh_models[index];
}

MeshModel& Scene::GetActiveModel() const
{
	return *mesh_models[active_model_index];
}

void Scene::AddCamera(const std::shared_ptr<Camera>& camera)
{
	cameras.push_back(camera);
}

int Scene::GetCameraCount() const
{
	return cameras.size();
}

Camera& Scene::GetCamera(int index)
{
	return *cameras[index];
}

Camera& Scene::GetActiveCamera() const
{
	return *cameras[active_camera_index];
}

void Scene::SetActiveCameraIndex(int index)
{
	active_camera_index = index;
    for(int i=0; i<this->GetModelCount(); i++)
        mesh_models[i]->activeCamera = &GetCamera(index);
}

int Scene::GetActiveCameraIndex() const
{
	return active_camera_index;
}

void Scene::SetActiveModelIndex(int index)
{
	active_model_index = index;
}

int Scene::GetActiveModelIndex() const
{
	return active_model_index;
}

void Scene::AddLight(const shared_ptr<Light>& light)
{
    lights.push_back(light);
}
int Scene::GetLightCount() const
{
    return this->lights.size();
}
Light& Scene::GetLight(int index)
{
    return *lights[index];
}

int Scene::GetActiveLightIndex() const
{
    return active_light_index;
}

void Scene::SetActiveLightIndex(int index)
{
    this->active_light_index = index;
}

Light& Scene::GetActiveLight() const
{
    return *lights[active_light_index];
}

vector<shared_ptr<Light>> Scene::getAllLight()
{
    return lights;
}
