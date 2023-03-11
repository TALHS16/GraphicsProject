#pragma once
#include "Scene.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include "MeshModel.h"
#include <GLFW/glfw3.h>

#define MAX_Float 3.40282e+038
#define MIN_Float 1.17549e-038
class Renderer
{
public:
	Renderer(int viewportWidth, int viewportHeight);
	virtual ~Renderer();
	void Render(const Scene& scene);
	void SwapBuffers();
	void ClearColorBuffer(const glm::vec3& color);
	int GetViewportWidth() const;
	int GetViewportHeight() const;
	void SetViewportWidth(int width) ;
	void SetViewportHeight(int height) ;
	void refreshBuff();
    void Norm(glm::vec2& v);
    void Norm3(glm::vec3& v);
    void draw_Z_buffer(GLFWwindow* window);
    void Rendrer_Z_buffer();
    bool& get_show_z_buffer();
    
private:
	void DrawCircle(const glm::vec2& center, double radius,const glm::vec3& color, int a);
	void PutPixel(const int i, const int j,const float z , const glm::vec3& color);
	void DrawLine(const glm::ivec2& p1, const glm::ivec2& p2, const glm::vec3& color,float z=MAX_Float);
	void Bresenham(const glm::ivec2& p1, const glm::ivec2& p2, const glm::vec3& color);
	void CreateBuffers(int w, int h);
	void CreateOpenglBuffer();
	void InitOpenglRendering();
    
	void DrawModel2D(MeshModel& model, int firstIndex, int secondIndex);
	void Draw_model_normals(MeshModel& model, int firstIndex, int secondIndex);
	void DrawTriangle(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& point2, const glm::vec3& color,bool flag);
	void FillTringel( glm::vec3 point0,  glm::vec3 point1,  glm::vec3 point2,  glm::vec3 color);

	float* color_buffer;
    float* zBuffer;
    float maxDepth;
    float minDepth;
	int viewport_width;
	int viewport_height;
    bool show_z_buffer;
	GLuint gl_screen_tex;
	GLuint gl_screen_vtc;
};
