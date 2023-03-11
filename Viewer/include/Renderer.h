#pragma once
#include "Scene.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include "MeshModel.h"
#include <GLFW/glfw3.h>
#include "ShaderProgram.h"
#include "Texture2D.h"

#define MAX_Float 3.40282e+038
#define MIN_Float 1.17549e-038
class Renderer
{
public:
    Renderer(int viewportWidth, int viewportHeight);
	virtual ~Renderer();
	void Render(Scene& scene);
	/*void SwapBuffers();
	void ClearColorBuffer(const glm::vec3& color);
	int GetViewportWidth() const;
	int GetViewportHeight() const;
	void SetViewportWidth(int width) ;
	void SetViewportHeight(int height) ;
	void refreshBuff();
    void Norm(glm::vec2& v);*/
    
    //void draw_Z_buffer(GLFWwindow* window);
    //void Rendrer_Z_buffer();
    //bool& get_show_z_buffer();
    /*glm::vec3 calculate_Ambient(glm::vec3 model_ambient, Scene& scene);
    glm::vec3 calculate_Diffuse(glm::vec3 model_ambient,glm::vec3 point_normal,glm::vec3 point, Scene& scene);
    glm::vec3 calculate_Specular(glm::vec3 model_ambient,glm::vec3 point_normal,glm::vec3 point,float alpha, Scene& scene);
    glm::vec3 calculate_normal_in_triangle(glm::vec3& point0,glm::vec3& point1,glm::vec3& point2,
                                           glm::vec2& point, glm::vec3& normal0, glm::vec3& normal1, glm::vec3& normal2);
    glm::vec3 reflect(glm::vec3 lightDirection, glm::vec3 normal);*/

    void LoadShaders();
    void LoadTextures();

private:
	/*void DrawCircle(const glm::vec2& center, double radius,const glm::vec3& color, int a);
	void PutPixel(const int i, const int j,const float z , const glm::vec3& color);
	void DrawLine(const glm::ivec2& p1, const glm::ivec2& p2, const glm::vec3& color,float z=MAX_Float);
	void Bresenham(const glm::ivec2& p1, const glm::ivec2& p2, const glm::vec3& color);
	void CreateBuffers(int w, int h);
	void CreateOpenglBuffer();
	void InitOpenglRendering();
    float calculate_area(glm::vec3& point0,glm::vec3& point1,glm::vec2& point2);
    float calculate_z_in_triangle(glm::vec3& point0,glm::vec3& point1,glm::vec3& point2, glm::vec2 point);
	void DrawModel2D(MeshModel& model, Scene& scene);
	void Draw_model_normals(MeshModel& model, int firstIndex, int secondIndex);
	void DrawTriangle(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& point2, const glm::vec3& color,bool flag);
    void FillTringel( glm::vec3 point0,  glm::vec3 point1,  glm::vec3 point2, glm::vec3 color);
    */
    
    //void Fill_tringel_by_shading_model(glm::vec3 point0,  glm::vec3 point1,  glm::vec3 point2, glm::vec3 normal0,glm::vec3 normal1, glm::vec3 normal2,Material model_materials,bool gouraud_shading, Scene& scene);

	/*float* color_buffer;
    float* zBuffer;
    float maxDepth;
    float minDepth;
	int viewport_width;
	int viewport_height;
    bool show_z_buffer;*/
    int viewport_width;
    int viewport_height;

	GLuint gl_screen_tex;
	GLuint gl_screen_vtc;

    ShaderProgram lightShader;
    ShaderProgram colorShader;
    ShaderProgram cubeShader;
    TextureCube myCube;
    int cubemapTexture;
    GLuint cubeVAO;
    Texture2D texture1;
    Texture2D texture2;
    Texture2D texture3;
    Texture2D enviroText;
};

static void Norm3(glm::vec3& v);

static void Norm3(glm::vec3& v)
{
    v[0] = (v[0] + 1) * (Camera::screen_width / 2);
    v[1] = (v[1] + 1) * (Camera::screen_hight / 2);
}
//
static glm::vec3 Norm3Copy(glm::vec3 v)
{
    v[0] = (v[0] + 1) * (Camera::screen_width / 2);
    v[1] = (v[1] + 1) * (Camera::screen_hight / 2);
    return v;
}
//
static glm::vec3 Rev_Norm3Copy(glm::vec3 v)
{
    v[0] = (v[0] / (Camera::screen_width / 2)) - 1 ;
    v[1] = (v[1] / (Camera::screen_hight / 2)) - 1;
    return v;
}
