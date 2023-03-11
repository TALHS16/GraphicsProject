#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>

#include "Renderer.h"
#include "InitShader.h"

#define INDEX(width,x,y,c) ((x)+(y)*(width))*3+(c)
#define Z_INDEX(width,x,y) ((x)+(y)*(width))

float minINF = MIN_Float;
float maxINF = MAX_Float;

Renderer::Renderer(int viewport_width, int viewport_height) :
	viewport_width(viewport_width),
	viewport_height(viewport_height)
{
	/*InitOpenglRendering();
	CreateBuffers(viewport_width, viewport_height);*/
	LoadShaders();
	LoadTextures();
	cubemapTexture = myCube.loadCubemap();

    float cubeVertices[] = {
        // positions          
        -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        -1.0f,  1.0f, -1.0f,
         1.0f,  1.0f, -1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f,  1.0f
    };

    // skybox VAO
    GLuint cubeVBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glBindVertexArray(cubeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), &cubeVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
}

Renderer::~Renderer()
{
	/*delete[] color_buffer;
    delete[] zBuffer;*/
}

//void Renderer::PutPixel(const int i, const int j,const float z , const glm::vec3& color)
//{
//	if (i < 0) return; if (i >= viewport_width) return;
//	if (j < 0) return; if (j >= viewport_height) return;
//    if (zBuffer[Z_INDEX(viewport_width, i, j)] > z)
//        return;
//    else{
//        zBuffer[Z_INDEX(viewport_width, i, j)] = z;
//        this->maxDepth = glm::max(z, maxDepth);
//        this->minDepth = glm::min(z, minDepth);
//    }
//        
//	color_buffer[INDEX(viewport_width, i, j, 0)] = color.x;
//	color_buffer[INDEX(viewport_width, i, j, 1)] = color.y;
//	color_buffer[INDEX(viewport_width, i, j, 2)] = color.z;
//}

//void Renderer::Rendrer_Z_buffer()
//{
//    for(int i=0;i<viewport_width;i++){
//        for(int j=0;j<viewport_height;j++){
//            float z = zBuffer[Z_INDEX(viewport_width, i, j)];
//            float color;
//            color = float((z-this->minDepth)/(this->maxDepth - this->minDepth));
//            color_buffer[INDEX(viewport_width, i, j, 0)] = color;
//            color_buffer[INDEX(viewport_width, i, j, 1)] = color;
//            color_buffer[INDEX(viewport_width, i, j, 2)] = color;
//        }
//    }
//}



//void Renderer::DrawLine(const glm::ivec2& p1, const glm::ivec2& p2, const glm::vec3& color, float z)
//{
//	double deltaP = p2[0] - p1[0];
//	double deltaQ = p2[1] - p1[1];
//	double x = p1[0];
//	double y = p1[1];
//	double x2 = p2[0];
//	double y2 = p2[1];
//	double m;
//    
//	if (deltaP != 0)
//		m = deltaQ / deltaP;
//	else
//		m = deltaQ > 0 ? 2 : -2;
//	double e;
//	int c = 1;
//	if (m < 0)
//	{
//		c = -1;
//	}
//	if (m <= 1 && m >= -1)
//	{
//
//		if (p1[0] > p2[0])
//		{
//			x = p2[0];
//			y = p2[1];
//			x2 = p1[0];
//			deltaP *= -1;
//			deltaQ *= -1;
//		}
//		e = -deltaP;
//		while (x <= x2)
//		{
//			if (e > 0)
//			{
//				y += c;
//				e -= 2 * deltaP;
//			}
//			PutPixel((int)x, (int)y,z, color);
//			x++;
//			e += c * 2 * deltaQ;
//		}
//	}
//	else
//	{
//		if (p1[1] > p2[1])
//		{
//			x = p2[0];
//			y = p2[1];
//			y2 = p1[1];
//			deltaP *= -1;
//			deltaQ *= -1;
//		}
//		e = -deltaQ;
//		while (y <= y2)
//		{
//			if (e > 0)
//			{
//				x += c;
//				e -= 2 * deltaQ;
//			}
//			PutPixel((int)x, (int)y,z, color);
//			y++;
//			e += c * 2 * deltaP;
//		}
//	}
//
//
//}
//
//bool comparePoints(glm::vec2 point0, glm::vec2 point1) {
//	return point0[1] < point1[1];
//}



//glm::vec3 Renderer::calculate_Ambient(glm::vec3 model_ambient, Scene& scene)
//{
//    for(int i=0; i< scene.GetLightCount();i++){
//        model_ambient += model_ambient * scene.GetLight(i).materials.get_ambient();
//    }
//    return model_ambient;
//}
//
//glm::vec3 Renderer::calculate_Diffuse(glm::vec3 model_ambient,glm::vec3 point_normal,glm::vec3 point, Scene& scene)
//{
//    glm::vec3 diffuse = {0,0,0};
//    for(int i=0; i< scene.GetLightCount();i++){
//        glm::vec3 l = glm::normalize(scene.GetLight(i).get_position() - point);
//        diffuse += model_ambient * scene.GetLight(i).materials.get_diffuse() * max(glm::dot(l,point_normal),0.0f);
//    }
//    return diffuse;
//    
//}
//
//glm::vec3 Renderer::calculate_Specular(glm::vec3 model_ambient,glm::vec3 point_normal,glm::vec3 point,float alpha, Scene& scene)
//{
//    glm::vec3 specular = {0,0,0};
//    glm::vec3 v = glm::normalize(Norm3Copy(scene.GetActiveCamera().get_eye()) - point);
//    for(int i=0; i< scene.GetLightCount();i++){
//        glm::vec3 l = glm::normalize(scene.GetLight(i).get_position() - point);
//        glm::vec3 r = reflect(l, point_normal);
//        specular += model_ambient * scene.GetLight(i).materials.get_specular() * pow(max(glm::dot(v,r), 0.0f), alpha);
//    }
//    return specular;
//}
//
//
//
//glm::vec3 Renderer::calculate_normal_in_triangle(glm::vec3& point0,glm::vec3& point1,glm::vec3& point2,
//                                                 glm::vec2& point, glm::vec3& normal0, glm::vec3& normal1, glm::vec3& normal2)
//{
//    float area12 = calculate_area(point1, point2, point);
//    float area02 = calculate_area(point0, point2, point);
//    float area01 = calculate_area(point0, point1, point);
//    
//    float sum = area01 + area02 + area12;
//    return normal0 * (area12/sum) + normal1*(area02/sum) + normal2*(area01/sum);
//}
//
//glm::vec3 Renderer::reflect(glm::vec3 lightDirection, glm::vec3 normal)
//{
//    return glm::normalize(2.0f * normal * glm::dot(lightDirection, normal) - lightDirection);
//}

//void Renderer::Fill_tringel_by_shading_model(glm::vec3 point0,  glm::vec3 point1,  glm::vec3 point2, glm::vec3 normal0,
//                              glm::vec3 normal1, glm::vec3 normal2,Material model_materials,bool gouraud_shading, Scene& scene)
//{
//    vector<glm::vec2> vec = {point0,point1,point2 };
//
//    // Sort points by y-coordinate
//    sort(vec.begin(), vec.end(), comparePoints);
//    
//    float z;
//    
//    float y0 = vec[0][1];
//    float y1 = vec[1][1];
//    float y2 = vec[2][1];
//    
//    float a0 = (vec[0][1] - vec[1][1]) / (vec[0][0] - vec[1][0]);
//    float a1 = (vec[1][1] - vec[2][1]) / (vec[1][0] - vec[2][0]);
//    float a2 = (vec[2][1] - vec[0][1]) / (vec[2][0] - vec[0][0]);
//    float b0 = vec[0][1] - (a0 * vec[0][0]);
//    float b1 = vec[1][1] - (a1 * vec[1][0]);
//    float b2 = vec[2][1] - (a2 * vec[2][0]);
//    int xl;
//    int xr;
//    glm::vec3 color;
//    glm::vec3 color0;
//    glm::vec3 color1;
//    glm::vec3 color2;
//    if(gouraud_shading == true){
//        glm::vec3 ambient = this->calculate_Ambient(model_materials.get_ambient(), scene);
//        glm::vec3 diffuse = this->calculate_Diffuse(model_materials.get_diffuse(), normal0,point0 ,scene);
//        glm::vec3 specular = this->calculate_Specular(model_materials.get_specular(), normal0,point0,*model_materials.get_shine_factor(), scene);
//        
//        color0 = ambient + diffuse + specular;
//        diffuse = this->calculate_Diffuse(model_materials.get_diffuse(), normal1, point1,scene);
//        specular = this->calculate_Specular(model_materials.get_specular(), normal1,point1,*model_materials.get_shine_factor(), scene);
//        
//        color1 = ambient + diffuse + specular;
//        
//        diffuse = this->calculate_Diffuse(model_materials.get_diffuse(), normal2,point2, scene);
//        specular = this->calculate_Specular(model_materials.get_specular(), normal2,point2,*model_materials.get_shine_factor(), scene);
//        
//        color2 = ambient + diffuse + specular;
//    }
//        // Rasterize bottom part of triangle
//    for (float y = y0; y < y1; y++) {
//        if ((vec[0][0] - vec[1][0]) == 0)
//            xl = vec[1][0];
//        else
//            xl = (y - b0) / a0;
//        if ((vec[2][0] - vec[0][0]) == 0)
//            xr = vec[0][0];
//        else
//            xr = (y - b2) / a2;
//        const glm::vec2& temp0 = { xl,(int)y };
//        const glm::vec2& temp1 = { xr,(int)y };
//        float xMax = glm::max(xl, xr);
//        float xMin = glm::min(xl, xr);
//        for(;xMin<xMax;xMin++){
//            glm::vec2 newPoint = glm::vec2(xMin,y);
//            z = calculate_z_in_triangle(point0, point1, point2, {xMin, y});
//            if(gouraud_shading == false){ // phong shading
//                glm::vec3 point_normal = calculate_normal_in_triangle(point0, point1, point2, newPoint, normal0, normal1, normal2);
//                glm::vec3 ambient = this->calculate_Ambient(model_materials.get_ambient(), scene);
//                glm::vec3 diffuse = this->calculate_Diffuse(model_materials.get_diffuse(), point_normal,Rev_Norm3Copy({xMin,y,z}), scene);
//                glm::vec3 specular = this->calculate_Specular(model_materials.get_specular(), point_normal,Rev_Norm3Copy({xMin,y,z}),*model_materials.get_shine_factor(), scene);
//                color = ambient + diffuse + specular;
//            }
//            else{ // gouraud shading
//                color = calculate_normal_in_triangle(point0, point1, point2, newPoint, color0, color1, color2);
//            }
//            
//            PutPixel(xMin, y, z, color);
//        }
//    }
//    // Rasterize top part of triangle
//    for (float y = y1; y < y2; y++) {
//        if ((vec[1][0] - vec[2][0]) == 0)
//            xl = vec[1][0];
//        else
//            xl = (y - b1) / a1;
//        if ((vec[2][0] - vec[0][0]) == 0)
//            xr = vec[0][0];
//        else
//            xr = (y - b2) / a2;
//        const glm::vec2& temp0 = { xl,(int)y };
//        const glm::vec2& temp1 = { xr,(int)y };
//        float xMax = glm::max(xl, xr);
//        float xMin = glm::min(xl, xr);
//        for(;xMin<xMax;xMin++){
//            glm::vec2 newPoint = glm::vec2(xMin,y);
//            z = calculate_z_in_triangle(point0, point1, point2, {xMin, y});
//            if(gouraud_shading == false){ // phong shading
//                glm::vec3 point_normal = calculate_normal_in_triangle(point0, point1, point2, newPoint, normal0, normal1, normal2);
//                glm::vec3 ambient = this->calculate_Ambient(model_materials.get_ambient(), scene);
//                glm::vec3 diffuse = this->calculate_Diffuse(model_materials.get_diffuse(), point_normal, Rev_Norm3Copy({xMin,y,z}), scene);
//                glm::vec3 specular = this->calculate_Specular(model_materials.get_specular(), point_normal,Rev_Norm3Copy({xMin,y,z}), *model_materials.get_shine_factor(), scene);
//                color = ambient + diffuse + specular;
//            }
//            else{ // gouraud shading
//                color = calculate_normal_in_triangle(point0, point1, point2, newPoint, color0, color1, color2);
//            }
//            PutPixel(xMin, y, z, color);
//        }
//    }
//}

//void Renderer::FillTringel(glm::vec3 point0,  glm::vec3 point1,  glm::vec3 point2, glm::vec3 color)
//{
//	vector<glm::vec2> vec = {point0,point1,point2 };
//
//	// Sort points by y-coordinate
//	sort(vec.begin(), vec.end(), comparePoints);
//    
//    float z;
//    
//	float y0 = vec[0][1];
//	float y1 = vec[1][1];
//	float y2 = vec[2][1];
//	
//	float a0 = (vec[0][1] - vec[1][1]) / (vec[0][0] - vec[1][0]);
//	float a1 = (vec[1][1] - vec[2][1]) / (vec[1][0] - vec[2][0]);
//	float a2 = (vec[2][1] - vec[0][1]) / (vec[2][0] - vec[0][0]);
//	float b0 = vec[0][1] - (a0 * vec[0][0]);
//	float b1 = vec[1][1] - (a1 * vec[1][0]);
//	float b2 = vec[2][1] - (a2 * vec[2][0]);
//	int xl;
//	int xr;
//		// Rasterize bottom part of triangle
//		for (float y = y0; y < y1; y++) {
//			if ((vec[0][0] - vec[1][0]) == 0)
//				xl = vec[1][0];
//			else
//				xl = (y - b0) / a0;
//			if ((vec[2][0] - vec[0][0]) == 0)
//				xr = vec[0][0];
//			else
//				xr = (y - b2) / a2;
//			const glm::vec2& temp0 = { xl,(int)y };
//			const glm::vec2& temp1 = { xr,(int)y };
//            float xMax = glm::max(xl, xr);
//            float xMin = glm::min(xl, xr);
//            for(;xMin<xMax;xMin++){
//                glm::vec2 newPoint = glm::vec2(xMin,y);
//                z = calculate_z_in_triangle(point0, point1, point2, newPoint);
//                PutPixel(xMin, y, z, color);
//            }
//		}
//	// Rasterize top part of triangle
//	for (float y = y1; y < y2; y++) {
//		if ((vec[1][0] - vec[2][0]) == 0)
//			xl = vec[1][0];
//		else
//			xl = (y - b1) / a1;
//		if ((vec[2][0] - vec[0][0]) == 0)
//			xr = vec[0][0];
//		else
//			xr = (y - b2) / a2;
//		const glm::vec2& temp0 = { xl,(int)y };
//		const glm::vec2& temp1 = { xr,(int)y };
//        float xMax = glm::max(xl, xr);
//        float xMin = glm::min(xl, xr);
//        for(;xMin<xMax;xMin++){
//            glm::vec2 newPoint = glm::vec2(xMin,y);
//            z = calculate_z_in_triangle(point0, point1, point2,newPoint);
//            PutPixel(xMin, y, z, color);
//        }
//	}
//}
//float Renderer::calculate_area(glm::vec3& point0,glm::vec3& point1,glm::vec2& point2)
//{
//    float a = point2.x*(point0.y-point1.y) + point0.x*(point1.y - point2.y) + point1.x*(point2.y - point0.y);
//    if(a>0)
//        return 0.5*a;
//    else
//        return -0.5*a;
//}
//
//float Renderer::calculate_z_in_triangle(glm::vec3& point0,glm::vec3& point1,glm::vec3& point2, glm::vec2 point)
//{
//    float area12 = calculate_area(point1, point2, point);
//    float area02 = calculate_area(point0, point2, point);
//    float area01 = calculate_area(point0, point1, point);
//    
//    float sum = area01 + area02 + area12;
//    return point0.z * (area12/sum) + point1.z*(area02/sum) + point2.z*(area01/sum);
//}


//void Renderer::DrawCircle(const glm::vec2& center, double radius,const glm::vec3& color, int a=10000)
//{
//	for(double i=0; i < a; i+=1){
//		glm::ivec2 p2(center[0] + radius*sin((2*M_PI*i)/a), center[1] + radius * cos((2 * M_PI * i) / a));
//		DrawLine(center, p2, color);
//	}
//}

//void Renderer::CreateBuffers(int w, int h)
//{
//	CreateOpenglBuffer(); //Do not remove this line.
//	color_buffer = new float[3 * w * h];
//    zBuffer = new float[w * h];
//	ClearColorBuffer(glm::vec3(0.0f, 0.0f, 0.0f));
//    
//    
//    
//}

//##############################
//##OpenGL stuff. Don't touch.##
//##############################

// Basic tutorial on how opengl works:
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-2-the-first-triangle/
// don't linger here for now, we will have a few tutorials about opengl later.
//void Renderer::InitOpenglRendering()
//{
//	// Creates a unique identifier for an opengl texture.
//	glGenTextures(1, &gl_screen_tex);
//
//	// Same for vertex array object (VAO). VAO is a set of buffers that describe a renderable object.
//	glGenVertexArrays(1, &gl_screen_vtc);
//
//	GLuint buffer;
//
//	// Makes this VAO the current one.
//	glBindVertexArray(gl_screen_vtc);
//
//	// Creates a unique identifier for a buffer.
//	glGenBuffers(1, &buffer);
//
//	// (-1, 1)____(1, 1)
//	//	     |\  |
//	//	     | \ | <--- The exture is drawn over two triangles that stretch over the screen.
//	//	     |__\|
//	// (-1,-1)    (1,-1)
//	const GLfloat vtc[]={
//		-1, -1,
//		 1, -1,
//		-1,  1,
//		-1,  1,
//		 1, -1,
//		 1,  1
//	};
//
//	const GLfloat tex[]={
//		0,0,
//		1,0,
//		0,1,
//		0,1,
//		1,0,
//		1,1};
//
//	// Makes this buffer the current one.
//	glBindBuffer(GL_ARRAY_BUFFER, buffer);
//
//	// This is the opengl way for doing malloc on the gpu. 
//	glBufferData(GL_ARRAY_BUFFER, sizeof(vtc)+sizeof(tex), NULL, GL_STATIC_DRAW);
//
//	// memcopy vtc to buffer[0,sizeof(vtc)-1]
//	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vtc), vtc);
//
//	// memcopy tex to buffer[sizeof(vtc),sizeof(vtc)+sizeof(tex)]
//	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vtc), sizeof(tex), tex);
//
//	// Loads and compiles a sheder.
//	GLuint program = InitShader( "vshader.glsl", "fshader.glsl" );
//
//	// Make this program the current one.
//	glUseProgram(program);
//
//	// Tells the shader where to look for the vertex position data, and the data dimensions.
//	GLint  vPosition = glGetAttribLocation( program, "vPosition" );
//	glEnableVertexAttribArray( vPosition );
//	glVertexAttribPointer( vPosition,2,GL_FLOAT,GL_FALSE,0,0 );
//
//	// Same for texture coordinates data.
//	GLint  vTexCoord = glGetAttribLocation( program, "vTexCoord" );
//	glEnableVertexAttribArray( vTexCoord );
//	glVertexAttribPointer( vTexCoord,2,GL_FLOAT,GL_FALSE,0,(GLvoid *)sizeof(vtc) );
//
//	//glProgramUniform1i( program, glGetUniformLocation(program, "texture"), 0 );
//
//	// Tells the shader to use GL_TEXTURE0 as the texture id.
//	glUniform1i(glGetUniformLocation(program, "texture"),0);
//}
//
//void Renderer::CreateOpenglBuffer()
//{
//	// Makes GL_TEXTURE0 the current active texture unit
//	glActiveTexture(GL_TEXTURE0);
//
//	// Makes glScreenTex (which was allocated earlier) the current texture.
//	glBindTexture(GL_TEXTURE_2D, gl_screen_tex);
//
//	// malloc for a texture on the gpu.
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, viewport_width, viewport_height, 0, GL_RGB, GL_FLOAT, NULL);
//	glViewport(0, 0, viewport_width, viewport_height);
//}
//
//void Renderer::SwapBuffers()
//{
//	// Makes GL_TEXTURE0 the current active texture unit
//	glActiveTexture(GL_TEXTURE0);
//
//	// Makes glScreenTex (which was allocated earlier) the current texture.
//	glBindTexture(GL_TEXTURE_2D, gl_screen_tex);
//
//	// memcopy's colorBuffer into the gpu.
//	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, viewport_width, viewport_height, GL_RGB, GL_FLOAT, color_buffer);
//
//	// Tells opengl to use mipmapping
//	glGenerateMipmap(GL_TEXTURE_2D);
//
//	// Make glScreenVtc current VAO
//	glBindVertexArray(gl_screen_vtc);
//
//	// Finally renders the data.
//	glDrawArrays(GL_TRIANGLES, 0, 6);
//}
//
//void Renderer::ClearColorBuffer(const glm::vec3& color)
//{
//	for (int i = 0; i < viewport_width; i++)
//	{
//		for (int j = 0; j < viewport_height; j++)
//		{
//            color_buffer[INDEX(viewport_width, i, j, 0)] = color.x;
//            color_buffer[INDEX(viewport_width, i, j, 1)] = color.y;
//            color_buffer[INDEX(viewport_width, i, j, 2)] = color.z;
//            zBuffer[Z_INDEX(viewport_width, i, j)] = minINF;
//            this->maxDepth = minINF;
//            this->minDepth = maxINF;
//		}
//	}
//}

//void Renderer::DrawTriangle(const glm::vec3& point0, const glm::vec3& point1, const glm::vec3& point2, const glm::vec3& color,bool flag)
//{
//    if(flag)
//        FillTringel(point0, point1, point2, color);
//    else{
//        DrawLine(point0, point1, color);
//        DrawLine(point1, point2, color);
//        DrawLine(point0, point2, color);
//    }
//
//}
//void Renderer::DrawModel2D(MeshModel& model, Scene& scene)
//{
//	// TODO: over all the vertices and draw the lines
//    int firstIndex = 0;
//    int secondIndex = 1;
//    
//	for( int i=0;i<model.GetFacesCount();i++){
//		int point0_index = model.GetFace(i).GetVertexIndex(0);
//		int point1_index = model.GetFace(i).GetVertexIndex(1);
//		int point2_index = model.GetFace(i).GetVertexIndex(2);
//
//		glm::vec3 point0_3 = model.get_vertex3_after_transform(point0_index-1);
//		glm::vec3 point1_3 = model.get_vertex3_after_transform(point1_index-1);
//		glm::vec3 point2_3 = model.get_vertex3_after_transform(point2_index-1);
//        
//        glm::vec2 point0 = glm::vec2(point0_3);
//        glm::vec2 point1 = glm::vec2(point1_3);
//        glm::vec2 point2 = glm::vec2(point2_3);
//        
//		
//        
//        if(model.get_show_object()){
//            const glm::vec3 model_color =  std::as_const(model.get_model_color());
//            Norm(point0);
//            Norm(point1);
//            Norm(point2);
//            Norm3(point0_3);
//            Norm3(point1_3);
//            Norm3(point2_3);
//            if(model.get_show_bound_triangle()){
//                float avgZ = (point0_3[2] + point1_3[2] + point2_3[2])/3.0;
//                float deph = (avgZ- model.zMin) / (model.zMax - model.zMin);
//                float left = std::min(point0.x,std::min(point1.x, point2.x));
//                float right = std::max(point0.x,std::max(point1.x, point2.x));
//                float down = std::min(point0.y,std::min(point1.y, point2.y));
//                float up = std::max(point0.y,std::max(point1.y, point2.y));
//                glm::vec3 color_deph = glm::vec3(deph,deph,deph);
//                DrawLine(glm::vec2(left, down), glm::vec2(left, up),color_deph);
//                DrawLine(glm::vec2(left, up), glm::vec2(right, up), color_deph);
//                DrawLine(glm::vec2(right, up), glm::vec2(right, down),color_deph);
//                DrawLine(glm::vec2(right, down), glm::vec2(left, down),color_deph);
//            }
//            if(model.get_fill_triangle()==false && (model.get_show_phong_shading() || model.get_show_show_gouraud_shading())){
//                int normal0_index = model.GetFace(i).GetNormalIndex(0);
//                int normal1_index = model.GetFace(i).GetNormalIndex(1);
//                int normal2_index = model.GetFace(i).GetNormalIndex(2);
//                glm::vec3 normal0 = model.get_normal(normal0_index-1);
//                glm::vec3 normal1 = model.get_normal(normal1_index-1);
//                glm::vec3 normal2 = model.get_normal(normal2_index-1);
//                Fill_tringel_by_shading_model(point0_3, point1_3, point2_3, normal0,
//                                            normal1, normal2,model.materials,model.get_show_show_gouraud_shading(), scene);
//            }
//
//            if(model.get_fill_triangle()==true && !(model.get_show_phong_shading() || model.get_show_show_gouraud_shading()))
//				DrawTriangle(point0_3, point1_3, point2_3, model.get_model_fill_color()[i], model.get_fill_triangle());
//			if(model.get_fill_triangle()==false && !(model.get_show_phong_shading() || model.get_show_show_gouraud_shading()))
//				DrawTriangle(point0_3, point1_3, point2_3, model.get_model_color(),model.get_fill_triangle());
//        }
//		
//        if(model.get_show_faces_normals()){
//			glm::vec2 normal_tranform = model.get_normal_of_face(i);
//            DrawLine((point0+point1+point2)/float(3.0), (point0+point1+point2)/float(3.0) + normal_tranform*float(model.get_scale_factor()*Camera::screen_hight/2), glm::vec3(0,255,0));
//        }
//        
//        if(model.get_show_refelections()){
//            glm::vec3 normal_tranform = model.get_normal_of_face(i);
//            glm::vec3 center = (point0_3+point1_3+point2_3)/float(3.0);
//            glm::vec3 l = glm::normalize(scene.GetActiveLight().get_position() - center);
//            glm::vec3 r = reflect(l, normal_tranform);
//            DrawLine((point0+point1+point2)/float(3.0), (point0+point1+point2)/float(3.0) + glm::vec2(r)*float(model.get_scale_factor()), glm::vec3(0,255,0));
//        }
//        
//        if(model.get_show_light_direction()){
//            glm::vec3 normal_tranform = model.get_normal_of_face(i);
//            glm::vec3 center = (point0_3+point1_3+point2_3)/float(3.0);
//            glm::vec3 l = glm::normalize(scene.GetActiveLight().get_position() - center);
//            DrawLine((point0+point1+point2)/float(3.0), (point0+point1+point2)/float(3.0) + glm::vec2(l)*float(model.get_scale_factor()), glm::vec3(0,255,0));
//        }
//        
//        if(model.get_show_verteces_normals()){
//            int normal0_index = model.GetFace(i).GetNormalIndex(0);
//            int normal1_index = model.GetFace(i).GetNormalIndex(1);
//            int normal2_index = model.GetFace(i).GetNormalIndex(2);
//            glm::vec2 normal0 = model.get_normal(normal0_index-1);
//            glm::vec2 normal1 = model.get_normal(normal1_index-1);
//            glm::vec2 normal2 = model.get_normal(normal2_index-1);
//            
//            DrawLine(point0,point0 + normal0*float(model.get_scale_factor()*Camera::screen_hight/4), glm::vec3(0,255,0));
//            DrawLine(point1,point1 + normal1*float(model.get_scale_factor()*Camera::screen_hight/4), glm::vec3(0,255,0));
//            DrawLine(point2,point2 + normal2*float(model.get_scale_factor()*Camera::screen_hight/4), glm::vec3(0,255,0));
//        }
//	}
//    
//    
//    if(model.get_show_axes()){
//        glm::vec2 local0 = model.world_transform_point(glm::vec3(0,0,0))*Camera::screen_width;
//        glm::vec2 localX = model.world_transform_point(glm::vec3(1,0,0))*Camera::screen_width;
//        glm::vec2 localY = model.world_transform_point(glm::vec3(0,1,0))*Camera::screen_width;
//        glm::vec2 localZ = model.world_transform_point(glm::vec3(0,0,1))*Camera::screen_width;
//        Norm(local0);
//        Norm(localX);
//        Norm(localY);
//        Norm(localZ);
//        DrawLine(local0, localX, glm::vec3(255,0,0));
//        DrawLine(local0, localY, glm::vec3(255,0,0));
//        DrawLine(local0, localZ, glm::vec3(255,0,0));
//        
//        glm::vec2 world0 = glm::vec2(10,10);
//        glm::vec2 worldX = glm::vec2(500,10);
//        glm::vec2 worldY = glm::vec2(10,500);
//
//        DrawLine(world0, worldX, glm::vec3(0,255,0));
//        DrawLine(world0, worldY, glm::vec3(0,255,0));
//        
//        //Draw local box
//        
//        glm::mat4 mul_mat = model.activeCamera->final_matrix() * model.get_world_transform_matrix() * glm::inverse(model.activeCamera->final_matrix());
//        
//        glm::vec4 XYZ4 =  mul_mat * glm::vec4(model.right, model.up, model.zMax, 1);
//        glm::vec4 XYz4 =  mul_mat * glm::vec4(model.right, model.up, model.zMin, 1);
//        glm::vec4 XyZ4 =  mul_mat * glm::vec4(model.right, model.down, model.zMax, 1);
//        glm::vec4 Xyz4 =  mul_mat * glm::vec4(model.right, model.down, model.zMin, 1);
//        glm::vec4 xyz4 =  mul_mat * glm::vec4(model.left, model.down, model.zMin, 1);
//        glm::vec4 xyZ4 =  mul_mat * glm::vec4(model.left, model.down, model.zMax, 1);
//        glm::vec4 xYz4 =  mul_mat * glm::vec4(model.left, model.up, model.zMin, 1);
//        glm::vec4 xYZ4 =  mul_mat * glm::vec4(model.left, model.up, model.zMax, 1);
//        
//        glm::vec2 XYZ = glm::vec2(XYZ4)/XYZ4[3];
//        glm::vec2 XYz = glm::vec2(XYz4)/XYz4[3];
//        glm::vec2 XyZ = glm::vec2(XyZ4)/XyZ4[3];
//        glm::vec2 Xyz = glm::vec2(Xyz4)/Xyz4[3];
//        glm::vec2 xYZ = glm::vec2(xYZ4)/xYZ4[3];
//        glm::vec2 xYz = glm::vec2(xYz4)/xYz4[3];
//        glm::vec2 xyZ = glm::vec2(xyZ4)/xyZ4[3];
//        glm::vec2 xyz = glm::vec2(xyz4)/xyz4[3];
//        
//        Norm(XYZ);
//        Norm(XYz);
//        Norm(XyZ);
//        Norm(Xyz);
//        Norm(xYZ);
//        Norm(xYz);
//        Norm(xyZ);
//        Norm(xyz);
//
//        
//        DrawLine(XYZ, XyZ, glm::vec3(255,0,0));
//        DrawLine(XyZ, xyZ, glm::vec3(255,0,0));
//        DrawLine(xyZ, xYZ, glm::vec3(255,0,0));
//        DrawLine(xYZ, XYZ, glm::vec3(255,0,0));
//        
//        DrawLine(XYz, Xyz, glm::vec3(255,0,0));
//        DrawLine(Xyz, xyz, glm::vec3(255,0,0));
//        DrawLine(xyz, xYz, glm::vec3(255,0,0));
//        DrawLine(xYz, XYz, glm::vec3(255,0,0));
//        
//        DrawLine(xYZ, xYz, glm::vec3(255,0,0));
//        DrawLine(xyZ, xyz, glm::vec3(255,0,0));
//        
//        DrawLine(XYZ, XYz, glm::vec3(255,0,0));
//        DrawLine(XyZ, Xyz, glm::vec3(255,0,0));
//    }
//    
//    
//
//}
//
//void Renderer::Draw_model_normals(MeshModel& model, int firstIndex, int secondIndex)
//{
//    // TODO: over all the vertices and draw the lines
//
//    for( int i=0;i<model.GetNormalsCount();i++){
//        int point0_index = model.GetFace(i).GetVertexIndex(0);
//        int point1_index = model.GetFace(i).GetVertexIndex(1);
//        int point2_index = model.GetFace(i).GetVertexIndex(2);
//
//        glm::vec2 point0 = model.get_vertex2(point0_index-1, firstIndex, secondIndex);
//        glm::vec2 point1 = model.get_vertex2(point1_index-1, firstIndex, secondIndex);
//        glm::vec2 point2 = model.get_vertex2(point2_index-1, firstIndex, secondIndex);
//
//        const glm::vec3 model_color =  std::as_const(model.get_model_color());
//
//        Norm(point0);
//        Norm(point1);
//        Norm(point2);
//
//        if ((point0[0] < GetViewportWidth() && GetViewportHeight() > point0[1]) || (point1[0] < GetViewportWidth() && GetViewportHeight() > point1[1]))
//            DrawLine(point0, point1, model_color);
//        if ((point1[0] < GetViewportWidth() && GetViewportHeight() > point1[1]) || (point2[0] < GetViewportWidth() && GetViewportHeight() > point2[1]))
//            DrawLine(point1, point2, model_color);
//        if ((point2[0] < GetViewportWidth() && GetViewportHeight() > point2[1]) || (point0[0] < GetViewportWidth() && GetViewportHeight() > point0[1]))
//            DrawLine(point0, point2, model_color);
//
//        //DrawLine(point0, point1, model_color);
//        //DrawLine(point1, point2, model_color);
//        //DrawLine(point0, point2, model_color);
//    }
//
//}
//
//bool& Renderer::get_show_z_buffer(){
//    return this->show_z_buffer;
//}

void Renderer::Render(Scene& scene)
{
	scene.GetActiveCamera().refresh_transform_matrix();
	glDepthMask(GL_FALSE);
	cubeShader.use();
	cubeShader.setUniform("skybox", 0);
	cubeShader.setUniform("view", glm::mat4(glm::mat3(scene.GetActiveCamera().GetViewTransformation())));
	cubeShader.setUniform("projection", scene.GetActiveCamera().GetProjectionTransformation());
	glBindVertexArray(cubeVAO);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
	glDrawArrays(GL_TRIANGLES, 0, 36);
	glBindVertexArray(0);
	glDepthFunc(GL_LESS);
	
	//draw all the models in the scene
	for (int i = 0; i < scene.GetModelCount(); i++) {
		MeshModel model = scene.GetModel(i);
		model.refresh_transform_matrix();
		//model.set_final_transform_matrix(scene.GetActiveCamera().final_matrix() * (model.get_world_transform_matrix() * model.get_local_transform_matrix()));
		//DrawModel2D(model, scene);
		colorShader.use();

		colorShader.setUniform("model", model.get_world_transform_matrix() * model.get_local_transform_matrix());
		colorShader.setUniform("view", scene.GetActiveCamera().GetViewTransformation());
		colorShader.setUniform("projection", scene.GetActiveCamera().GetProjectionTransformation());
		colorShader.setUniform("material.textureMap", 0);
		//colorShader.setUniform("enviroMap", 1);
		//colorShader.setUniform("diffuseMap", 1);
		//colorShader.setUniform("normalMap", 2);
		colorShader.setUniformInt("diffuseMap", 0);
		colorShader.setUniformInt("normalMap", 1);
		colorShader.setUniform("texture_type", model.get_curr_text_type());
		colorShader.setUniform("mapping_type", model.get_curr_mapping());

		colorShader.setUniform("camera_pos", glm::inverse(scene.GetActiveCamera().get_world_transform_matrix() * scene.GetActiveCamera().get_local_transform_matrix()));
		colorShader.setUniform("material.ambient", model.materials.get_ambient());
		colorShader.setUniform("material.diffuse", model.materials.get_diffuse());
		colorShader.setUniform("material.specular", model.materials.get_specular());
		colorShader.setUniform("material.shine_factor", *model.materials.get_shine_factor());

		int light_count = scene.getAllLight().size();
		colorShader.setUniform("lights_size", light_count);


		for (int i = 0; i < scene.getAllLight().size(); i++)
		{
			string name = "lights[" + std::to_string(i) + "]";

			colorShader.setUniform((name + ".position").c_str(), scene.GetLight(i).get_position());
			colorShader.setUniform((name + ".ambient").c_str(), scene.GetLight(i).materials.get_ambient());
			colorShader.setUniform((name + ".diffuse").c_str(), scene.GetLight(i).materials.get_diffuse());
			colorShader.setUniform((name + ".specular").c_str(), scene.GetLight(i).materials.get_specular());
		}

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
		glBindVertexArray(1);

		//// Create a buffer object to hold the data for the lights.
		//GLuint lightsBuffer;
		//glGenBuffers(1, &lightsBuffer);

		//// Get the size of the lights vector in bytes.
		//size_t lightsSize = scene.getAllLight().size() * sizeof(Light);

		//glBindBuffer(GL_UNIFORM_BUFFER, lightsBuffer);
		//glBufferData(GL_UNIFORM_BUFFER, lightsSize, scene.getAllLight().data(), GL_DYNAMIC_DRAW);

		// Allocate memory for the buffer object.
		//glBufferData(lightsBuffer, lightsSize, nullptr, GL_DYNAMIC_DRAW);

		// Copy the data from the lights vector to the buffer object.
		//glBufferSubData(lightsBuffer, 0, lightsSize, scene.getAllLight().data());

		// Bind the buffer object to a uniform block in the shader.
		/*GLuint lightsBindingPoint = 0;
		GLintptr lightsOffset = 0;*/

		//glBindBufferRange(GL_UNIFORM_BUFFER, lightsBindingPoint, lightsBuffer, lightsOffset, lightsSize);

		// Set 'texture1' as the active texture at slot #0
		texture1.bind(0);

		// Drag our model's faces (triangles) in fill mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBindVertexArray(model.GetVAO());
		glDrawArrays(GL_TRIANGLES, 0, model.GetModelVertices().size());
		glBindVertexArray(0);

		// Unset 'texture1' as the active texture at slot #0
		texture1.unbind(0);

		//colorShader.setUniform("color", glm::vec3(0, 0, 0));

		//// Drag our model's faces (triangles) in line mode (wireframe)
		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		//glBindVertexArray(model.GetVAO());
		//glDrawArrays(GL_TRIANGLES, 0, model.GetModelVertices().size());
		//glBindVertexArray(0);
	}
    /*if(this->show_z_buffer){
        Rendrer_Z_buffer();
    }*/
        

}

//void Renderer::Norm(glm::vec2& v)
//{
//	v[0] = (v[0] + 1) * (Camera::screen_width / 2);
//	v[1] = (v[1] + 1) * (Camera::screen_hight / 2);
//}



//int Renderer::GetViewportWidth() const
//{
//	return viewport_width;
//}
//
//void Renderer::refreshBuff()
//{
//	delete[] color_buffer;
//	InitOpenglRendering();
//	CreateBuffers(viewport_width, viewport_height);
//}
//
//void Renderer::SetViewportWidth(int width)
//{
//	viewport_width = width;
//}
//
//int Renderer::GetViewportHeight() const
//{
//	return viewport_height;
//}
//
//void Renderer::SetViewportHeight(int height)
//{
//	viewport_height = height;
//}

void Renderer::LoadShaders()
{
	colorShader.loadShaders("vshader.glsl", "fshader.glsl");
	cubeShader.loadShaders("cubevshader.glsl", "cubefshader.glsl");
}

void Renderer::LoadTextures()
{
	//if (!texture1.loadTexture("../Data/crate.jpg" , true))
	//{
		//texture1.loadTexture("../Data/crate.jpg", true);
		texture1.loadTexture("../Data/normal.jpeg", true);
		texture1.loadTexture("../Data/diffuse.jpeg", true);

	//}
}


