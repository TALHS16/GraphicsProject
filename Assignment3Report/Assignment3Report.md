## 1. How this texture is displayed?
```
1. Create a unique ID for the texture and for the  vertex array object.
The texture ID is stored in the gl_screen_tex and the vertex array object ID is stored in the gl_screen_vtc.
2. Generate a buffer.
3. Make the vertex array object the current selection.
4. Create an ID for the buffer.
5. Create 2 const GLfloat arrays:
	a. vtc- create 2 triangles that stretch over the screen.
	b. tex - The texture that is drawn over the two triangles.
6. Bind the buffer.
7. Allocate memory to the 2 arrays (vtc and tex) and transfer the data from them to the baffer.
8. The vshader.glsl (vertex shader) and the fshader.glsl (fragment shader) are loaded and compiled, then made the current program.
9. Tells the vertex shader where to look for the vertex position data, and the data dimensions.
10. Tells the fragment shader where to look for the texture coordinates data.
11. GL_TEXTURE0 is used as the texture ID.
```
## 3. Vertex Shader
```cpp
#version 330 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 texCoords;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 orig_fragPos;


void main()
{
	orig_fragPos = vec3(vec4(pos, 1.0f));
	gl_Position = projection * view *  model * vec4(pos, 1.0f);
}
```

## 4. Fragment Shader

The fragment shader we wrote only outputs a constant color:

![color_assigment3_part4](https://user-images.githubusercontent.com/100534541/220434191-bf1e84e2-e9cb-43cd-81be-befd2ec760ce.png) `(0,0.63,0.90,1)`

```cpp
#version 330 core

in vec3 orig_fragPos;

out vec4 frag_color;

void main()
{
	frag_color = vec4(0,0.63,0.90,1);
}
```

## 5. Renderer

We updated the attributes and uniforms in the renderer as you can see in the code below:

```cpp
void Renderer::Render(Scene& scene)
{
	scene.GetActiveCamera().refresh_transform_matrix();
	//draw all the models in the scene
    for(int i=0; i<scene.GetModelCount();i++){

		MeshModel model = scene.GetModel(i);
		model.refresh_transform_matrix();

		colorShader.use();

		colorShader.setUniform("model", model.get_world_transform_matrix() * model.get_local_transform_matrix());
		colorShader.setUniform("view", scene.GetActiveCamera().GetViewTransformation());
		colorShader.setUniform("projection", scene.GetActiveCamera().GetProjectionTransformation());
		colorShader.setUniform("material.textureMap", 0);

		// Set 'texture1' as the active texture at slot #0
		texture1.bind(0);

		// Drag our model's faces (triangles) in fill mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBindVertexArray(model.GetVAO());
		glDrawArrays(GL_TRIANGLES, 0, model.GetModelVertices().size());
		glBindVertexArray(0);

		// Unset 'texture1' as the active texture at slot #0
		texture1.unbind(0);

	}
}
```
and now we can view the mesh  in our application. Here is some screenshots:

![randerer_example1](https://user-images.githubusercontent.com/100534541/220434252-e583013e-923b-4840-bfbd-929acb5e2a24.png)
![randerer_example2](https://user-images.githubusercontent.com/100534541/220434275-ea93571e-768f-4d0a-a9f2-225fa32a7568.png)

## 6. Phong Shading
We implemented the phong shading in the fragment shader. To see the implemention clearly from several angels we made a video as display below:




https://user-images.githubusercontent.com/100534541/220435470-7804692d-a59d-4542-9fc9-11304628a84c.mp4



## 7. Texture Mapping
### Plane Texture Mapping
#### cow
![assgment7_planar](https://user-images.githubusercontent.com/100534541/220434391-5c3d4f0b-b9e8-4334-bda8-b0f9bbd6b010.png)

#### spongebob
![assgment7_spongebob_planar](https://user-images.githubusercontent.com/100534541/220434403-93091d6b-f8b8-4edd-9c38-c170f2faf49f.png)

#### teapot
![assgment7_teapot_planar](https://user-images.githubusercontent.com/100534541/220434423-0c27825b-4cab-43c2-9a06-3c6b67868614.png)

### cube
![assgment7_crate_planar](https://user-images.githubusercontent.com/100534541/220434447-b27ce6b4-4a3a-4aea-b4be-db05e399410d.png)

### Cylinder Texture Mapping
#### cow
![assgment7_cylinder](https://user-images.githubusercontent.com/100534541/220434473-fdc8516c-4c54-4726-83b7-97a4cf0a86dc.png)

#### spongebob
![assgment7_spongebob_cylinder](https://user-images.githubusercontent.com/100534541/220434490-10081c4c-f86e-481d-bc3b-751c8a34d712.png)

#### teapot
![assgment7_teapot_cylinder](https://user-images.githubusercontent.com/100534541/220434528-a7becf60-3b94-4b90-98b0-e07128b8254c.png)

### cube
![assgment7_crate_cylinder](https://user-images.githubusercontent.com/100534541/220434556-66e138e4-3c0d-4972-b1cc-908c34b606c3.png)

### Sphere Texture Mapping
#### cow
![assgment7_sphere](https://user-images.githubusercontent.com/100534541/220434605-0497ac4a-9173-4e92-bd31-5e89ee8e0dea.png)

#### spongebob
![assgment7_spongebob_sphere](https://user-images.githubusercontent.com/100534541/220434630-5f2aa7ef-c71b-4277-bdc8-00a783aa3ad0.png)

#### teapot
![assgment7_teapot_sphere](https://user-images.githubusercontent.com/100534541/220434649-8d35ca86-9841-4dcd-ae94-adf6b889fb73.png)

### cube
![assgment7_crate_sphere](https://user-images.githubusercontent.com/100534541/220434674-be9959d3-ce5b-44c5-939b-2ec8fee1935f.png)

## 8. Normal mapping, Environment mapping & Toon shading (without the silhouette)
The fshader code:
```cpp
#version 330 core
# define MAX_LIGHTS 1000

struct Material
{
    sampler2D textureMap;
    vec3 ambient;
	vec3 diffuse;
	vec3 specular;
    float shine_factor;
};

struct Light {
        vec3 position;
        vec3 ambient;
	    vec3 diffuse;
	    vec3 specular;
};

uniform Light lights[100];
uniform Material material;
uniform int texture_type;
uniform int mapping_type;
uniform sampler2D diffuseMap;
uniform sampler2D normalMap;
uniform samplerCube myCube;

uniform int lights_size;

uniform vec4 camera_pos;

in vec3 orig_fragPos;
in vec3 fragPos;
in vec3 fragNormal;
in vec2 fragTexCoords;
in mat3 TBN;

out vec4 frag_color;

void main()
{
    vec3 point_normal = normalize(fragNormal);
    vec3 v = normalize(camera_pos.xyz - fragPos);
    vec3 reflectedVector = reflect(v,point_normal);
    vec3 color = vec3(0.0f, 0.0f, 0.0f);
    vec2 TextureColorType;
    vec3 textureDiffuse = vec3(texture(diffuseMap, fragTexCoords));
   
    // plane
    if(texture_type ==  0)
        TextureColorType = fragTexCoords;
    // cylinder
    else if(texture_type ==  1)
        TextureColorType = vec2(atan(fragPos.z, fragPos.x) / (2.0 * 3.14159) + 0.5, fragTexCoords.y);
    // sphere
    else if(texture_type ==  2) 
        TextureColorType = vec2(0.5 + atan(orig_fragPos.z, orig_fragPos.x) / (2.0 * 3.14159), 0.5 - atan(-orig_fragPos.y, length(orig_fragPos.xz)) / 3.14159);
    
    vec4 texture_color = texture(material.textureMap, TextureColorType);
    // if it's normal mapping
    if(mapping_type == 0)
    {
         point_normal = normalize(texture(normalMap, fragTexCoords)).xyz;
        point_normal = normalize(TBN * (2.0 * point_normal - 1.0));

    }
    for (int i = 0; i < lights_size; i++)
    {
        vec3 l = normalize(lights[i].position - fragPos);
        vec3 r = normalize(2.0f * point_normal * dot(l, point_normal) - l);

        vec3 ambient = material.ambient * lights[i].ambient;
        vec3 diffuse = material.diffuse * lights[i].diffuse * max(dot(point_normal, l), 0.0f);
        vec3 specular = material.specular * lights[i].specular * pow(max(dot(v, r), 0.0f), material.shine_factor);
        color += (ambient + diffuse + specular);
        vec3 mapping_color = vec3(1.0f, 1.0f, 1.0f);
        // if it's toon shading
        if (mapping_type == 2)
        {
	        float intensity  = dot(l, fragNormal);
	        if (intensity > 0.95)
		        mapping_color = vec3(1.0,0.5,0.5);
	        else if (intensity > 0.5)
		        mapping_color = vec3(0.6,0.3,0.3);
	        else if (intensity > 0.25)
		        mapping_color = vec3(0.4,0.2,0.2);
	        else
		        mapping_color = vec3(0.2,0.1,0.1);
            color = color * mapping_color;
        }

    }
    // if it's normal mapping
    if(mapping_type == 0)
    {
         frag_color = vec4(color * textureDiffuse, 1.0f);
    }
    // if it's environment mapping
    else if(mapping_type == 1)
    {
         frag_color = vec4(texture(myCube, reflectedVector).rgb, 1.0);
    }
    else
    {
        frag_color = vec4(texture_color.rgb * color , 1.0f);
    }
    
    
}
```
### Normal mapping
![assigment3_teapot_normal_mapping_8](https://user-images.githubusercontent.com/100534541/221507151-e5edb8de-76bd-45c1-b973-4d62cb1ecd1a.png)

### Environment mapping
![assgment8_envie](https://user-images.githubusercontent.com/100534541/221508637-72f21522-3946-4100-b3eb-21ef3f2e884b.png)

### Toon shading (without the silhouette)
![assigment3_horse_toon_shading_8](https://user-images.githubusercontent.com/100534541/220434711-1bff023e-ac66-4629-85a3-9c0faf46f077.png)
