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
    if(mapping_type == 0)
    {
         frag_color = vec4(color * textureDiffuse, 1.0f);
    }
    else if(mapping_type == 1)
    {
         //reflection mode
         frag_color = vec4(texture(myCube, reflectedVector).rgb, 1.0);
         //float ratio = 1.00 / 1.52;
         //vec3 R = refract(v, point_normal, ratio);
         //frag_color = vec4(texture(myCube, R).rgb, 1.0);
    }
    else
    {
        frag_color = vec4(texture_color.rgb * color , 1.0f);
    }
    
    
}