#version 330 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 texCoords;
layout(location = 3) in vec3 tangent;
layout(location = 4) in vec3 bitangent;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 orig_fragPos;
out vec3 fragPos;
out vec3 fragNormal;
out vec2 fragTexCoords;
out mat3 TBN;


void main()
{
	orig_fragPos = vec3(vec4(pos, 1.0f));
    fragPos = vec3(model * vec4(pos, 1.0f));
    fragNormal = mat3(view * model) * normal;
    fragTexCoords = texCoords;

    // transform the tangent and bitangent vectors into view space
    vec3 T = normalize(vec3(view * vec4(tangent, 0.0)));
    vec3 B = normalize(vec3(view * vec4(bitangent, 0.0)));
    vec3 N = normalize(vec3(view * vec4(normal, 0.0)));

    // create a matrix that transforms normal vectors in texture space to object space
    TBN = mat3(T, B, N);

	gl_Position = projection * view *  model * vec4(pos, 1.0f);
}