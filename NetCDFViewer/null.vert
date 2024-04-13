#version 130

uniform mat4x4 transform;

in  vec3 inPosition;
in  vec2 inTexCoord;

out vec2 texCoord;

void main(void)
{
    gl_Position = transform * vec4(inPosition, 1.0);
    texCoord = inTexCoord;
}
