#version 130

in vec2 texCoord;
out vec4 outColor;

uniform sampler2D tex;
uniform float min, max;

void main(void)
{
    float v = texture(tex, texCoord).r;    
    float value = 1.0 - clamp(v / 5000.0, 0.0, 1.0);
    
    outColor = vec4(vec3(0.0, 0.5, 1.0) * value, 1.0);
}
