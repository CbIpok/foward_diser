#version 130

in vec2 texCoord;
out vec4 outColor;

uniform sampler2D tex;
uniform float min, max;

void main(void)
{
    float v = texture(tex, texCoord).r;
    v = (v - min) / (max - min);
    
    float sat = 1.0 - clamp(2.0 * (v - 0.5), 0.0, 1.0);
    float value = clamp(2.0 * v, 0.0, 1.0);
    
    outColor = vec4(((vec3(0.0, 0.0, 1.0) - 1.0) * sat + 1.0) * value, 0.75);
}
