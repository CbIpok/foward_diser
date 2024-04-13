#version 130

in vec2 texCoord;
out vec4 outColor;

uniform sampler2D tex;
uniform float min, max;

void main(void)
{
    //outColor = vec4(inColor, 1.0);
    float v = texture(tex, texCoord).r;
    v = (v - min) / (max - min);
    outColor = vec4(v, v, v, 1.0);
}
