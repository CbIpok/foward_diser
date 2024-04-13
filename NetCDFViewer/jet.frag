#version 130

in vec2 texCoord;
out vec4 outColor;

uniform sampler2D tex;
uniform float min, max;

float jetInterpolate(float val, float y0, float x0, float y1, float x1)
{
    return (val - x0) * (y1 - y0) / (x1 - x0) + y0;
}

float jetBase(float val)
{
    if (val <= -0.75)
        return 0.0;
    else if (val <= -0.25)
        return jetInterpolate(val, 0.0, -0.75, 1.0, -0.25);
    else if (val <= 0.25)
        return 1.0;
    else if (val <= 0.75)
        return jetInterpolate(val, 1.0, 0.25, 0.0, 0.75);
    else
        return 0.0;
}

vec3 jetColor(float value)
{
    return vec3(jetBase(value - 0.5), jetBase(value), jetBase(value + 0.5));
}

void main(void)
{
    float v = texture(tex, texCoord).r;
    v = (v - min) / (max - min);
    outColor = vec4(jetColor(v), 1.0);
}
