#version 130

in vec2 texCoord;
out vec4 outColor;

uniform sampler2D tex;
uniform float min, max;

vec3 gradient[5] = vec3[](
    vec3(0.290, 0.611, 0.058),
    vec3(1.000, 0.980, 0.407),
    vec3(1.000, 0.701, 0.149),
    vec3(0.573, 0.392, 0.117),
    vec3(1.000, 1.000, 1.000)
);

vec3 geo(float val)
{
    int bin = int(clamp(val * 4.0, 0.0, 3.0));
    float ratio = val * 4.0 - bin;
    return mix(gradient[bin], gradient[bin + 1], ratio);
}

void main(void)
{
    float v = texture(tex, texCoord).r;
    outColor = vec4(geo(v / -5000.0), (v > 0) ? 0.0 : 1.0);
}
