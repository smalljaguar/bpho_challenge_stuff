#version 330

uniform vec2 resolution;
out vec4 fragColor;

void main()
{
    vec2 uv = gl_FragCoord.xy / resolution;
    vec3 color = mix(vec3(0.5, 0.7, 1.0), vec3(1.0), uv.y);
    fragColor = vec4(color, 1.0);
}