#define M_PI 3.1415926535897932384626433832795
#define R_SQRT_2 0.7071067811865475
#define DEG_TO_RAD (M_PI/180.0)
#define SQ(x) ((x)*(x))

#define ROT_Y(a) mat3(cos(a), 0, sin(a), 0, 1, 0, -sin(a), 0, cos(a))

uniform vec2 resolution;
uniform float time;
uniform vec3 cam_pos;
uniform vec3 cam_x;
uniform vec3 cam_y;
uniform vec3 cam_z;

void main() {
    vec2 uv = gl_FragCoord.xy / resolution.xy;
    gl_FragColor = vec4(uv, 0.5, 1.0);
}