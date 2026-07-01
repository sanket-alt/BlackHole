#define M_PI 3.1415926535897932384626433832795
#define R_SQRT_2 0.7071067811865475
#define DEG_TO_RAD (M_PI/180.0)
#define SQ(x) ((x)*(x))

#define ROT_Y(a) mat3(cos(a), 0, sin(a), 0, 1, 0, -sin(a), 0, cos(a))

const float BLACK_BODY_TEXTURE_COORD = 1.0;
const float SINGLE_WAVELENGTH_TEXTURE_COORD = 0.5;
const float TEMPERATURE_LOOKUP_RATIO_TEXTURE_COORD = 0.0;
const float SPECTRUM_TEX_TEMPERATURE_RANGE = 65504.0;
const float SPECTRUM_TEX_WAVELENGTH_RANGE = 2048.0;
const float SPECTRUM_TEX_RATIO_RANGE = 6.48053329012;

#define BLACK_BODY_COLOR(t) texture2D(spectrum_texture, vec2((t) / SPECTRUM_TEX_TEMPERATURE_RANGE, BLACK_BODY_TEXTURE_COORD))
#define SINGLE_WAVELENGTH_COLOR(lambda) texture2D(spectrum_texture, vec2((lambda) / SPECTRUM_TEX_WAVELENGTH_RANGE, SINGLE_WAVELENGTH_TEXTURE_COORD))
#define TEMPERATURE_LOOKUP(ratio) (texture2D(spectrum_texture, vec2((ratio) / SPECTRUM_TEX_RATIO_RANGE, TEMPERATURE_LOOKUP_RATIO_TEXTURE_COORD)).r * SPECTRUM_TEX_TEMPERATURE_RANGE)

uniform vec2 resolution;
uniform float time;
uniform vec3 cam_pos;
uniform vec3 cam_x;
uniform vec3 cam_y;
uniform vec3 cam_z;
uniform vec3 cam_vel;

uniform float planet_distance, planet_radius;
uniform float noise_scale;
uniform float noise_speed;

uniform float magnetic_strength;
uniform float magnetic_dipole_tilt;
uniform float field_lines_density;
uniform vec3 field_color;
uniform int show_plasma_effects;

uniform bool hawking_radiation_enabled;
uniform int hawking_mode;
uniform float hawking_temperature;
uniform float bh_radius;

uniform sampler2D galaxy_texture, star_texture, planet_texture, spectrum_texture;

const int NSTEPS = 256;
const float MAX_REVOLUTIONS = 2.0;

const float ACCRETION_MIN_R_FACTOR = 1.3;
const float ACCRETION_WIDTH_FACTOR = 6.0;
const float ACCRETION_BRIGHTNESS = 1.0;
const float ACCRETION_TEMPERATURE = 12000.0;

const float STAR_MIN_TEMPERATURE = 4000.0;
const float STAR_MAX_TEMPERATURE = 15000.0;
const float STAR_BRIGHTNESS = 1.0;
const float GALAXY_BRIGHTNESS = 0.4;
const float PLANET_AMBIENT = 0.1;
const float PLANET_LIGHTNESS = 1.5;

mat3 BG_COORDS = ROT_Y(45.0 * DEG_TO_RAD);
const float PLANET_AXIAL_TILT = 30.0 * DEG_TO_RAD;
mat3 PLANET_COORDS = ROT_Y(PLANET_AXIAL_TILT);
const float FOV_ANGLE_DEG = 90.0;
float FOV_MULT = 1.0 / tan(DEG_TO_RAD * FOV_ANGLE_DEG*0.5);

vec3 hash33(vec3 p) {
    p = fract(p * vec3(443.897, 441.423, 437.195));
    p += dot(p, p.yzx + 19.19);
    return fract((p.xxy + p.yxx)*p.zyx);
}

float rand(vec2 n) {
    return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float noise(vec2 p) {
    vec2 ip = floor(p);
    vec2 u = fract(p);
    u = u*u*(3.0-2.0*u);
    float res = mix(
        mix(rand(ip), rand(ip + vec2(1.0, 0.0)), u.x),
        mix(rand(ip + vec2(0.0, 1.0)), rand(ip + vec2(1.0, 1.0)), u.x),
        u.y);
    return res*res;
}

float fbm(vec2 p) {
    float f = 0.0;
    mat2 m = mat2(1.6, 1.2, -1.2, 1.6);
    f += 0.5000 * noise(p); p = m * p;
    f += 0.2500 * noise(p); p = m * p;
    f += 0.1250 * noise(p); p = m * p;
    f += 0.0625 * noise(p);
    return f / 0.9375;
}

vec4 get_hawking_radiation_color(vec3 pos_at_horizon, float u) {
    if (!hawking_radiation_enabled) return vec4(0.0);
    float horizon_proximity = 1.0 - smoothstep(1.0/bh_radius, 1.0/bh_radius + 0.5, u);
    if (horizon_proximity < 0.01) return vec4(0.0);
    vec4 color = vec4(0.0);
    vec3 p_hash = hash33(pos_at_horizon * 50.0);
    if (hawking_mode == 0) {
        float time_phase = fract(time * 0.5 + p_hash.x);
        if (time_phase > 0.95 && p_hash.y > 0.5) {
            float lifetime = smoothstep(0.95, 0.96, time_phase) * smoothstep(1.0, 0.99, time_phase);
            if (lifetime > 0.01) {
                vec3 particle1_pos = pos_at_horizon + normalize(p_hash - 0.5) * 0.05;
                vec3 particle2_pos = pos_at_horizon - normalize(p_hash - 0.5) * 0.05;
                float d1 = length(pos_at_horizon - particle1_pos) - 0.02;
                float d2 = length(pos_at_horizon - particle2_pos) - 0.02;
                float escape_factor = smoothstep(1.0/bh_radius, 1.0/bh_radius + 0.02, u);
                color += vec4(0.5, 0.8, 1.0, 1.0) * (1.0 - smoothstep(0.0, 0.01, d1)) * lifetime * escape_factor;
                color += vec4(1.0, 0.5, 0.5, 1.0) * (1.0 - smoothstep(0.0, 0.01, d2)) * lifetime * (1.0 - escape_factor);
            }
        }
    } else {
        float visual_temp = hawking_temperature * 1e-6;
        float particle_speed = 5.0 + visual_temp * 0.1;
        float particle_brightness = 0.5 + visual_temp * 0.01;
        float particle_size = max(0.005, 0.05 - visual_temp * 0.001);
        float time_phase = fract(time * particle_speed + p_hash.x);
        vec3 particle_pos = pos_at_horizon + normalize(p_hash - 0.5) * time_phase * 2.0;
        float d = length(pos_at_horizon - particle_pos) - particle_size;
        float particle_glow = (1.0 - smoothstep(0.0, 0.1, d));
        if (particle_glow > 0.01 && p_hash.y > 0.8) {
            vec4 particle_color = BLACK_BODY_COLOR(visual_temp);
            if (hawking_mode == 2) {
                particle_brightness *= 5.0;
                float jet_focus = 1.0 - abs(dot(normalize(pos_at_horizon), vec3(p_hash.x-0.5, p_hash.y-0.5, 1.0)));
                particle_brightness *= pow(jet_focus, 2.0);
            }
            color += particle_color * particle_glow * particle_brightness * horizon_proximity;
        }
    }
    return color;
}

vec2 sphere_map(vec3 p) {
    return vec2(atan(p.x,p.y)/M_PI*0.5+0.5, asin(p.z)/M_PI+0.5);
}

float smooth_step_custom(float x, float threshold, float steepness) {
    return 1.0 / (1.0 + exp(-(x-threshold)*steepness));
}

void main() {
    vec2 uv = gl_FragCoord.xy / resolution;
    gl_FragColor = vec4(uv.x, uv.y, 0.5, 1.0);
}
