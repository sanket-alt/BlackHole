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

// Magnetic field uniforms
uniform float magnetic_strength;
uniform float magnetic_dipole_tilt;
uniform float field_lines_density;
uniform vec3 field_color;
uniform int show_plasma_effects;

// --- HAWKING RADIATION UNIFORMS ---
uniform bool hawking_radiation_enabled;
uniform int hawking_mode; // 0: Quantum Source, 1: Energy Spectrum, 2: Evaporation
uniform float hawking_temperature;
uniform float bh_radius; // Schwarzschild radius, now dynamic

uniform sampler2D galaxy_texture, star_texture,
    planet_texture, spectrum_texture;

// stepping parameters
const int NSTEPS = {{n_steps}};
const float MAX_REVOLUTIONS = 2.0;

const float ACCRETION_MIN_R_FACTOR = 1.3; // Factor of bh_radius - smaller disk
const float ACCRETION_WIDTH_FACTOR = 6.0; // Factor of bh_radius - more compact disk
const float ACCRETION_BRIGHTNESS = 1.0; // Reduced brightness
const float ACCRETION_TEMPERATURE = 12000.0; // Higher base temperature

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

float PLANET_RADIUS,
    PLANET_DISTANCE,
    PLANET_ORBITAL_ANG_VEL,
    PLANET_ROTATION_ANG_VEL,
    PLANET_GAMMA;

// --- START PROCEDURAL NOISE FUNCTIONS ---
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
// --- END PROCEDURAL NOISE FUNCTIONS ---

// --- HAWKING RADIATION VISUALIZATION (USER-FRIENDLY) ---
vec4 get_hawking_radiation_color(vec3 pos_at_horizon, float u) {
    if (!hawking_radiation_enabled) return vec4(0.0);

    float horizon_proximity = 1.0 - smoothstep(1.0/bh_radius, 1.0/bh_radius + 0.5, u);
    if (horizon_proximity < 0.01) return vec4(0.0);

    vec4 color = vec4(0.0);
    vec3 p_hash = hash33(pos_at_horizon * 50.0); // Use a larger scale for more distinct particles

    // --- Mode 0: Quantum Source (Particle Pair Analogy) ---
    if (hawking_mode == 0) {
        float time_phase = fract(time * 0.5 + p_hash.x);

        // Create a short-lived particle pair
        if (time_phase > 0.95 && p_hash.y > 0.5) {
            float lifetime = smoothstep(0.95, 0.96, time_phase) * smoothstep(1.0, 0.99, time_phase);
            if (lifetime > 0.01) {
                // Visualize two particles, one blue, one red
                vec3 particle1_pos = pos_at_horizon + normalize(p_hash - 0.5) * 0.05;
                vec3 particle2_pos = pos_at_horizon - normalize(p_hash - 0.5) * 0.05;

                // Simple sphere SDF to draw the particles
                float d1 = length(pos_at_horizon - particle1_pos) - 0.02;
                float d2 = length(pos_at_horizon - particle2_pos) - 0.02;

                // One particle is captured (fades inside horizon), one escapes
                float escape_factor = smoothstep(1.0/bh_radius, 1.0/bh_radius + 0.02, u);

                // Blue particle (matter)
                color += vec4(0.5, 0.8, 1.0, 1.0) * (1.0 - smoothstep(0.0, 0.01, d1)) * lifetime * escape_factor;
                // Red particle (anti-matter)
                color += vec4(1.0, 0.5, 0.5, 1.0) * (1.0 - smoothstep(0.0, 0.01, d2)) * lifetime * (1.0 - escape_factor);
            }
        }
    }
    // --- Mode 1 & 2: Energy Spectrum & Evaporation (Intuitive Particles) ---
    else {
        // Visual properties based on temperature
        float visual_temp = hawking_temperature * 1e-6;
        float particle_speed = 5.0 + visual_temp * 0.1;
        float particle_brightness = 0.5 + visual_temp * 0.01;
        float particle_size = max(0.005, 0.05 - visual_temp * 0.001);

        // Animate particles moving away from the black hole
        float time_phase = fract(time * particle_speed + p_hash.x);
        vec3 particle_pos = pos_at_horizon + normalize(p_hash - 0.5) * time_phase * 2.0;

        // Draw the particle as a sphere
        float d = length(pos_at_horizon - particle_pos) - particle_size;
        float particle_glow = (1.0 - smoothstep(0.0, 0.1, d));

        if (particle_glow > 0.01 && p_hash.y > 0.8) { // Only draw some particles
            vec4 particle_color = BLACK_BODY_COLOR(visual_temp);

            // In evaporation mode, make the effect a more intense "wind"
            if (hawking_mode == 2) {
                particle_brightness *= 5.0; // More intense
                // Create a "jet" or "wind" effect by shaping the particle distribution
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

vec3 lorentz_velocity_transformation(vec3 moving_v, vec3 frame_v) {
    float v = length(frame_v);
    if (v > 0.0) {
        vec3 v_axis = -frame_v / v;
        float gamma = 1.0/sqrt(1.0 - v*v);

        float moving_par = dot(moving_v, v_axis);
        vec3 moving_perp = moving_v - v_axis*moving_par;

        float denom = 1.0 + v*moving_par;
        return (v_axis*(moving_par+v)+moving_perp/gamma)/denom;
    }
    return moving_v;
}

vec3 contract(vec3 x, vec3 d, float mult) {
    float par = dot(x,d);
    return (x-par*d) + d*par*mult;
}

vec4 planet_intersection(vec3 old_pos, vec3 ray, float t, float dt,
        vec3 planet_pos0, float ray_doppler_factor) {

    vec4 ret = vec4(0,0,0,0);
    vec3 ray0 = ray;
    ray = ray/dt;

    vec3 planet_dir = vec3(planet_pos0.y, -planet_pos0.x, 0.0) / PLANET_DISTANCE;

    {{#light_travel_time}}
    float planet_ang1 = (t-dt) * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos1 = vec3(cos(planet_ang1), sin(planet_ang1), 0)*PLANET_DISTANCE;
    vec3 planet_vel = (planet_pos1-planet_pos0)/dt;

    ray = ray - planet_vel;
    {{/light_travel_time}}
    {{^light_travel_time}}
    vec3 planet_vel = planet_dir * PLANET_ORBITAL_ANG_VEL * PLANET_DISTANCE;
    {{/light_travel_time}}

    vec3 d = old_pos - planet_pos0;

    {{#lorentz_contraction}}
    ray = contract(ray, planet_dir, PLANET_GAMMA);
    d = contract(d, planet_dir, PLANET_GAMMA);
    {{/lorentz_contraction}}

    float dotp = dot(d,ray);
    float c_coeff = dot(d,d) - SQ(PLANET_RADIUS);
    float ray2 = dot(ray, ray);
    float discr = dotp*dotp - ray2*c_coeff;

    if (discr < 0.0) return ret;
    float isec_t = (-dotp - sqrt(discr)) / ray2;

    float MIN_ISEC_DT = 0.0;
    {{#lorentz_contraction}}
    MIN_ISEC_DT = -dt;
    {{/lorentz_contraction}}

    if (isec_t < MIN_ISEC_DT || isec_t > dt) return ret;

    vec3 surface_point = (d + isec_t*ray) / PLANET_RADIUS;

    isec_t = isec_t/dt;

    vec3 light_dir = planet_pos0;
    float rot_phase = t;

    {{#light_travel_time}}
    light_dir += planet_vel*isec_t*dt;
    rot_phase -= isec_t*dt;
    {{/light_travel_time}}

    rot_phase = rot_phase * PLANET_ROTATION_ANG_VEL*0.5/M_PI;
    light_dir = light_dir / PLANET_DISTANCE;

    {{#light_travel_time}}
    light_dir = light_dir - planet_vel;
    {{/light_travel_time}}

    vec3 surface_normal = surface_point;
    {{#lorentz_contraction}}
    light_dir = contract(light_dir, planet_dir, PLANET_GAMMA);
    {{/lorentz_contraction}}
    light_dir = normalize(light_dir);

    vec2 tex_coord = sphere_map(surface_point * PLANET_COORDS);
    tex_coord.x = mod(tex_coord.x + rot_phase, 1.0);

    float diffuse = max(0.0, dot(surface_normal, -light_dir));
    float lightness = ((1.0-PLANET_AMBIENT)*diffuse + PLANET_AMBIENT) *
        PLANET_LIGHTNESS;

    float light_temperature = ACCRETION_TEMPERATURE;
    {{#doppler_shift}}
    float doppler_factor = SQ(PLANET_GAMMA) *
        (1.0 + dot(planet_vel, light_dir)) *
        (1.0 - dot(planet_vel, normalize(ray)));
    light_temperature /= doppler_factor * ray_doppler_factor;
    {{/doppler_shift}}

    vec4 light_color = BLACK_BODY_COLOR(light_temperature);
    ret = texture2D(planet_texture, tex_coord) * lightness * light_color;
    if (isec_t < 0.0) isec_t = 0.5;
    ret.w = isec_t;

    return ret;
}

vec4 galaxy_color(vec2 tex_coord, float doppler_factor) {

    vec4 color = texture2D(galaxy_texture, tex_coord);
    {{^observerMotion}}
    return color;
    {{/observerMotion}}

    {{#observerMotion}}
    vec4 ret = vec4(0.0,0.0,0.0,0.0);
    float red = max(0.0, color.r - color.g);

    const float H_ALPHA_RATIO = 0.1;
    const float TEMPERATURE_BIAS = 0.95;

    color.r -= red*H_ALPHA_RATIO;

    float i1 = max(color.r, max(color.g, color.b));
    float ratio = (color.g+color.b) / color.r;

    if (i1 > 0.0 && color.r > 0.0) {

        float temperature = TEMPERATURE_LOOKUP(ratio) * TEMPERATURE_BIAS;
        color = BLACK_BODY_COLOR(temperature);

        float i0 = max(color.r, max(color.g, color.b));
        if (i0 > 0.0) {
            temperature /= doppler_factor;
            ret = BLACK_BODY_COLOR(temperature) * max(i1/i0,0.0);
        }
    }

    ret += SINGLE_WAVELENGTH_COLOR(656.28 * doppler_factor) * red / 0.214 * H_ALPHA_RATIO;

    return ret;
    {{/observerMotion}}
}

vec3 get_camera_ray_dir(vec2 screen_pos) {
    return cam_x * screen_pos.x + cam_y * screen_pos.y + cam_z;
}

// Physics-based magnetic dipole field calculation
vec3 magneticField(vec3 pos) {
    float r = length(pos);
    if (r < 1.05) return vec3(0.0); // Inside event horizon

    // Convert to spherical coordinates centered at black hole
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);

    // Magnetic dipole moment (aligned with z-axis, then tilted)
    float r3 = r * r * r;
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);

    // Classical magnetic dipole field in spherical coordinates
    float B_r = (2.0 * magnetic_strength * cos_theta) / r3;
    float B_theta = (magnetic_strength * sin_theta) / r3;
    float B_phi = 0.0; // Axisymmetric dipole

    // Convert to Cartesian coordinates
    vec3 B_field;
    B_field.x = B_r * sin_theta * cos(phi) + B_theta * cos_theta * cos(phi);
    B_field.y = B_r * sin_theta * sin(phi) + B_theta * cos_theta * sin(phi);
    B_field.z = B_r * cos_theta - B_theta * sin_theta;

    // Apply dipole tilt rotation
    float tilt_cos = cos(magnetic_dipole_tilt);
    float tilt_sin = sin(magnetic_dipole_tilt);

    // Rotate around y-axis
    vec3 tilted_field;
    tilted_field.x = B_field.x * tilt_cos + B_field.z * tilt_sin;
    tilted_field.y = B_field.y;
    tilted_field.z = -B_field.x * tilt_sin + B_field.z * tilt_cos;

    return tilted_field;
}

// Calculate distance to nearest magnetic field line with animation
float distanceToFieldLine(vec3 pos) {
    float r = length(pos);
    if (r < 1.1) return 1000.0;

    vec3 B = magneticField(pos);
    float B_mag = length(B);
    if (B_mag < 0.001) return 1000.0;

    vec3 B_normalized = B / B_mag;

    // Calculate field line parameter using magnetic flux conservation
    // For a dipole: flux ~ sin^2(theta) / r
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);
    float sin_theta = sin(theta);
    float flux_param = sin_theta * sin_theta / r;

    // Animate field lines by rotating them around the magnetic axis
    float animation_speed = 0.3;
    float animated_phi = phi + time * animation_speed;

    // Add magnetic field line twist animation (simulates field line reconnection)
    float twist_factor = sin(time * 0.5 + r * 1.5) * 0.1;
    animated_phi += twist_factor;

    // Reduce line density significantly for better visibility
    float reduced_density = field_lines_density * 0.3; // Reduce by 70%
    float line_spacing = 1.2 / reduced_density;
    float nearest_line = floor(flux_param / line_spacing) * line_spacing;
    float dist_to_line = abs(flux_param - nearest_line);

    // Add animated polar jet field lines with pulsing
    float polar_dist = abs(cos(theta));
    if (polar_dist > 0.8) {
        // Animate jet rotation
        float jet_animation = time * 0.2;
        float jet_line = floor((animated_phi + jet_animation) * 6.0 / (2.0 * M_PI)) * 2.0 * M_PI / 6.0;
        float angular_dist = abs(mod(animated_phi + M_PI, 2.0 * M_PI) - M_PI - jet_line);

        // Add pulsing effect to jets
        float jet_pulse = 0.8 + 0.2 * sin(time * 2.0);
        dist_to_line = min(dist_to_line, angular_dist * r * 0.08 * jet_pulse);
    }

    // Add equatorial current sheet animation
    float equatorial_dist = abs(sin(theta));
    if (equatorial_dist < 0.2) {
        float current_wave = sin(animated_phi * 4.0 + time * 1.5) * 0.05;
        dist_to_line = min(dist_to_line, (equatorial_dist + current_wave) * 20.0);
    }

    return dist_to_line * 60.0; // Slightly increased scale for better visibility
}

// Field line visualization with proper physics and animation
float fieldLineIntensity(vec3 pos, vec3 ray_dir) {
    float r = length(pos);
    if (r < 1.1 || r > 20.0) return 0.0;

    vec3 B = magneticField(pos);
    float B_mag = length(B);
    if (B_mag < 0.001) return 0.0;

    float dist_to_line = distanceToFieldLine(pos);
    float line_width = 0.06 + 0.03 * (r - 1.0); // Slightly thicker lines for better visibility

    // Smooth field line visualization with animation
    float line_intensity = 1.0 - smoothstep(0.0, line_width, dist_to_line);

    // Add flowing energy animation along field lines
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);

    // Create flowing effect along field lines
    float flow_speed = 2.0;
    float flow_pattern = sin(r * 3.0 - time * flow_speed) * 0.3 + 0.7;
    line_intensity *= flow_pattern;

    // Enhanced polar jets with pulsing animation
    float polar_enhancement = 1.0 + 3.0 * exp(-12.0 * sin(theta) * sin(theta));

    // Add jet pulsing
    float jet_pulse = 0.7 + 0.3 * sin(time * 1.8);
    if (abs(cos(theta)) > 0.7) {
        polar_enhancement *= jet_pulse;
    }

    // Field strength falloff with animation
    float field_strength = B_mag * magnetic_strength;
    float strength_factor = smoothstep(0.0, 0.12, field_strength);

    // Add magnetic activity fluctuations
    float activity = 0.8 + 0.2 * sin(time * 0.8 + r);
    strength_factor *= activity;

    // Distance-based intensity falloff
    float distance_factor = exp(-0.04 * (r - 1.0));

    // Add occasional field line brightening (magnetic reconnection events)
    float reconnection_event = smoothstep(0.95, 1.0, sin(time * 0.3 + phi * 2.0)) * 2.0;

    return line_intensity * polar_enhancement * strength_factor * distance_factor * (1.0 + reconnection_event);
}

// Enhanced plasma effects with physics-based heating
vec3 plasmaGlow(vec3 pos, float field_strength) {
    if (show_plasma_effects == 0) return vec3(0.0);

    float r = length(pos);
    vec3 B = magneticField(pos);
    float B_mag = length(B);

    // Plasma heating from magnetic reconnection and field compression
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float reconnection_heating = exp(-5.0 * abs(sin(theta))); // Heating at equatorial current sheet
    float compression_heating = B_mag * 0.5; // Heating proportional to field strength

    float total_heating = (reconnection_heating + compression_heating) * field_strength;
    float plasma_temp = 5000.0 + total_heating * 10000.0; // Temperature in Kelvin

    // Convert temperature to color (simplified blackbody)
    vec3 plasma_color;
    if (plasma_temp < 3000.0) {
        plasma_color = vec3(1.0, 0.3, 0.1); // Red hot
    } else if (plasma_temp < 6000.0) {
        plasma_color = vec3(1.0, 0.8, 0.4); // Orange-yellow
    } else if (plasma_temp < 10000.0) {
        plasma_color = vec3(1.0, 1.0, 0.8); // White hot
    } else {
        plasma_color = vec3(0.8, 0.9, 1.0); // Blue hot
    }

    float intensity = total_heating * exp(-0.1 * r);
    return plasma_color * intensity;
}

// Get different colors for different magnetic field components with animation
// Improved magnetic field color visualization
vec3 getMagneticFieldColor(vec3 pos, float field_strength) {
    float r = length(pos);
    vec3 B = magneticField(pos);
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);

    // Color definitions for different field line types
    vec3 dipole_color = vec3(0.0, 0.8, 1.0); // Cyan for regular dipole lines
    vec3 twisted_color = vec3(1.0, 0.2, 0.8); // Magenta for twisted/complex lines
    vec3 polar_jet_color = vec3(0.0, 1.0, 0.2); // Green for polar jets

    // Animate the twist factor for dynamic field evolution
    float twist_animation = time * 0.15 + r * 1.5;
    float twist_factor = sin(twist_animation) * 0.5 + 0.5;

    // Add magnetic activity zones that change over time
    float activity_zone = sin(time * 0.4 + phi * 3.0) * 0.5 + 0.5;
    vec3 final_color;

    // Polar jets (strong vertical field near poles) with pulsing
    if (abs(cos(theta)) > 0.65) {
        float jet_intensity = 0.8 + 0.2 * sin(time * 2.5);
        final_color = polar_jet_color * jet_intensity;
    }
    // Twisted field lines with animated regions
    else if ((twist_factor > 0.55 && abs(sin(theta)) > 0.25) || activity_zone > 0.75) {
        float twist_intensity = 0.7 + 0.3 * sin(time * 1.2 + r * 0.8);
        final_color = twisted_color * twist_intensity;
    }
    // Regular dipole field lines with flowing energy
    else {
        float flow_intensity = 0.6 + 0.4 * sin(time * 1.0 - r * 2.0);
        final_color = dipole_color * flow_intensity;
    }

    // Add occasional brightness flares (magnetic reconnection)
    float flare = smoothstep(0.98, 1.0, sin(time * 0.25 + phi * 4.0 + theta * 2.0));
    final_color += vec3(0.8, 0.8, 0.2) * flare * 0.5;

    return final_color * field_strength;
}

// Enhanced polar jets with advanced physics-based features
// New and improved polar jet emission with enhanced physics
vec3 polarJetEmission(vec3 pos, vec3 ray_dir) {
    float r = length(pos);
    // Jets only exist outside the event horizon and within a certain range
    if (r < 1.2 || r > 80.0) return vec3(0.0);
    
    // Calculate spherical coordinates for jet axis alignment
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);

    float polar_factor = abs(cos(theta));
    
    // Cylindrical coordinates for the jet itself
    float z_jet = abs(pos.z) - 1.2;
    float rho_jet = sqrt(pos.x * pos.x + pos.y * pos.y);

    // === PHYSICS-BASED JET COLLIMATION ===
    // More realistic magnetic collimation profile: r(z) ∝ z^0.4
    // This is a standard model for astrophysical jets
    float collimation_exponent = 0.4;
    float jet_base_radius = 0.15; // Initial jet radius at launch
    float expected_jet_radius = jet_base_radius * pow(max(z_jet / jet_base_radius, 1.0), collimation_exponent);
    
    // Magnetic pinching factor - stronger near base to collimate the jet
    float magnetic_pinch_strength = 2.0 / (1.0 + 0.1 * z_jet);
    expected_jet_radius /= magnetic_pinch_strength;
    
    // Determine if the ray is inside the jet cone
    if (rho_jet > expected_jet_radius * 2.0) return vec3(0.0);

    // === CORE-SHEATH STRUCTURE ===
    // Model jets with a fast-moving core and a slower, turbulent sheath
    float core_radius = expected_jet_radius * 0.3;
    float sheath_radius = expected_jet_radius;
    
    bool in_core = rho_jet < core_radius;
    bool in_sheath = rho_jet >= core_radius && rho_jet < sheath_radius;
    
    if (!in_core && !in_sheath) return vec3(0.0);

    // === VELOCITY STRUCTURE ===
    float core_velocity = 0.95; // Highly relativistic core velocity
    float sheath_velocity = 0.6; // Slower sheath velocity
    
    float jet_velocity = in_core ? core_velocity : mix(sheath_velocity, core_velocity, smoothstep(sheath_radius, core_radius, rho_jet));

    // Improved helical magnetic field visualization
    float helical_pitch = 3.0 + 1.0 * sin(time * 0.5);
    float helix_phase = phi + helical_pitch * z_jet / expected_jet_radius + time * 1.2;
    vec3 jet_direction = normalize(vec3(0.0, 0.0, pos.z > 0.0 ? 1.0 : -1.0));
    vec3 jet_bulk_velocity = jet_direction * jet_velocity;
    float helical_amplitude = in_core ? 0.05 : 0.15; // Stronger helical motion in the sheath
    vec3 helical_velocity = vec3(
        cos(helix_phase) * helical_amplitude * jet_velocity,
        sin(helix_phase) * helical_amplitude * jet_velocity,
        0.0
    );
    vec3 total_jet_velocity = jet_bulk_velocity + helical_velocity;

    // === RELATIVISTIC EFFECTS ===
    float gamma_jet = 1.0 / sqrt(1.0 - dot(total_jet_velocity, total_jet_velocity));
    float doppler_factor = gamma_jet * (1.0 + dot(normalize(ray_dir), total_jet_velocity));
    float beaming_factor = pow(max(doppler_factor, 0.1), 2.5);

    // === KNOTS AND INSTABILITIES ===
    // Create moving knots (overdense regions)
    float knot_speed = 0.8 * jet_velocity;
    float knot_spacing = 5.0;
    float knot_position = mod(z_jet - time * knot_speed * 20.0, knot_spacing);
    float knot1 = exp(-8.0 * pow(knot_position - knot_spacing * 0.2, 2.0));
    float knot2 = exp(-12.0 * pow(mod(z_jet - time * knot_speed * 15.0 + 2.0, knot_spacing) - knot_spacing * 0.3, 2.0));
    float knot3 = exp(-10.0 * pow(mod(z_jet - time * knot_speed * 25.0 + 4.0, knot_spacing) - knot_spacing * 0.5, 2.0));
    float knot_enhancement = 1.0 + 1.5 * (knot1 + knot2 + knot3);
    knot_enhancement *= 0.8 + 0.2 * sin(time * 3.0 + z_jet * 0.5);

    // Kelvin-Helmholtz instabilities at core-sheath boundary
    vec3 turb_coord = pos * 1.2 + vec3(time * 3.0, time * 2.5, time * 1.8);
    float kh_instability = 0.0;
    if (abs(rho_jet - core_radius) < core_radius * 0.5) {
        for (int i = 0; i < 3; i++) {
            float freq = pow(2.0, float(i)) * 0.8;
            float amp = 0.4 / pow(2.0, float(i));
            kh_instability += amp * (fbm(turb_coord.xy * freq + turb_coord.z * 0.3) - 0.5);
        }
        float boundary_factor = exp(-5.0 * abs(rho_jet - core_radius) / core_radius);
        kh_instability *= boundary_factor;
    }
    
    // Rayleigh-Taylor instabilities (density fingers)
    float rt_instability = 0.0;
    if (in_sheath) {
        rt_instability = 0.3 * sin(z_jet * 2.0 + time * 4.0 + phi * 6.0) * exp(-0.5 * rho_jet / expected_jet_radius);
    }

    // === SYNCHROTRON EMISSION GRADIENT ===
    float magnetic_field_strength = magnetic_pinch_strength * (in_core ? 3.0 : 1.5);
    float cooling_length = 15.0;
    float cooling_factor = exp(-z_jet / cooling_length);
    float base_temp_core = 150000.0;
    float base_temp_sheath = 80000.0;
    
    float base_temperature = in_core ? base_temp_core : base_temp_sheath;
    base_temperature *= cooling_factor;
    
    float magnetic_heating = 1.0 + magnetic_field_strength * 0.3;
    float shock_heating = 1.0 + 0.5 * (abs(kh_instability) + abs(rt_instability));
    float knot_heating = knot_enhancement * 0.8 + 0.2;
    float plasma_temp = base_temperature * magnetic_heating * shock_heating * knot_heating;
    plasma_temp /= doppler_factor;

    // === SYNCHROTRON COLOR GRADIENT ===
    vec3 synchrotron_color;
    if (z_jet < 5.0) {
        synchrotron_color = mix(vec3(1.0, 1.0, 1.0), vec3(0.8, 0.9, 1.0), z_jet / 5.0);
    }
    else if (z_jet < 15.0) {
        float t = (z_jet - 5.0) / 10.0;
        synchrotron_color = mix(vec3(0.8, 0.9, 1.0), vec3(0.5, 1.0, 0.7), t);
    }
    else {
        float t = min((z_jet - 15.0) / 20.0, 1.0);
        synchrotron_color = mix(vec3(0.5, 1.0, 0.7), vec3(1.0, 0.6, 0.3), t);
    }
    
    if (in_core) {
        synchrotron_color *= 1.4;
    }

    // === HELICAL MAGNETIC FIELD VISUALIZATION ===
    float helix_vis_strength = 0.0;
    float helix_line_dist = abs(sin(helix_phase * 8.0)) * abs(cos(helix_phase * 4.0));
    helix_vis_strength = exp(-20.0 * helix_line_dist) * 0.3;
    
    vec3 poloidal_color = vec3(1.0, 0.7, 0.3);
    vec3 toroidal_color = vec3(0.3, 0.7, 1.0);
    
    float field_mix = 0.5 + 0.5 * sin(helix_phase * 2.0);
    vec3 helix_color = mix(poloidal_color, toroidal_color, field_mix);

    // === TOTAL EMISSION CALCULATION ===
    float distance_falloff = 1.0 / (1.0 + 0.005 * z_jet * z_jet);
    float density_factor = in_core ? 2.0 : 1.0;
    
    float confinement = exp(-2.0 * pow(rho_jet / expected_jet_radius, 2.0));
    float base_intensity = density_factor * confinement * distance_falloff * magnetic_field_strength * beaming_factor;
    base_intensity *= (1.0 + 0.5 * abs(kh_instability) + 0.3 * abs(rt_instability));
    base_intensity *= knot_enhancement;
    
    vec4 thermal_emission = BLACK_BODY_COLOR(plasma_temp);
    vec3 final_color = thermal_emission.rgb * 0.4 + synchrotron_color * 0.6;
    final_color += helix_color * helix_vis_strength * 0.8;
    if (doppler_factor > 1.2) {
        final_color = final_color * vec3(0.8, 0.9, 1.3);
    } else if (doppler_factor < 0.8) {
        final_color = final_color * vec3(1.3, 0.9, 0.7);
    }
    
    float opacity = exp(-0.02 * z_jet) * cooling_factor;
    float edge_fade = smoothstep(sheath_radius, sheath_radius * 0.8, rho_jet);
    opacity *= edge_fade;
    if (z_jet > 40.0) {
        float shock_dist = z_jet - 40.0;
        float shock_strength = exp(-shock_dist * 0.1) * sin(shock_dist * 0.5 + time * 2.0) * 0.5 + 0.5;
        final_color += vec3(1.0, 0.8, 0.6) * shock_strength * 0.4;
        base_intensity += shock_strength;
    }
    
    float activity_cycle = 0.9 + 0.1 * sin(time * 0.3);
    float precession = 1.0 + 0.05 * sin(time * 0.8 + phi * 2.0);
    base_intensity *= activity_cycle * precession;
    
    return final_color * base_intensity * opacity * 0.8;
}
void main() {

    {{#planetEnabled}}
    PLANET_RADIUS = planet_radius;
    PLANET_DISTANCE = max(planet_distance,planet_radius + 1.5 * bh_radius);
    PLANET_ORBITAL_ANG_VEL = -1.0 / sqrt(2.0*(PLANET_DISTANCE-bh_radius)) / PLANET_DISTANCE;
    float MAX_PLANET_ROT = max((1.0 + PLANET_ORBITAL_ANG_VEL*PLANET_DISTANCE) / PLANET_RADIUS,0.0);
    PLANET_ROTATION_ANG_VEL = -PLANET_ORBITAL_ANG_VEL + MAX_PLANET_ROT * 0.5;
    PLANET_GAMMA = 1.0/sqrt(1.0-SQ(PLANET_ORBITAL_ANG_VEL*PLANET_DISTANCE));
    {{/planetEnabled}}

    vec2 p = -1.0 + 2.0 * gl_FragCoord.xy / resolution.xy;
    p.y *= resolution.y / resolution.x;

    vec3 pos = cam_pos;
    vec3 ray = normalize(p.x*cam_x + p.y*cam_y + FOV_MULT*cam_z);

    {{#aberration}}
    ray = lorentz_velocity_transformation(ray, cam_vel);
    {{/aberration}}

    float ray_intensity = 1.0;
    float ray_doppler_factor = 1.0;

    float gamma = 1.0/sqrt(1.0-dot(cam_vel,cam_vel));
    ray_doppler_factor = gamma*(1.0 + dot(ray,-cam_vel));
    {{#beaming}}
    ray_intensity /= ray_doppler_factor*ray_doppler_factor*ray_doppler_factor;
    {{/beaming}}
    {{^doppler_shift}}
    ray_doppler_factor = 1.0;
    {{/doppler_shift}}

    float step = 0.01;
    vec4 color = vec4(0.0,0.0,0.0,1.0);

    float u = 1.0 / length(pos), old_u;
    float u0 = u;

    vec3 normal_vec = normalize(pos);
    vec3 tangent_vec = normalize(cross(cross(normal_vec, ray), normal_vec));

    float du = -dot(ray,normal_vec) / dot(ray,tangent_vec) * u;
    float du0 = du;

    float phi = 0.0;
    float t = time;
    float dt = 1.0;

    {{^light_travel_time}}
    float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
    vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
    {{/light_travel_time}}

    vec3 old_pos;

    // The event horizon is at r = bh_radius, so u = 1/bh_radius
    float HORIZON_U = 1.0 / bh_radius;

    for (int j=0; j < NSTEPS; j++) {

        step = MAX_REVOLUTIONS * 2.0*M_PI / float(NSTEPS);

        float max_rel_u_change = (1.0-log(u))*10.0 / float(NSTEPS);
        if ((du > 0.0 || (du0 < 0.0 && u0/u < 5.0)) && abs(du) > abs(max_rel_u_change*u) / step)
            step = max_rel_u_change*u/abs(du);

        old_u = u;

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        dt = sqrt(du*du + u*u*(1.0-u*bh_radius))/(u*u*(1.0-u*bh_radius))*step;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        u += du*step;
        float ddu = -u*(1.0 - 1.5*u*u*bh_radius*bh_radius);
        du += ddu*step;

        if (u < 0.0) break;

        phi += step;

        old_pos = pos;
        pos = (cos(phi)*normal_vec + sin(phi)*tangent_vec)/u;

        ray = pos-old_pos;
        float solid_isec_t = 2.0;
        float ray_l = length(ray);

        {{#light_travel_time}}
        {{#gravitational_time_dilation}}
        float mix = smooth_step_custom(1.0/u, 8.0 * bh_radius, 1.0 / bh_radius);
        dt = mix*ray_l + (1.0-mix)*dt;
        {{/gravitational_time_dilation}}
        {{^gravitational_time_dilation}}
        dt = ray_l;
        {{/gravitational_time_dilation}}
        {{/light_travel_time}}

        {{#planetEnabled}}
        if (
            (
                old_pos.z * pos.z < 0.0 ||
                min(abs(old_pos.z), abs(pos.z)) < PLANET_RADIUS
            ) &&
            max(u, old_u) > 1.0/(PLANET_DISTANCE+PLANET_RADIUS) &&
            min(u, old_u) < 1.0/(PLANET_DISTANCE-PLANET_RADIUS)
        ) {

            {{#light_travel_time}}
            float planet_ang0 = t * PLANET_ORBITAL_ANG_VEL;
            vec3 planet_pos0 = vec3(cos(planet_ang0), sin(planet_ang0), 0)*PLANET_DISTANCE;
            {{/light_travel_time}}

            vec4 planet_isec = planet_intersection(old_pos, ray, t, dt,
                    planet_pos0, ray_doppler_factor);
            if (planet_isec.w > 0.0) {
                solid_isec_t = planet_isec.w;
                planet_isec.w = 1.0;
                color += planet_isec;
            }
        }
        {{/planetEnabled}}

        {{#magneticFieldEnabled}}
        // Magnetic field visualization
        float field_line_intensity = fieldLineIntensity(pos, ray);
        if (field_line_intensity > 0.01) {
            // Use new magnetic field color system instead of uniform field_color
            vec3 field_contribution = getMagneticFieldColor(pos, field_line_intensity);
            vec3 plasma_contribution = plasmaGlow(pos, field_line_intensity);

            // Combine field and plasma contributions
            vec3 combined_contribution = field_contribution + plasma_contribution;
            float blend_factor = min(field_line_intensity * 0.3, 1.0);

            // Blend with existing color using proper vector addition
            color.rgb += combined_contribution * blend_factor;
        }
        
        // Enhanced polar jets emission
        vec3 jet_emission = polarJetEmission(pos, ray);
        if (length(jet_emission) > 0.01) {
            color.rgb += jet_emission;
        }
        {{/magneticFieldEnabled}}

        {{#accretion_disk}}
        if (old_pos.z * pos.z < 0.0) {
            float acc_isec_t = -old_pos.z / ray.z;
            if (acc_isec_t < solid_isec_t && acc_isec_t > 0.0 && acc_isec_t < 1.0) {
                vec3 isec = old_pos + ray*acc_isec_t;
                float r = length(isec);

                float accretion_min_r = ACCRETION_MIN_R_FACTOR * bh_radius;
                float accretion_width = ACCRETION_WIDTH_FACTOR * bh_radius;

                if (r > accretion_min_r && r < accretion_min_r + accretion_width) {
                    // Enhanced 3D turbulence with multiple octaves
                    vec3 noise_coord_3d = isec * noise_scale;

                    // Orbital velocity and shear
                    float orbital_velocity = 1.0 / sqrt(r * r * r);
                    float shear_rate = 1.5 * orbital_velocity / r;

                    // Time-dependent rotation with proper Keplerian motion
                    float angle = time * orbital_velocity * noise_speed * 2.0;
                    mat2 rot = mat2(cos(angle), -sin(angle), sin(angle), cos(angle));

                    // Apply rotation to 2D components
                    noise_coord_3d.xy = rot * noise_coord_3d.xy;

                    // Multi-scale turbulence for realistic structure
                    float density = 0.0;
                    float amplitude = 1.0;
                    float frequency = 1.0;

                    // Large-scale spiral structure
                    vec2 spiral_coord = noise_coord_3d.xy;
                    float spiral_angle = atan(spiral_coord.y, spiral_coord.x);
                    float spiral_r = length(spiral_coord);
                    spiral_angle += log(spiral_r) * 0.3; // Logarithmic spiral
                    spiral_coord = vec2(cos(spiral_angle), sin(spiral_angle)) * spiral_r;

                    for (int i = 0; i < 4; i++) {
                        // 3D noise for vertical structure
                        vec3 sample_coord = noise_coord_3d * frequency;
                        sample_coord.z += time * noise_speed * 0.5;

                        // Simulate turbulent eddies
                        float noise_val = fbm(sample_coord.xy + sample_coord.z * 0.1);

                        // Add vertical turbulence
                        float vertical_noise = fbm(vec2(spiral_r * 2.0, sample_coord.z * 3.0));
                        noise_val += vertical_noise * 0.3;

                        // Shear instabilities
                        float shear_coord = spiral_coord.x + shear_rate * time * noise_speed;
                        float shear_noise = fbm(vec2(shear_coord, spiral_coord.y) * frequency * 0.7);
                        noise_val += shear_noise * 0.4;

                        density += noise_val * amplitude;
                        amplitude *= 0.6;
                        frequency *= 2.1;
                    }

                    // Normalize and enhance contrast
                    density = density * 0.5 + 0.5;
                    density = pow(density, 1.5);

                    // Radial density profile with more realistic falloff
                    float radial_profile = exp(-2.0 * (r - accretion_min_r) / accretion_width);
                    radial_profile *= (1.0 - smoothstep(accretion_min_r + accretion_width * 0.8, accretion_min_r + accretion_width, r));

                    // Add inner edge enhancement (hot inner rim)
                    float inner_rim = exp(-10.0 * abs(r - accretion_min_r * 1.2));
                    radial_profile += inner_rim * 0.5;

                    density *= radial_profile;

                    // Improved visibility threshold
                    if (density > 0.02) {
                        // More realistic intensity scaling
                        float accretion_intensity = ACCRETION_BRIGHTNESS * density * 3.0;

                        // Temperature gradient with turbulent heating
                        float base_temperature = ACCRETION_TEMPERATURE / pow(r/bh_radius, 0.75);

                        // Turbulent heating increases temperature in high-density regions
                        float turbulent_heating = 1.0 + density * density * 2.0;

                        // Magnetic heating near the inner edge
                        float magnetic_heating = 1.0 + 0.5 * exp(-3.0 * (r - accretion_min_r));

                        float temperature = base_temperature * turbulent_heating * magnetic_heating;

                        // Realistic velocity field with turbulent components
                        vec3 accretion_v = vec3(-isec.y, isec.x, 0.0) * orbital_velocity;

                        // Add turbulent velocity perturbations
                        vec3 turb_v = vec3(
                            (fbm(noise_coord_3d.yz + time) - 0.5) * 0.1,
                            (fbm(noise_coord_3d.xz + time + 100.0) - 0.5) * 0.1,
                            (fbm(noise_coord_3d.xy + time + 200.0) - 0.5) * 0.05
                        ) * orbital_velocity;

                        accretion_v += turb_v;

                        gamma = 1.0/sqrt(1.0-dot(accretion_v,accretion_v));
                        float doppler_factor = gamma*(1.0+dot(normalize(ray),accretion_v));

                        {{#beaming}}
                        // Reduce the intensity of Doppler beaming to prevent excessive brightness
                        float reduced_doppler = 1.0 + (doppler_factor - 1.0) * 0.6;
                        accretion_intensity /= reduced_doppler*reduced_doppler;
                        {{/beaming}}
                        {{#doppler_shift}}
                        temperature /= ray_doppler_factor*doppler_factor;
                        {{/doppler_shift}}

                        // Enhanced color mixing for more realistic appearance
                        vec4 disk_color = BLACK_BODY_COLOR(temperature) * accretion_intensity;

                        // Add some scattering effects for realism
                        float scattering = 1.0 + 0.2 * density;
                        disk_color.rgb *= scattering;

                        color += disk_color;
                    }
                }
            }
        }
        {{/accretion_disk}}

        {{#light_travel_time}}
        t -= dt;
        {{/light_travel_time}}

        if (solid_isec_t <= 1.0) u = HORIZON_U + 1.0; // break
        if (u > HORIZON_U) break;
    }

    if (u > HORIZON_U) {
        // Ray hit the event horizon
        color += get_hawking_radiation_color(pos, u);
    } else {
        // Ray escaped to infinity
        {{^hawkingRadiationEnabled}}
            ray = normalize(pos - old_pos);
            vec2 tex_coord = sphere_map(ray * BG_COORDS);
            float t_coord;

            vec4 star_color = texture2D(star_texture, tex_coord);
            if (star_color.r > 0.0) {
                t_coord = (STAR_MIN_TEMPERATURE +
                    (STAR_MAX_TEMPERATURE-STAR_MIN_TEMPERATURE) * star_color.g)
                     / ray_doppler_factor;

                color += BLACK_BODY_COLOR(t_coord) * star_color.r * STAR_BRIGHTNESS;
            }

            color += galaxy_color(tex_coord, ray_doppler_factor) * GALAXY_BRIGHTNESS;
        {{/hawkingRadiationEnabled}}
    }

    gl_FragColor = color*ray_intensity;
}
