
// Add these uniforms to your raytracer.glsl fragment shader:
uniform float magnetic_strength;
uniform float magnetic_dipole_tilt;
uniform float field_lines_density;
uniform vec3 field_color;
uniform int show_plasma_effects;

// Magnetic dipole field calculation
vec3 magneticField(vec3 pos) {
    float r = length(pos);
    if (r < 1.1) return vec3(0.0); // Inside event horizon
    
    // Convert to spherical coordinates
    float theta = acos(pos.z / r);
    float phi = atan(pos.y, pos.x);
    
    // Apply dipole tilt
    float tilt_cos = cos(magnetic_dipole_tilt);
    float tilt_sin = sin(magnetic_dipole_tilt);
    
    // Magnetic dipole field components
    float r3 = r * r * r;
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    
    // Dipole field strength (simplified)
    float B_r = (2.0 * magnetic_strength * cos_theta) / r3;
    float B_theta = (magnetic_strength * sin_theta) / r3;
    
    // Convert back to Cartesian coordinates
    vec3 B_field;
    B_field.x = B_r * sin_theta * cos(phi) + B_theta * cos_theta * cos(phi);
    B_field.y = B_r * sin_theta * sin(phi) + B_theta * cos_theta * sin(phi);
    B_field.z = B_r * cos_theta - B_theta * sin_theta;
    
    return B_field * tilt_cos + cross(B_field, vec3(0, 0, 1)) * tilt_sin;
}

// Field line visualization
float fieldLineIntensity(vec3 pos, vec3 ray_dir) {
    vec3 B = magneticField(pos);
    float B_mag = length(B);
    
    if (B_mag < 0.001) return 0.0;
    
    vec3 B_normalized = B / B_mag;
    
    // Create field line pattern using sin/cos functions
    float line_param = dot(pos, B_normalized) * field_lines_density;
    float line_intensity = sin(line_param * 6.28318) * cos(line_param * 3.14159);
    line_intensity = pow(abs(line_intensity), 4.0);
    
    // Fade with distance and make lines thinner closer to black hole
    float r = length(pos);
    float fade = exp(-r * 0.1) * (r - 1.0) / r;
    
    return line_intensity * fade * magnetic_strength;
}

// Plasma effects around magnetic field lines
vec3 plasmaGlow(vec3 pos, float field_strength) {
    if (show_plasma_effects == 0) return vec3(0.0);
    
    float r = length(pos);
    float plasma_intensity = field_strength * exp(-r * 0.5);
    
    // Color temperature effect - hotter near field lines
    vec3 hot_color = vec3(1.0, 0.8, 0.6);  // Hot white-yellow
    vec3 cool_color = vec3(0.6, 0.8, 1.0); // Cool blue
    
    return mix(cool_color, hot_color, plasma_intensity) * plasma_intensity;
}

// Add this to your main raytracing loop in raytracer.glsl:
{{#magneticFieldEnabled}}
    // Magnetic field visualization
    float field_line_intensity = fieldLineIntensity(ray_pos, ray_dir);
    if (field_line_intensity > 0.01) {
        vec3 field_contribution = field_color * field_line_intensity;
        vec3 plasma_contribution = plasmaGlow(ray_pos, field_line_intensity);
        
        // Blend with existing color
        color = mix(color, field_contribution + plasma_contribution, field_line_intensity * 0.3);
    }
{{/magneticFieldEnabled}}
