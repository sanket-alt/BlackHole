/**
 * MagneticFieldTopology.js - RK4 integration of magnetic field lines
 * Traces field lines using: dx/dt = B(x) / |B(x)|
 */

class MagneticFieldTopology {
  constructor(kerrMetric) {
    this.metric = kerrMetric;
    this.maxSteps = 5000;
    this.tolerance = 1e-6;
  }

  /**
   * Dipole magnetic field in Boyer-Lindquist coordinates
   */
  magneticFieldVector(r, theta, phi, strength = 1.0, tilt = 0.0) {
    const r3 = r * r * r;
    const sin_t = Math.sin(theta);
    const cos_t = Math.cos(theta);
    
    let B_r = (2.0 * strength * cos_t) / r3;
    let B_theta = (strength * sin_t) / r3;
    const B_phi = 0.0;
    
    if (Math.abs(tilt) > 1e-6) {
      const cos_tilt = Math.cos(tilt);
      const sin_tilt = Math.sin(tilt);
      
      const B_r_old = B_r;
      const B_z_old = B_theta;
      
      B_r = B_r_old * cos_tilt + B_z_old * sin_tilt;
      B_theta = -B_r_old * sin_tilt + B_z_old * cos_tilt;
    }
    
    return { B_r, B_theta, B_phi, magnitude: Math.sqrt(B_r*B_r + B_theta*B_theta) };
  }

  /**
   * Mixed magnetic field: Poloidal + Toroidal components
   */
  mixedMagneticField(r, theta, phi, strength = 1.0, mixParam = 0.5, tilt = 0.0) {
    const poloidal = this.magneticFieldVector(r, theta, phi, strength * (1 - mixParam), tilt);
    const B_phi_toroidal = (strength * mixParam) / (r * r);
    
    return {
      B_r: poloidal.B_r,
      B_theta: poloidal.B_theta,
      B_phi: B_phi_toroidal,
      magnitude: Math.sqrt(
        poloidal.B_r * poloidal.B_r +
        poloidal.B_theta * poloidal.B_theta +
        B_phi_toroidal * B_phi_toroidal
      ),
    };
  }

  /**
   * RK4 step for field line integration
   */
  rk4Step(r, theta, phi, stepSize, strength, mixParam, tilt) {
    const B1 = this.mixedMagneticField(r, theta, phi, strength, mixParam, tilt);
    const mag1 = Math.max(B1.magnitude, 1e-10);
    const k1_r = B1.B_r / mag1;
    const k1_theta = B1.B_theta / mag1;
    const k1_phi = B1.B_phi / mag1;
    
    const r2 = r + k1_r * stepSize * 0.5;
    const theta2 = theta + k1_theta * stepSize * 0.5;
    const phi2 = phi + k1_phi * stepSize * 0.5;
    const B2 = this.mixedMagneticField(r2, theta2, phi2, strength, mixParam, tilt);
    const mag2 = Math.max(B2.magnitude, 1e-10);
    const k2_r = B2.B_r / mag2;
    const k2_theta = B2.B_theta / mag2;
    const k2_phi = B2.B_phi / mag2;
    
    const r3 = r + k2_r * stepSize * 0.5;
    const theta3 = theta + k2_theta * stepSize * 0.5;
    const phi3 = phi + k2_phi * stepSize * 0.5;
    const B3 = this.mixedMagneticField(r3, theta3, phi3, strength, mixParam, tilt);
    const mag3 = Math.max(B3.magnitude, 1e-10);
    const k3_r = B3.B_r / mag3;
    const k3_theta = B3.B_theta / mag3;
    const k3_phi = B3.B_phi / mag3;
    
    const r4 = r + k3_r * stepSize;
    const theta4 = theta + k3_theta * stepSize;
    const phi4 = phi + k3_phi * stepSize;
    const B4 = this.mixedMagneticField(r4, theta4, phi4, strength, mixParam, tilt);
    const mag4 = Math.max(B4.magnitude, 1e-10);
    const k4_r = B4.B_r / mag4;
    const k4_theta = B4.B_theta / mag4;
    const k4_phi = B4.B_phi / mag4;
    
    return {
      r: r + (k1_r + 2*k2_r + 2*k3_r + k4_r) * stepSize / 6.0,
      theta: theta + (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta) * stepSize / 6.0,
      phi: phi + (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi) * stepSize / 6.0,
    };
  }

  /**
   * Trace a single magnetic field line
   */
  traceFieldLine(r0, theta0, strength = 1.0, mixParam = 0.5, tilt = 0.0, direction = 1.0) {
    const fieldLine = [{ r: r0, theta: theta0, x: 0, y: 0, z: 0 }];
    
    let r = r0;
    let theta = theta0;
    let phi = 0;
    
    const stepSize = direction * 0.05;
    
    for (let step = 0; step < this.maxSteps; step++) {
      if (r < this.metric.r_plus * 1.05 || r > 100) break;
      if (theta < 0.01 || theta > Math.PI - 0.01) break;
      
      const newPos = this.rk4Step(r, theta, phi, stepSize, strength, mixParam, tilt);
      
      r = newPos.r;
      theta = newPos.theta;
      phi = newPos.phi;
      
      const cart = this.metric.toCartesian(r, theta, phi);
      
      fieldLine.push({
        r, theta,
        x: cart.x, y: cart.y, z: cart.z,
      });
    }
    
    return fieldLine;
  }

  /**
   * Trace multiple field lines in a grid pattern
   */
  traceFieldLineGrid(numLines = 32, strength = 1.0, mixParam = 0.5, tilt = 0.0) {
    const fieldLines = [];
    const r_start = this.metric.r_isco;
    
    for (let i = 0; i < numLines; i++) {
      const theta0 = Math.PI / 2;
      const line_forward = this.traceFieldLine(r_start, theta0, strength, mixParam, tilt, 1.0);
      const line_backward = this.traceFieldLine(r_start, theta0, strength, mixParam, tilt, -1.0);
      
      fieldLines.push({
        forward: line_forward,
        backward: line_backward,
        index: i,
      });
    }
    
    return fieldLines;
  }

  /**
   * Extract field line data for WebGL visualization
   */
  getFieldLineVertices(fieldLines) {
    const vertices = [];
    const colors = [];
    
    for (const line of fieldLines) {
      for (const point of line.forward) {
        vertices.push(point.x, point.y, point.z);
        colors.push(0.0, 0.8, 1.0, 1.0);
      }
      
      for (const point of line.backward) {
        vertices.push(point.x, point.y, point.z);
        colors.push(1.0, 0.2, 0.8, 1.0);
      }
    }
    
    return { vertices: new Float32Array(vertices), colors: new Float32Array(colors) };
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = MagneticFieldTopology;
}
