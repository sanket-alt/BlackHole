/**
 * GeodesicIntegrator.js - RK4 integrator for null geodesics in Kerr spacetime
 * Solves the geodesic equation: du^μ/dλ + Γ^μ_νρ u^ν u^ρ = 0
 */

class GeodesicIntegrator {
  constructor(kerrMetric) {
    this.metric = kerrMetric;
    this.maxSteps = 10000;
    this.tolerance = 1e-8;
  }

  /**
   * RK4 step for null geodesic integration
   */
  rk4Step(state, stepSize) {
    const k1 = this.geodesicDerivative(state);
    const midState1 = this.addScaledState(state, k1, stepSize * 0.5);
    const k2 = this.geodesicDerivative(midState1);
    const midState2 = this.addScaledState(state, k2, stepSize * 0.5);
    const k3 = this.geodesicDerivative(midState2);
    const endState = this.addScaledState(state, k3, stepSize);
    const k4 = this.geodesicDerivative(endState);
    
    const combined = this.addStates(
      k1,
      this.scaleState(this.addStates(k2, k3), 2.0),
      k4
    );
    
    return this.addScaledState(state, combined, stepSize / 6.0);
  }

  /**
   * Geodesic equation: du^μ/dλ + Γ^μ_νρ u^ν u^ρ = 0
   */
  geodesicDerivative(state) {
    const r = state.r;
    const theta = state.theta;
    const u_r = state.u_r;
    const u_theta = state.u_theta;
    const u_phi = state.u_phi;
    const u_t = state.u_t;
    
    const u = [u_t, u_r, u_theta, u_phi];
    
    const drdt = u_r;
    const dtheta_dt = u_theta;
    const dphi_dt = u_phi;
    
    const dudt = [0, 0, 0, 0];
    
    for (let mu = 0; mu < 4; mu++) {
      let accel = 0;
      for (let nu = 0; nu < 4; nu++) {
        for (let rho = 0; rho < 4; rho++) {
          const gamma = this.metric.christoffelSymbol(r, theta, mu, nu, rho);
          accel -= gamma * u[nu] * u[rho];
        }
      }
      dudt[mu] = accel;
    }
    
    return {
      r: drdt,
      theta: dtheta_dt,
      phi: dphi_dt,
      u_t: dudt[0],
      u_r: dudt[1],
      u_theta: dudt[2],
      u_phi: dudt[3],
    };
  }

  /**
   * Integrate a null geodesic from observer to black hole
   */
  integrateRay(initialState, maxAffineParam = 100.0) {
    const trajectory = [{ ...initialState }];
    let state = { ...initialState };
    let affineParam = 0;
    let stepSize = 0.01;
    
    for (let step = 0; step < this.maxSteps; step++) {
      const horizonDist = state.r - this.metric.r_plus;
      if (horizonDist < this.metric.r_plus * 0.2) {
        stepSize = 0.001;
      } else if (horizonDist < this.metric.r_plus * 0.5) {
        stepSize = 0.005;
      } else {
        stepSize = 0.01;
      }
      
      const newState = this.rk4Step(state, stepSize);
      
      if (newState.r < 0 || newState.r > 1000) break;
      if (isNaN(newState.r) || isNaN(newState.theta)) break;
      
      affineParam += stepSize;
      if (affineParam > maxAffineParam) break;
      
      trajectory.push({ ...newState });
      state = newState;
    }
    
    return {
      trajectory,
      affineParameter: affineParam,
      finalState: state,
      stepsTaken: trajectory.length,
    };
  }

  /**
   * Compute initial four-velocity for a ray
   */
  computeInitialFourVelocity(r, theta, rayDir) {
    const u_r = rayDir.r || 0.1;
    const u_theta = rayDir.theta || 0;
    const u_phi = rayDir.phi || 0;
    
    return this.metric.normalizeNullGeodesic(r, theta, u_r, u_theta, u_phi);
  }

  addStates(...states) {
    const sum = {
      r: 0, theta: 0, phi: 0,
      u_t: 0, u_r: 0, u_theta: 0, u_phi: 0,
    };
    
    for (const state of states) {
      sum.r += state.r;
      sum.theta += state.theta;
      sum.phi += state.phi;
      sum.u_t += state.u_t;
      sum.u_r += state.u_r;
      sum.u_theta += state.u_theta;
      sum.u_phi += state.u_phi;
    }
    return sum;
  }

  scaleState(state, scalar) {
    return {
      r: state.r * scalar,
      theta: state.theta * scalar,
      phi: state.phi * scalar,
      u_t: state.u_t * scalar,
      u_r: state.u_r * scalar,
      u_theta: state.u_theta * scalar,
      u_phi: state.u_phi * scalar,
    };
  }

  addScaledState(baseState, deltaState, scale) {
    return {
      r: baseState.r + deltaState.r * scale,
      theta: baseState.theta + deltaState.theta * scale,
      phi: baseState.phi + deltaState.phi * scale,
      u_t: baseState.u_t + deltaState.u_t * scale,
      u_r: baseState.u_r + deltaState.u_r * scale,
      u_theta: baseState.u_theta + deltaState.u_theta * scale,
      u_phi: baseState.u_phi + deltaState.u_phi * scale,
    };
  }

  checkNullConstraint(state) {
    const g = this.metric.metricTensor(state.r, state.theta);
    const u = [state.u_t, state.u_r, state.u_theta, state.u_phi];
    
    let sum = 0;
    for (let mu = 0; mu < 4; mu++) {
      for (let nu = 0; nu < 4; nu++) {
        sum += g[mu][nu] * u[mu] * u[nu];
      }
    }
    return Math.abs(sum);
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = GeodesicIntegrator;
}
