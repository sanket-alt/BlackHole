/**
 * KerrMetric.js - Comprehensive Kerr geometry calculations
 * Handles Boyer-Lindquist coordinates, metric tensor, Christoffel symbols,
 * and derived quantities (event horizon, ergosphere, ISCO)
 */

class KerrMetric {
  constructor(mass = 1.0, spin = 0.0) {
    this.M = mass;
    // Enforce Kerr bound: |a| < M
    this.a = Math.min(Math.abs(spin), 0.998 * mass);
    
    // Pre-compute horizon and other critical radii
    this.updateCriticalRadii();
  }

  updateCriticalRadii() {
    const M = this.M;
    const a = this.a;
    
    // Event horizon radius: r₊ = M + √(M² - a²)
    const discriminant = M * M - a * a;
    this.r_plus = M + Math.sqrt(Math.max(0, discriminant));
    
    // Cauchy horizon: r₋ = M - √(M² - a²)
    this.r_minus = M - Math.sqrt(Math.max(0, discriminant));
    
    // ISCO (prograde): r_ISCO = M[3 + Z₂ - √((3 - Z₁)(3 + Z₁ + 2Z₂))]
    this.computeISCO();
    
    // Frame-dragging radius (where ω = Ω_F)
    this.computeFrameDragging();
  }

  /**
   * ISCO radius for prograde orbits
   * r_ISCO = M(3 + Z₂ - √((3 - Z₁)(3 + Z₁ + 2Z₂)))
   * where Z₁, Z₂ involve the spin parameter q = a/M
   */
  computeISCO() {
    const M = this.M;
    const q = this.a / M;  // Dimensionless spin
    
    const term1 = (1 - q * q) ** (1/3);
    const term2 = (1 + q) ** (1/3) + (1 - q) ** (1/3);
    
    const Z1 = 1 + term1 * term2;
    const Z2 = Math.sqrt(Z1 * Z1 - 4 * q * q);
    
    this.r_isco = M * (3 + Z2 - Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)));
  }

  /**
   * Frame-dragging effect: angular velocity of inertial frame
   * Ω_F = 2ar₊ / [r₊² + a² + 2a²r₊/r₊² + a²]
   */
  computeFrameDragging() {
    const r_plus = this.r_plus;
    const a = this.a;
    
    const numerator = 2 * a * r_plus;
    const denominator = (r_plus * r_plus + a * a) * (r_plus * r_plus + a * a + 2 * a * a * r_plus / (r_plus * r_plus + a * a));
    
    this.omega_F = numerator / Math.max(denominator, 1e-10);
  }

  /**
   * Boyer-Lindquist to Cartesian coordinates
   */
  toCartesian(r, theta, phi) {
    const a = this.a;
    const sigma = r * r + a * a * Math.cos(theta) * Math.cos(theta);
    const sqrtSigma = Math.sqrt(sigma);
    
    return {
      x: sqrtSigma * Math.sin(theta) * Math.cos(phi),
      y: sqrtSigma * Math.sin(theta) * Math.sin(phi),
      z: r * Math.cos(theta),
    };
  }

  /**
   * Cartesian to Boyer-Lindquist coordinates
   */
  toBoyer Lindquist(x, y, z) {
    const a = this.a;
    const rho_sq = x * x + y * y;
    const r_sq = rho_sq + z * z;
    
    // Solve for r using quadratic formula
    // r² - 2ar_perp*sin(phi) + (rho² + (z² - a²))
    // Approximate: r ≈ √(rho² + z²) for small a
    
    const r = Math.sqrt(rho_sq + z * z);
    const theta = Math.atan2(Math.sqrt(rho_sq), z);
    const phi = Math.atan2(y, x);
    
    return { r, theta, phi };
  }

  /**
   * Metric tensor components g_μν in Boyer-Lindquist coordinates
   * Returns a 4x4 matrix
   */
  metricTensor(r, theta) {
    const a = this.a;
    const M = this.M;
    
    const sigma = r * r + a * a * Math.cos(theta) * Math.cos(theta);
    const delta = r * r - 2 * M * r + a * a;
    const sin_t = Math.sin(theta);
    const sin2_t = sin_t * sin_t;
    const cos_t = Math.cos(theta);
    const cos2_t = cos_t * cos_t;
    
    // Metric components
    const g_00 = -(1 - 2 * M * r / sigma);
    const g_03 = -4 * M * a * r * sin2_t / sigma;
    const g_11 = sigma / delta;
    const g_22 = sigma;
    const g_33 = (r * r + a * a + 2 * M * a * a * r * sin2_t / sigma) * sin2_t;
    
    return [
      [g_00, 0, 0, g_03],
      [0, g_11, 0, 0],
      [0, 0, g_22, 0],
      [g_03, 0, 0, g_33],
    ];
  }

  /**
   * Inverse metric tensor g^μν
   */
  inverseMetricTensor(r, theta) {
    const g = this.metricTensor(r, theta);
    
    // For this specific block-diagonal structure, compute inverse directly
    const a = this.a;
    const M = this.M;
    const sigma = r * r + a * a * Math.cos(theta) * Math.cos(theta);
    const delta = r * r - 2 * M * r + a * a;
    const sin_t = Math.sin(theta);
    const sin2_t = sin_t * sin_t;
    
    const g_00 = -(1 - 2 * M * r / sigma);
    const g_03 = -4 * M * a * r * sin2_t / sigma;
    const g_11 = sigma / delta;
    const g_22 = sigma;
    const g_33 = (r * r + a * a + 2 * M * a * a * r * sin2_t / sigma) * sin2_t;
    
    // Determinant of 2x2 block [g_00, g_03; g_03, g_33]
    const det_block = g_00 * g_33 - g_03 * g_03;
    
    const g_inv_00 = g_33 / det_block;
    const g_inv_03 = -g_03 / det_block;
    const g_inv_11 = 1 / g_11;
    const g_inv_22 = 1 / g_22;
    const g_inv_33 = g_00 / det_block;
    
    return [
      [g_inv_00, 0, 0, g_inv_03],
      [0, g_inv_11, 0, 0],
      [0, 0, g_inv_22, 0],
      [g_inv_03, 0, 0, g_inv_33],
    ];
  }

  /**
   * Christoffel symbol Γ^μ_νρ
   * Returns only non-zero components for efficiency
   * This is critical for geodesic integration
   */
  christoffelSymbol(r, theta, mu, nu, rho) {
    const a = this.a;
    const M = this.M;
    
    const sigma = r * r + a * a * Math.cos(theta) * Math.cos(theta);
    const delta = r * r - 2 * M * r + a * a;
    const sin_t = Math.sin(theta);
    const sin2_t = sin_t * sin_t;
    const cos_t = Math.cos(theta);
    const cos2_t = cos_t * cos_t;
    
    // Only compute non-zero Christoffel symbols
    // Due to symmetry of the metric, many are zero
    
    // Γ^r_tt
    if (mu === 0 && nu === 0 && rho === 0) {
      return (M * (sigma - 2 * a * a * cos2_t)) / (sigma * sigma * delta);
    }
    
    // Γ^r_θθ
    if (mu === 0 && nu === 2 && rho === 2) {
      return -r / delta;
    }
    
    // Γ^r_φφ
    if (mu === 0 && nu === 3 && rho === 3) {
      const A = r * r + a * a;
      const num = -sin2_t * r * A * (A + 2 * M * a * a * sin2_t / sigma) + a * a * M * delta * sin2_t;
      return num / (sigma * sigma * delta);
    }
    
    // Γ^r_tφ
    if ((mu === 0 && nu === 0 && rho === 3) || (mu === 0 && nu === 3 && rho === 0)) {
      return (2 * M * a * r) / (sigma * sigma * delta);
    }
    
    // Γ^θ_tt
    if (mu === 2 && nu === 0 && rho === 0) {
      return (-2 * M * a * a * sin_t * cos_t) / (sigma * sigma);
    }
    
    // Γ^θ_rr
    if (mu === 2 && nu === 1 && rho === 1) {
      return -(a * a * sin_t * cos_t) / sigma;
    }
    
    // Γ^θ_φφ
    if (mu === 2 && nu === 3 && rho === 3) {
      return -sin_t * cos_t * (r * r + a * a + 2 * M * a * a * r * sin2_t / sigma);
    }
    
    // Γ^θ_θθ
    if (mu === 2 && nu === 2 && rho === 2) {
      return (a * a * sin_t * cos_t) / sigma;
    }
    
    // Γ^θ_tφ
    if ((mu === 2 && nu === 0 && rho === 3) || (mu === 2 && nu === 3 && rho === 0)) {
      return (4 * M * a * r * sin_t * cos_t) / (sigma * sigma);
    }
    
    // Γ^φ_rφ
    if ((mu === 3 && nu === 1 && rho === 3) || (mu === 3 && nu === 3 && rho === 1)) {
      return 1 / r;
    }
    
    // Γ^φ_θθ
    if ((mu === 3 && nu === 2 && rho === 2) || (mu === 3 && nu === 2 && rho === 2)) {
      return cos_t / sin_t;
    }
    
    // Γ^t_rφ
    if ((mu === 1 && nu === 0 && rho === 3) || (mu === 1 && nu === 3 && rho === 0)) {
      const num = 2 * M * a * r * (a * a + r * r);
      const denom = sigma * sigma * delta;
      return num / denom;
    }
    
    return 0.0;
  }

  /**
   * Ergosphere radius as function of polar angle θ
   * r_ergo(θ) = M + √(M² - a² cos²θ)
   */
  ergosphereRadius(theta) {
    const cos_t = Math.cos(theta);
    const discriminant = this.M * this.M - this.a * this.a * cos_t * cos_t;
    return this.M + Math.sqrt(Math.max(0, discriminant));
  }

  /**
   * Normalize a null geodesic: ensure g_μν u^μ u^ν = 0
   */
  normalizeNullGeodesic(r, theta, u_r, u_theta, u_phi) {
    const sigma = r * r + this.a * this.a * Math.cos(theta) * Math.cos(theta);
    const delta = r * r - 2 * this.M * r + this.a * this.a;
    const sin2_t = Math.sin(theta) * Math.sin(theta);
    
    // Constraint: g_μν u^μ u^ν = 0
    // Solve for u_t (time component)
    
    const g_00 = -(1 - 2 * this.M * r / sigma);
    const g_03 = -4 * this.M * this.a * r * sin2_t / sigma;
    const g_11 = sigma / delta;
    const g_22 = sigma;
    const g_33 = (r * r + this.a * this.a + 2 * this.M * this.a * this.a * r * sin2_t / sigma) * sin2_t;
    
    const spatial_part = g_11 * u_r * u_r + g_22 * u_theta * u_theta + g_33 * u_phi * u_phi;
    
    // Quadratic formula: g_00 (u_t)² + 2 g_03 u_t u_phi + spatial_part = 0
    const A = g_00;
    const B = 2 * g_03 * u_phi;
    const C = spatial_part;
    
    const discriminant = B * B - 4 * A * C;
    const u_t = discriminant >= 0 
      ? (-B + Math.sqrt(discriminant)) / (2 * A)
      : 0;
    
    return { u_t, u_r, u_theta, u_phi };
  }

  /**
   * Retrograde ISCO (for prograde, multiply by -1)
   */
  iscoRetrograde() {
    const M = this.M;
    const q = this.a / M;
    
    const term1 = (1 - q * q) ** (1/3);
    const term2 = -(1 + q) ** (1/3) - (1 - q) ** (1/3);
    
    const Z1 = 1 + term1 * term2;
    const Z2 = Math.sqrt(Z1 * Z1 - 4 * q * q);
    
    return M * (3 + Z2 + Math.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)));
  }

  /**
   * Photon sphere radius (unstable circular orbits for photons)
   * For Kerr: r_ph ≈ 2M [1 + cos(2/3 * arccos(±a/M))]  (prograde/retrograde)
   */
  photonSphereRadius(prograde = true) {
    const M = this.M;
    const a = this.a;
    
    // Approximate formula (exact requires solving cubic)
    const sign = prograde ? 1 : -1;
    const q = Math.min(Math.abs(a) / M, 0.998);
    
    const angle = Math.acos(sign * q);
    const factor = Math.cos((2/3) * angle);
    
    return 2 * M * (1 + factor);
  }

  /**
   * Get spin parameter as percentage (0-100)
   */
  getSpinPercentage() {
    return (this.a / (0.998 * this.M)) * 100;
  }

  /**
   * Check if a point is inside the event horizon
   */
  isInsideHorizon(r) {
    return r < this.r_plus;
  }

  /**
   * Check if a point is in the ergosphere
   */
  isInErgosphere(r, theta) {
    return r < this.ergosphereRadius(theta) && r > this.r_plus;
  }
}

// Export for use in both Node.js and browser
if (typeof module !== 'undefined' && module.exports) {
  module.exports = KerrMetric;
}
