/**
 * Physical and mathematical constants for Kerr metric calculations
 */

const PhysicsConstants = {
  // Fundamental constants (SI units)
  G: 6.67430e-11,           // Gravitational constant (m³ kg⁻¹ s⁻²)
  c: 299792458.0,           // Speed of light (m/s)
  hbar: 1.0545718e-34,      // Reduced Planck constant (J·s)
  kB: 1.380649e-23,         // Boltzmann constant (J/K)
  sigma_SB: 5.670374e-8,    // Stefan-Boltzmann constant (W m⁻² K⁻⁴)
  
  // Astronomical constants
  M_solar: 1.989e30,        // Solar mass (kg)
  AU: 1.496e11,             // Astronomical unit (m)
  pc: 3.086e16,             // Parsec (m)
  
  // Black hole parameters (natural units where c = G = 1)
  // For 1 solar mass black hole:
  schwarzschild_radius_1Msun: 2954.0,  // meters
  
  // Mathematical constants
  PI: Math.PI,
  TWO_PI: 2 * Math.PI,
  SQRT_2: Math.sqrt(2),
  
  /**
   * Convert solar masses to meters (c = G = 1)
   */
  solMassToMeters(M_solar_count) {
    return (2 * this.G * M_solar_count * this.M_solar) / (this.c * this.c);
  },
  
  /**
   * Convert meters to solar masses (c = G = 1)
   */
  metersToSolMass(r_meters) {
    return (r_meters * this.c * this.c) / (2 * this.G * this.M_solar);
  },
  
  /**
   * Hawking temperature in Kelvin
   * T_H = (ℏc³) / (8πkB G M)
   */
  hawkingTemperature(M_kg) {
    const numerator = this.hbar * Math.pow(this.c, 3);
    const denominator = 8 * Math.PI * this.kB * this.G * M_kg;
    return numerator / denominator;
  },
  
  /**
   * Hawking luminosity in Watts
   * L_H = (ℏc⁶) / (15360πG²M²)
   */
  hawkingLuminosity(M_kg) {
    const numerator = this.hbar * Math.pow(this.c, 6);
    const denominator = 15360 * Math.PI * Math.pow(this.G, 2) * Math.pow(M_kg, 2);
    return numerator / denominator;
  },
  
  /**
   * Black hole evaporation time in seconds
   * t_evap = (5120 π G² M₀³) / (ℏc⁴)
   */
  evaporationTime(M_kg) {
    const numerator = 5120 * Math.PI * Math.pow(this.G, 2) * Math.pow(M_kg, 3);
    const denominator = this.hbar * Math.pow(this.c, 4);
    return numerator / denominator;
  },
  
  /**
   * Schwarzschild radius in meters
   * r_s = 2GM/c²
   */
  schwarzschildRadius(M_kg) {
    return (2 * this.G * M_kg) / (this.c * this.c);
  },
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = PhysicsConstants;
}
