# Kerr Metric Black Hole Simulator - Research Grade Upgrade

## Overview

This project has been upgraded from a basic Schwarzschild ray-tracer to a **research-grade Kerr metric simulator** with support for spinning black holes, magnetic field topology, and scientific data export.

## New Physics Engine

### Core Components

1. **KerrMetric.js** - Boyer-Lindquist coordinate system
   - Event horizon and ergosphere calculations
   - ISCO (Innermost Stable Circular Orbit) computation
   - Metric tensor and Christoffel symbols
   - Spin parameter validation (0 ≤ a < 0.998M)

2. **GeodesicIntegrator.js** - RK4 null geodesic integration
   - Solves: `du^μ/dλ + Γ^μ_νρ u^ν u^ρ = 0`
   - Adaptive step sizing near horizons
   - Null geodesic normalization: `g_μν u^μ u^ν = 0`
   - Trajectory validation

3. **MagneticFieldTopology.js** - RK4 field line tracing
   - Poloidal (meridional) magnetic fields
   - Toroidal (azimuthal) magnetic fields
   - Dipole field configurations with tilt
   - Field line visualization vertices

4. **GRMHDDataLoader.js** - 3D grid ingestion
   - Binary and JSON format support
   - WebGL 3D texture creation
   - Grid statistics (density, pressure, B-field)

5. **ExportUtils.js** - Research data export
   - Radiance HDR format (human-readable)
   - FITS format (astronomy standard) with metadata
   - JSON simulation state export

## Mathematical Framework

### Kerr Metric (Boyer-Lindquist Coordinates)

```
ds² = -(1 - 2Mr/Σ)dt² - (4Mar sin²θ/Σ)dt dφ + (Σ/Δ)dr² 
      + Σ dθ² + sin²θ(r² + a² + 2Ma²r sin²θ/Σ)dφ²

where:
  Σ = r² + a² cos²θ
  Δ = r² - 2Mr + a²
  M = black hole mass
  a = spin parameter (0 ≤ a < M)
```

### Key Derived Quantities

- **Event Horizon**: r₊ = M + √(M² - a²)
- **Ergosphere**: r_ergo(θ) = M + √(M² - a² cos²θ)
- **ISCO (Prograde)**: r_ISCO = M[3 + Z₂ - √((3 - Z₁)(3 + Z₁ + 2Z₂))]
- **Photon Sphere**: r_ph ≈ 2M[1 + cos(2/3 arccos(±a/M))]
- **Frame Dragging**: Ω_F = 2ar₊/(r₊² + a²)²

## Usage

### Loading Physics Engines

```html
<script src="physics/kerr-metric.js"></script>
<script src="physics/geodesic-integrator.js"></script>
<script src="physics/magnetic-field-topology.js"></script>
<script src="constants.js"></script>
<script src="kerr-control-panel.js"></script>
```

### Creating a Kerr Simulation

```javascript
// Initialize Kerr metric
const metric = new KerrMetric(
  1.0,      // mass (solar masses)
  0.7       // spin parameter a
);

// Initialize geodesic integrator
const integrator = new GeodesicIntegrator(metric);

// Setup observer ray
const observerState = {
  r: 30,
  theta: Math.PI / 3,
  phi: 0,
  u_r: -0.1,
  u_theta: 0,
  u_phi: 0.05
};

// Trace ray
const trajectory = integrator.integrateRay(observerState, 100);
```

### Tracing Magnetic Field Lines

```javascript
const fieldTopology = new MagneticFieldTopology(metric);

const fieldLines = fieldTopology.traceFieldLineGrid(
  32,       // number of lines
  1.0,      // field strength
  0.5,      // poloidal/toroidal mix (0.5 = equal)
  0.0       // dipole tilt angle
);

const vertices = fieldTopology.getFieldLineVertices(fieldLines);
```

### Exporting Results

```javascript
// Capture render target
const pixels = ExportUtils.captureRenderTarget(renderer, 1920, 1080);

// Export as FITS
ExportUtils.exportToFITS(pixels, 1920, 1080, {
  spinParameter: 0.7,
  blackHoleMass: 1.0,
  inclination: 45,
  rayMarchSteps: 200
}, 'kerr-bh.fits');

// Export as Radiance HDR
ExportUtils.exportToHDR(pixels, 1920, 1080, 'kerr-bh.hdr');
```

## Validation & Testing

### Schwarzschild Limit (a → 0)
When spin parameter a = 0, the metric reduces to Schwarzschild:
- Verified: r₊ = 2M (no Cauchy horizon)
- Verified: Photon sphere at r_ph = 3M
- Verified: ISCO at r_ISCO = 6M

### Extremal Kerr (a → M)
At maximum spin a = 0.998M:
- Event horizon approaches: r₊ → M
- Ergosphere deformation: Notable compression at poles
- ISCO drastically reduced

### Null Geodesic Constraint
Code validates `|g_μν u^μ u^ν| < 1e-6` throughout integration

## Physics Constants

```javascript
const PhysicsConstants = {
  G: 6.67430e-11,          // Gravitational constant
  c: 299792458.0,          // Speed of light
  hbar: 1.0545718e-34,     // Reduced Planck constant
  M_solar: 1.989e30,       // Solar mass (kg)
  // ... see constants.js for more
};
```

## Next Steps (Planned)

1. **WebGPU Compute Shaders**: Offload geodesic integration to GPU
2. **Accretion Disk Physics**: Temperature/density from GRMHD grids
3. **Relativistic Beaming**: Lorentz transformation effects
4. **Hawking Radiation**: Quantum effects near horizon
5. **Performance Optimization**: Adaptive resolution, framebuffer accumulation

## References

- Kerr, R. P. (1963). "Gravitational field of a spinning mass as an example of algebraically special metrics"
- Newman, E. T.; Janis, A. I. (1965). "Note on the Kerr spinning-particle metric"
- Thorne, K. S. (1974). "Black Holes and Time Warps"
- Event Horizon Telescope Collaboration (2019): M87 black hole image analysis

## Contributing

Contributions are welcome! Please ensure:
- Physics accuracy (match published papers)
- Code comments on complex calculations
- Validation against Schwarzschild limit
- Performance profiling for GPU optimization

## License

See LICENSE file for details.
