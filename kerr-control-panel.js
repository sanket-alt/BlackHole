/**
 * ControlPanel.js - UI for Kerr metric parameters
 * Manages spin, mass, observer inclination, and visualization overlays
 */

class KerrControlPanel {
  constructor(containerId, onUpdate) {
    this.container = document.getElementById(containerId);
    this.onUpdate = onUpdate;
    
    this.state = {
      spinParameter: 0.7,
      blackHoleMass: 1.0,
      inclination: 45,
      accretionDiskInner: 6.0,
      accretionDiskOuter: 40.0,
      showEventHorizon: true,
      showErgosphere: true,
      showISCO: true,
      magneticFieldEnabled: true,
      magneticFieldStrength: 1.0,
      poloidalVsToroidal: 0.5,
      rayMarchSteps: 200,
    };
    
    this.initializeUI();
  }

  initializeUI() {
    this.container.innerHTML = `
      <div class="kerr-control-panel">
        <h2>Kerr Black Hole Simulator</h2>
        
        <div class="panel-section">
          <h3>Physical Parameters</h3>
          <div class="control-group">
            <label>Spin Parameter (a): <span id="spin-value">0.700</span></label>
            <input type="range" id="spin-slider" min="0" max="0.998" step="0.01" value="0.7">
            <p class="help-text">0 = Schwarzschild | 0.998 = Extremal Kerr</p>
          </div>
          
          <div class="control-group">
            <label>Black Hole Mass (M☉): <span id="mass-value">1.0</span></label>
            <input type="range" id="mass-slider" min="0.1" max="10" step="0.1" value="1">
          </div>
          
          <div class="control-group">
            <label>Observer Inclination (°): <span id="inclination-value">45</span></label>
            <input type="range" id="inclination-slider" min="0" max="90" step="1" value="45">
          </div>
          
          <div class="info-box">
            <h4>Kerr Geometry</h4>
            <p>Event Horizon: r₊ = <span id="info-horizon">1.000</span></p>
            <p>ISCO: r_ISCO = <span id="info-isco">6.000</span></p>
            <p>Photon Sphere: r_ph ≈ <span id="info-photon">2.598</span></p>
            <p>Spin Energy: <span id="info-spin-energy">70%</span></p>
          </div>
        </div>
        
        <div class="panel-section">
          <h3>Overlays</h3>
          <div class="control-group">
            <label><input type="checkbox" id="show-horizon" checked> Event Horizon</label>
            <label><input type="checkbox" id="show-ergosphere" checked> Ergosphere</label>
            <label><input type="checkbox" id="show-isco" checked> ISCO</label>
          </div>
        </div>
        
        <div class="panel-section">
          <h3>Magnetic Field</h3>
          <div class="control-group">
            <label><input type="checkbox" id="magnetic-enabled" checked> Enable Magnetic Field</label>
          </div>
          
          <div class="control-group">
            <label>Field Strength: <span id="magnetic-strength-value">1.0</span></label>
            <input type="range" id="magnetic-strength" min="0.1" max="10" step="0.1" value="1">
          </div>
          
          <div class="control-group">
            <label>Poloidal ← → Toroidal: <span id="mix-value">50%</span></label>
            <input type="range" id="poloidal-toroidal" min="0" max="1" step="0.05" value="0.5">
          </div>
        </div>
        
        <div class="panel-section">
          <h3>Rendering</h3>
          <div class="control-group">
            <label>Ray March Steps: <span id="steps-value">200</span></label>
            <input type="range" id="ray-march-steps" min="50" max="500" step="10" value="200">
          </div>
        </div>
        
        <div class="panel-section">
          <h3>Export</h3>
          <button id="export-hdr">Export as HDR</button>
          <button id="export-fits">Export as FITS</button>
          <button id="export-state">Export State (JSON)</button>
        </div>
      </div>
    `;
    
    this.attachEventListeners();
  }

  attachEventListeners() {
    // Spin parameter
    document.getElementById('spin-slider').addEventListener('input', (e) => {
      this.state.spinParameter = parseFloat(e.target.value);
      document.getElementById('spin-value').textContent = e.target.value;
      this.updateInfo();
      this.onUpdate(this.state);
    });
    
    // Mass
    document.getElementById('mass-slider').addEventListener('input', (e) => {
      this.state.blackHoleMass = parseFloat(e.target.value);
      document.getElementById('mass-value').textContent = e.target.value;
      this.updateInfo();
      this.onUpdate(this.state);
    });
    
    // Inclination
    document.getElementById('inclination-slider').addEventListener('input', (e) => {
      this.state.inclination = parseFloat(e.target.value);
      document.getElementById('inclination-value').textContent = e.target.value;
      this.onUpdate(this.state);
    });
    
    // Overlays
    document.getElementById('show-horizon').addEventListener('change', (e) => {
      this.state.showEventHorizon = e.target.checked;
      this.onUpdate(this.state);
    });
    
    document.getElementById('show-ergosphere').addEventListener('change', (e) => {
      this.state.showErgosphere = e.target.checked;
      this.onUpdate(this.state);
    });
    
    document.getElementById('show-isco').addEventListener('change', (e) => {
      this.state.showISCO = e.target.checked;
      this.onUpdate(this.state);
    });
    
    // Magnetic field
    document.getElementById('magnetic-enabled').addEventListener('change', (e) => {
      this.state.magneticFieldEnabled = e.target.checked;
      this.onUpdate(this.state);
    });
    
    document.getElementById('magnetic-strength').addEventListener('input', (e) => {
      this.state.magneticFieldStrength = parseFloat(e.target.value);
      document.getElementById('magnetic-strength-value').textContent = e.target.value;
      this.onUpdate(this.state);
    });
    
    document.getElementById('poloidal-toroidal').addEventListener('input', (e) => {
      this.state.poloidalVsToroidal = parseFloat(e.target.value);
      document.getElementById('mix-value').textContent = (e.target.value * 100).toFixed(0) + '%';
      this.onUpdate(this.state);
    });
    
    // Ray march steps
    document.getElementById('ray-march-steps').addEventListener('input', (e) => {
      this.state.rayMarchSteps = parseInt(e.target.value);
      document.getElementById('steps-value').textContent = e.target.value;
      this.onUpdate(this.state);
    });
  }

  updateInfo() {
    const metric = new KerrMetric(this.state.blackHoleMass, this.state.spinParameter);
    
    document.getElementById('info-horizon').textContent = metric.r_plus.toFixed(3);
    document.getElementById('info-isco').textContent = metric.r_isco.toFixed(3);
    document.getElementById('info-photon').textContent = metric.photonSphereRadius().toFixed(3);
    document.getElementById('info-spin-energy').textContent = metric.getSpinPercentage().toFixed(0) + '%';
  }

  getState() {
    return { ...this.state };
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = KerrControlPanel;
}
