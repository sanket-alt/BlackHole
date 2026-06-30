/**
 * GRMHDDataLoader.js - Load and manage GRMHD simulation data
 * Handles 3D grids and volumetric rendering
 */

class GRMHDDataLoader {
  constructor() {
    this.gridData = null;
    this.dimensions = [0, 0, 0];
    this.metadata = {};
    this.texture3D = null;
  }

  /**
   * Parse binary GRMHD data (assumes float32 array)
   * Format: [density, pressure, Bx, By, Bz] per cell
   */
  loadBinaryGrid(buffer, nx, ny, nz) {
    const floats = new Float32Array(buffer);
    const expectedLength = nx * ny * nz * 5;
    
    if (floats.length !== expectedLength) {
      console.warn(`Expected ${expectedLength} floats, got ${floats.length}`);
    }
    
    this.gridData = floats;
    this.dimensions = [nx, ny, nz];
    
    return {
      data: floats,
      dimensions: [nx, ny, nz],
      metadata: this.metadata,
    };
  }

  /**
   * Load from File input
   */
  async loadFromFile(file) {
    const buffer = await file.arrayBuffer();
    
    if (file.name.endsWith('.bin')) {
      // Parse binary with assumed grid dimensions from filename or metadata
      // Example: "grmhd_128x128x128.bin"
      const match = file.name.match(/(\d+)x(\d+)x(\d+)/);
      if (match) {
        const [, nx, ny, nz] = match.map(Number);
        return this.loadBinaryGrid(buffer, nx, ny, nz);
      }
    } else if (file.name.endsWith('.csv')) {
      return this.parseCSV(buffer);
    } else if (file.name.endsWith('.json')) {
      return this.parseJSON(buffer);
    }
    
    throw new Error('Unsupported file format');
  }

  /**
   * Create WebGL 3D texture from grid data
   */
  createWebGL3DTexture(gl, channel = 0) {
    const [nx, ny, nz] = this.dimensions;
    
    // Extract single channel data
    const channelData = new Float32Array(nx * ny * nz);
    for (let i = 0; i < nx * ny * nz; i++) {
      channelData[i] = this.gridData[i * 5 + channel];
    }
    
    // Create 3D texture
    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_3D, texture);
    
    // Configure texture
    gl.texImage3D(
      gl.TEXTURE_3D,
      0,
      gl.R32F,
      nx, ny, nz,
      0,
      gl.RED,
      gl.FLOAT,
      channelData
    );
    
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
    
    this.texture3D = texture;
    return texture;
  }

  /**
   * Get grid statistics
   */
  getStatistics() {
    if (!this.gridData) return null;
    
    const [nx, ny, nz] = this.dimensions;
    const numCells = nx * ny * nz;
    
    const stats = {
      density: { min: Infinity, max: -Infinity, mean: 0 },
      pressure: { min: Infinity, max: -Infinity, mean: 0 },
      B_field: { min: Infinity, max: -Infinity, mean: 0 },
    };
    
    let densitySum = 0, pressureSum = 0, bSum = 0;
    
    for (let i = 0; i < numCells; i++) {
      const density = this.gridData[i * 5];
      const pressure = this.gridData[i * 5 + 1];
      const bx = this.gridData[i * 5 + 2];
      const by = this.gridData[i * 5 + 3];
      const bz = this.gridData[i * 5 + 4];
      const bmag = Math.sqrt(bx * bx + by * by + bz * bz);
      
      // Density stats
      stats.density.min = Math.min(stats.density.min, density);
      stats.density.max = Math.max(stats.density.max, density);
      densitySum += density;
      
      // Pressure stats
      stats.pressure.min = Math.min(stats.pressure.min, pressure);
      stats.pressure.max = Math.max(stats.pressure.max, pressure);
      pressureSum += pressure;
      
      // B-field stats
      stats.B_field.min = Math.min(stats.B_field.min, bmag);
      stats.B_field.max = Math.max(stats.B_field.max, bmag);
      bSum += bmag;
    }
    
    stats.density.mean = densitySum / numCells;
    stats.pressure.mean = pressureSum / numCells;
    stats.B_field.mean = bSum / numCells;
    
    return stats;
  }

  parseCSV(buffer) {
    // CSV parsing implementation
    return null;
  }

  parseJSON(buffer) {
    // JSON parsing implementation
    const text = new TextDecoder().decode(buffer);
    const data = JSON.parse(text);
    return data;
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = GRMHDDataLoader;
}
