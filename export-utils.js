/**
 * ExportUtils.js - Export rendered data to research-grade formats
 * Supports HDR and FITS with metadata
 */

class ExportUtils {
  /**
   * Extract raw float32 pixels from WebGL framebuffer
   */
  static captureRenderTarget(renderer, targetWidth, targetHeight) {
    const gl = renderer.getContext();
    const pixels = new Float32Array(targetWidth * targetHeight * 4);
    
    gl.readPixels(
      0, 0, targetWidth, targetHeight,
      gl.RGBA, gl.FLOAT, pixels
    );
    
    return pixels;
  }

  /**
   * Export to Radiance HDR format
   */
  static exportToHDR(pixels, width, height, filename) {
    let header = '#?RADIANCE\n';
    header += 'FORMAT=32-bit_rle_rgbe\n';
    header += `RESOLUTION X=${width} Y=${height}\n`;
    
    const rgbeData = this.floatToRGBE(pixels, width, height);
    
    const headerBytes = new TextEncoder().encode(header);
    const fullData = new Uint8Array(headerBytes.length + rgbeData.length);
    fullData.set(headerBytes);
    fullData.set(rgbeData, headerBytes.length);
    
    this.downloadFile(fullData, filename, 'image/vnd.radiance');
  }

  /**
   * Convert float RGBA to RGBE
   */
  static floatToRGBE(pixels, width, height) {
    const rgbeData = new Uint8Array(width * height * 4);
    
    for (let i = 0; i < pixels.length; i += 4) {
      const r = pixels[i];
      const g = pixels[i + 1];
      const b = pixels[i + 2];
      
      const maxc = Math.max(r, g, b);
      if (maxc < 1e-32) {
        rgbeData[i] = 0;
        rgbeData[i + 1] = 0;
        rgbeData[i + 2] = 0;
        rgbeData[i + 3] = 0;
      } else {
        const exp = Math.ceil(Math.log2(maxc));
        const scale = Math.pow(2, -exp + 8);
        
        rgbeData[i] = Math.round(r * scale);
        rgbeData[i + 1] = Math.round(g * scale);
        rgbeData[i + 2] = Math.round(b * scale);
        rgbeData[i + 3] = exp + 128;
      }
    }
    
    return rgbeData;
  }

  /**
   * Create minimal FITS header with physics metadata
   */
  static createFITSHeader(metadata) {
    const cards = [];
    
    cards.push({ keyword: 'SIMPLE', value: true, comment: 'Conforms to FITS standard' });
    cards.push({ keyword: 'BITPIX', value: -32, comment: '32-bit floating point' });
    cards.push({ keyword: 'NAXIS', value: 2, comment: 'Number of axes' });
    cards.push({ keyword: 'NAXIS1', value: metadata.width, comment: 'Image width' });
    cards.push({ keyword: 'NAXIS2', value: metadata.height, comment: 'Image height' });
    cards.push({ keyword: 'BH_SPIN', value: metadata.spinParameter, comment: 'Black hole spin parameter (a)' });
    cards.push({ keyword: 'BH_MASS', value: metadata.blackHoleMass, comment: 'Black hole mass (solar masses)' });
    cards.push({ keyword: 'INCLINATION', value: metadata.inclination, comment: 'Observer inclination angle (degrees)' });
    cards.push({ keyword: 'RAYMARCHING', value: metadata.rayMarchSteps, comment: 'Number of ray-marching steps' });
    cards.push({ keyword: 'TIMESTAMP', value: new Date().toISOString(), comment: 'Generation timestamp' });
    cards.push({ keyword: 'ORIGIN', value: 'Kerr Black Hole Ray Tracer', comment: 'Generation software' });
    cards.push({ keyword: 'END', value: '', comment: '' });
    
    return this.formatFITSCards(cards);
  }

  /**
   * Format FITS header cards (80 character fixed width)
   */
  static formatFITSCards(cards) {
    const headerBuffer = new Uint8Array(80 * Math.ceil(cards.length / 36) * 36);
    let offset = 0;
    
    for (const card of cards) {
      let line = '';
      
      if (card.keyword === 'END') {
        line = 'END' + ' '.repeat(77);
      } else {
        line = card.keyword.padEnd(8);
        line += '= ';
        
        if (typeof card.value === 'string') {
          line += `'${card.value.padEnd(8)}'`;
        } else {
          line += String(card.value).padStart(20);
        }
        
        if (card.comment) {
          line += ` / ${card.comment}`;
        }
        
        line = line.padEnd(80);
      }
      
      const bytes = new TextEncoder().encode(line);
      headerBuffer.set(bytes, offset);
      offset += 80;
    }
    
    return headerBuffer;
  }

  /**
   * Export to FITS format with data
   */
  static exportToFITS(pixels, width, height, metadata, filename) {
    const header = this.createFITSHeader({ width, height, ...metadata });
    
    const dataSize = width * height * 4 * 4;
    const totalSize = header.length + dataSize;
    const fitsData = new Uint8Array(totalSize);
    
    fitsData.set(header);
    fitsData.set(new Uint8Array(pixels.buffer), header.length);
    
    const paddedSize = Math.ceil(totalSize / 2880) * 2880;
    const paddedData = new Uint8Array(paddedSize);
    paddedData.set(fitsData);
    
    this.downloadFile(paddedData, filename, 'image/fits');
  }

  /**
   * Download file to user's computer
   */
  static downloadFile(data, filename, mimeType) {
    const blob = new Blob([data], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }

  /**
   * Export simulation state as JSON
   */
  static exportStateAsJSON(state, filename) {
    const data = {
      version: '1.0',
      timestamp: new Date().toISOString(),
      parameters: {
        spinParameter: state.spinParameter,
        blackHoleMass: state.blackHoleMass,
        inclination: state.inclination,
        rayMarchSteps: state.rayMarchSteps,
      },
    };
    
    const jsonString = JSON.stringify(data, null, 2);
    const blob = new Blob([jsonString], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = ExportUtils;
}
