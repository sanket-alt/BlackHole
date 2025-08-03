if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

function Observer() {
    this.position = new THREE.Vector3(10,0,0);
    this.velocity = new THREE.Vector3(0,1,0);
    this.orientation = new THREE.Matrix3();
    this.time = 0.0;
}

Observer.prototype.orbitalFrame = function() {

    var orbital_y = (new THREE.Vector3())
        .subVectors(observer.velocity.clone().normalize().multiplyScalar(4.0),
            observer.position).normalize();

    var orbital_z = (new THREE.Vector3())
        .crossVectors(observer.position, orbital_y).normalize();
    var orbital_x = (new THREE.Vector3()).crossVectors(orbital_y, orbital_z);


    return (new THREE.Matrix4()).makeBasis(
        orbital_x,
        orbital_y,
        orbital_z
    ).linearPart();
};

Observer.prototype.move = function(dt) {

    dt *= shader.parameters.time_scale;

    var r;
    var v = 0;

    if (shader.parameters.observer.motion) {

        r = shader.parameters.observer.distance;
        v =  1.0 / Math.sqrt(2.0*(r-1.0));
        var ang_vel = v / r;
        var angle = this.time * ang_vel;

        var s = Math.sin(angle), c = Math.cos(angle);

        this.position.set(c*r, s*r, 0);
        this.velocity.set(-s*v, c*v, 0);

        var alpha = degToRad(shader.parameters.observer.orbital_inclination);
        var orbit_coords = (new THREE.Matrix4()).makeRotationY(alpha);

        this.position.applyMatrix4(orbit_coords);
        this.velocity.applyMatrix4(orbit_coords);
    }
    else {
        r = this.position.length();
    }

    if (shader.parameters.gravitational_time_dilation) {
        dt = Math.sqrt((dt*dt * (1.0 - v*v)) / (1-1.0/r));
    }

    this.time += dt;
};

var container, stats;
var camera, scene, renderer, cameraControls, shader = null;
var observer = new Observer();

function Shader(mustacheTemplate) {

    this.parameters = {
        n_steps: 100,
        quality: 'medium',
        accretion_disk: true,
        noise_scale: 3.5,
        noise_speed: 0.15,
        planet: {
            enabled: true,
            distance: 18.0,
            radius: 0.4
        },
  
        hawking_radiation: {
            enabled: false,

            log_mass_solar: 1, 
            mode: 'Energy Spectrum',
            evaporation_speed: 1.0,
        },
        magnetic_field: {
            enabled: true,
            strength: 10.0,
            field_lines: true,
            field_lines_density: 40, 
            field_color: [1.0, 0.5, 0.0], 
            dipole_tilt: 0.0, // Tilt angle of magnetic dipole in degrees
            show_plasma_effects: true
        },
        lorentz_contraction: true,
        gravitational_time_dilation: true,
        aberration: true,
        beaming: true,
        doppler_shift: true,
        light_travel_time: true,
        time_scale: 1.0,
        observer: {
            motion: false,
            distance: 14.0,
            orbital_inclination: -10
        },

        planetEnabled: function() {
            return this.planet.enabled && this.quality !== 'fast';
        },

        observerMotion: function() {
            return this.observer.motion;
        },
        // New function for template
        hawkingRadiationEnabled: function() {
            return this.hawking_radiation.enabled;
        },
        
        magneticFieldEnabled: function() {
            return this.magnetic_field.enabled;
        }
    };
    var that = this;
    this.needsUpdate = false;

    this.hasMovingParts = function() {
        return this.parameters.planet.enabled || this.parameters.observer.motion || this.parameters.accretion_disk || this.parameters.hawking_radiation.enabled || this.parameters.magnetic_field.enabled;
    };

    this.compile = function() {
        return Mustache.render(mustacheTemplate, that.parameters);
    };
}

function degToRad(a) { return Math.PI * a / 180.0; }

(function(){
    var textures = {};

    function whenLoaded() {
        init(textures);
        $('#loader').hide();
        $('.initially-hidden').removeClass('initially-hidden');
        animate();
    }

    function checkLoaded() {
        if (shader === null) return;
        for (var key in textures) if (textures[key] === null) return;
        whenLoaded();
    }

    SHADER_LOADER.load(function(shaders) {
        shader = new Shader(shaders.raytracer.fragment);
        checkLoaded();
    });

    var texLoader = new THREE.TextureLoader();
    function loadTexture(symbol, filename, interpolation) {
        textures[symbol] = null;
        texLoader.load(filename, function(tex) {
            tex.magFilter = interpolation;
            tex.minFilter = interpolation;
            textures[symbol] = tex;
            checkLoaded();
        });
    }

    loadTexture('galaxy', 'img/milkyway.jpg', THREE.NearestFilter);
    loadTexture('spectra', 'img/spectra.png', THREE.LinearFilter);
    loadTexture('moon', 'img/beach-ball.png', THREE.LinearFilter);
    loadTexture('stars', 'img/stars.png', THREE.LinearFilter);
})();

var updateUniforms;

// --- Physics Constants ---
const G = 6.67430e-11; // Gravitational constant
const H_BAR = 1.0545718e-34; // Reduced Planck constant
const C = 299792458.0; // Speed of light
const K_B = 1.380649e-23; // Boltzmann constant
const STEFAN_BOLTZMANN = 5.670374e-8;
const M_SOLAR = 1.989e30; // Mass of the sun in kg

var black_hole_mass_kg = M_SOLAR; 

function init(textures) {

    container = document.createElement( 'div' );
    document.body.appendChild( container );

    scene = new THREE.Scene();

    var geometry = new THREE.PlaneBufferGeometry( 2, 2 );

    var uniforms = {
        time: { type: "f", value: 0 },
        resolution: { type: "v2", value: new THREE.Vector2() },
        cam_pos: { type: "v3", value: new THREE.Vector3() },
        cam_x: { type: "v3", value: new THREE.Vector3() },
        cam_y: { type: "v3", value: new THREE.Vector3() },
        cam_z: { type: "v3", value: new THREE.Vector3() },
        cam_vel: { type: "v3", value: new THREE.Vector3() },

        planet_distance: { type: "f" },
        planet_radius: { type: "f" },

        noise_scale: { type: "f" },
        noise_speed: { type: "f" },

        // --- HAWKING RADIATION UNIFORMS ---
        hawking_mode: { type: "i", value: 0 },
        hawking_temperature: { type: "f", value: 0.0 },
        bh_radius: { type: "f", value: 1.0 }, // Schwarzschild radius in simulation units

        magnetic_strength: { type: "f", value: 1.0 },
        magnetic_dipole_tilt: { type: "f", value: 0.0 },
        field_lines_density: { type: "f", value: 50.0 },
        field_color: { type: "v3", value: new THREE.Vector3(1.0, 0.5, 0.0) },
        show_plasma_effects: { type: "i", value: 0 },

        star_texture: { type: "t", value: textures.stars },
        galaxy_texture: { type: "t", value: textures.galaxy },
        planet_texture: { type: "t", value: textures.moon },
        spectrum_texture: { type: "t", value: textures.spectra }
    };

    updateUniforms = function(dt) {
        uniforms.planet_distance.value = shader.parameters.planet.distance;
        uniforms.planet_radius.value = shader.parameters.planet.radius;

        uniforms.noise_scale.value = shader.parameters.noise_scale;
        uniforms.noise_speed.value = shader.parameters.noise_speed;

        uniforms.resolution.value.x = renderer.domElement.width;
        uniforms.resolution.value.y = renderer.domElement.height;

        uniforms.time.value = observer.time;
        uniforms.cam_pos.value = observer.position;

        var e = observer.orientation.elements;

        uniforms.cam_x.value.set(e[0], e[1], e[2]);
        uniforms.cam_y.value.set(e[3], e[4], e[5]);
        uniforms.cam_z.value.set(e[6], e[7], e[8]);

        function setVec(target, value) {
            uniforms[target].value.set(value.x, value.y, value.z);
        }

        setVec('cam_pos', observer.position);
        setVec('cam_vel', observer.velocity);

        // --- HAWKING RADIATION LOGIC ---
        if (shader.parameters.hawking_radiation.enabled) {
            const p = shader.parameters.hawking_radiation;

            // Mode selection
            const modes = ['Quantum Source', 'Energy Spectrum', 'Evaporation'];
            uniforms.hawking_mode.value = modes.indexOf(p.mode);

            // Evaporation logic
            if (p.mode === 'Evaporation' && black_hole_mass_kg > 1e5) { // Stop evaporation at a small mass
                // Power radiated by a black hole (Stefan-Boltzmann for BHs)
                // P = (h_bar * c^6) / (15360 * pi * G^2 * M^2)
                const numerator = H_BAR * Math.pow(C, 6);
                const denominator = 15360 * Math.PI * Math.pow(G, 2) * Math.pow(black_hole_mass_kg, 2);
                const power = numerator / denominator;

                const mass_loss_rate = power / Math.pow(C, 2);

                black_hole_mass_kg -= mass_loss_rate * dt * p.evaporation_speed;

                p.log_mass_solar = Math.log10(black_hole_mass_kg / M_SOLAR);
    
                for (var i in gui.__folders['Hawking Radiation'].__controllers) {
                    gui.__folders['Hawking Radiation'].__controllers[i].updateDisplay();
                }
            } else {
       
                 black_hole_mass_kg = Math.pow(10, p.log_mass_solar) * M_SOLAR;
            }

 
            black_hole_mass_kg = Math.max(1e5, black_hole_mass_kg);


            // Calculate Hawking Temperature: T = (h_bar * c^3) / (8 * pi * G * M * k_B)
            const temp_numerator = H_BAR * Math.pow(C, 3);
            const temp_denominator = 8 * Math.PI * G * black_hole_mass_kg * K_B;
            uniforms.hawking_temperature.value = temp_numerator / temp_denominator;

            const base_radius_for_1_solar_mass = 1.0; 
            const base_mass_kg = M_SOLAR;

            // Schwarzschild Radius: Rs = 2GM/c^2
            const real_radius_m = (2 * G * black_hole_mass_kg) / Math.pow(C, 2);
            const base_real_radius_m = (2 * G * base_mass_kg) / Math.pow(C, 2);

            uniforms.bh_radius.value = (real_radius_m / base_real_radius_m) * base_radius_for_1_solar_mass;
        }

        if (shader.parameters.magnetic_field.enabled) {
            const mf = shader.parameters.magnetic_field;
            
            uniforms.magnetic_strength.value = mf.strength;
            uniforms.magnetic_dipole_tilt.value = degToRad(mf.dipole_tilt);
            uniforms.field_lines_density.value = mf.field_lines_density;
            uniforms.field_color.value.set(mf.field_color[0], mf.field_color[1], mf.field_color[2]);
            uniforms.show_plasma_effects.value = mf.show_plasma_effects ? 1 : 0;
        }
    };

    var material = new THREE.ShaderMaterial( {
        uniforms: uniforms,
        vertexShader: document.getElementById('vertex-shader').textContent,
    });

    scene.updateShader = function() {
        material.fragmentShader = shader.compile();
        material.needsUpdate = true;
        shader.needsUpdate = true;
    };

    scene.updateShader();

    var mesh = new THREE.Mesh( geometry, material );
    scene.add( mesh );

    renderer = new THREE.WebGLRenderer();
    renderer.setPixelRatio( window.devicePixelRatio );
    container.appendChild( renderer.domElement );

    stats = new Stats();
    stats.domElement.style.position = 'absolute';
    stats.domElement.style.top = '0px';
    container.appendChild( stats.domElement );
    $(stats.domElement).addClass('hidden-phone');

    camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 1, 80000 );
    initializeCamera(camera);

    cameraControls = new THREE.OrbitControls( camera, renderer.domElement );
    cameraControls.target.set( 0, 0, 0 );
    cameraControls.addEventListener( 'change', updateCamera );
    updateCamera();

    onWindowResize();

    window.addEventListener( 'resize', onWindowResize, false );

    setupGUI();
}

var gui; // Make gui global to access it later

function setupGUI() {

    var hint = $('#hint-text');
    var p = shader.parameters;

    function updateShader() {
        hint.hide();
        scene.updateShader();
    }

    gui = new dat.GUI();

    gui.add(p, 'quality', ['fast', 'medium', 'high']).onChange(function (value) {
        $('.planet-controls').show();
        switch(value) {
        case 'fast':
            p.n_steps = 40;
            $('.planet-controls').hide();
            break;
        case 'medium':
            p.n_steps = 100;
            break;
        case 'high':
            p.n_steps = 200;
            break;
        }

        updateShader();
    });
    gui.add(p, 'accretion_disk').onChange(updateShader);

    var diskFolder = gui.addFolder('Accretion Disk');
    diskFolder.add(p, 'noise_scale', 0.1, 10.0).name('Noise Scale');
    diskFolder.add(p, 'noise_speed', 0.0, 1.0).name('Animation Speed');
    diskFolder.open();

    // --- HAWKING RADIATION GUI ---
    var hawkingFolder = gui.addFolder('Hawking Radiation');
    hawkingFolder.add(p.hawking_radiation, 'enabled').name('Enable Simulation').onChange(updateShader);
    hawkingFolder.add(p.hawking_radiation, 'mode', ['Quantum Source', 'Energy Spectrum', 'Evaporation']).name('Visualization Mode');
    // Using log scale for mass, as it spans many orders of magnitude
    hawkingFolder.add(p.hawking_radiation, 'log_mass_solar', -15, 10).name('log10(Mass / M☉)');
    hawkingFolder.add(p.hawking_radiation, 'evaporation_speed', 1e15, 1e25).name('Evaporation Speed');
    hawkingFolder.open();

    // --- MAGNETIC FIELD GUI ---
    var magneticFolder = gui.addFolder('Magnetic Field');
    magneticFolder.add(p.magnetic_field, 'enabled').name('Enable Magnetic Field').onChange(function(enabled) {
        updateShader();
        updateMagneticLegendVisibility();
    });
    magneticFolder.add(p.magnetic_field, 'strength', 0.1, 10.0).name('Field Strength');
    magneticFolder.add(p.magnetic_field, 'field_lines').name('Show Field Lines').onChange(function(enabled) {
        updateShader();
        updateMagneticLegendVisibility();
    });
    magneticFolder.add(p.magnetic_field, 'field_lines_density', 10, 100).name('Field Lines Density');
    magneticFolder.add(p.magnetic_field, 'dipole_tilt', -90, 90).name('Dipole Tilt (degrees)');
    magneticFolder.addColor(p.magnetic_field, 'field_color').name('Field Color');
    magneticFolder.add(p.magnetic_field, 'show_plasma_effects').name('Plasma Effects').onChange(updateShader);
    magneticFolder.open();

    function updateMagneticLegendVisibility() {
        var legend = $('#magnetic-legend');
        if (p.magnetic_field.enabled && p.magnetic_field.field_lines) {
            legend.show();
        } else {
            legend.hide();
        }
    }
    updateMagneticLegendVisibility();


    var folder = gui.addFolder('Observer');
    folder.add(p.observer, 'distance').min(1.5).max(30).onChange(updateCamera);
    folder.open();

    folder = gui.addFolder('Planet');
    folder.add(p.planet, 'enabled').onChange(function(enabled) {
        updateShader();
        var controls = $('.indirect-planet-controls').show();
        if (enabled) controls.show();
        else controls.hide();
    });
    folder.add(p.planet, 'distance').min(1.5);
    folder.add(p.planet, 'radius').min(0.01).max(2.0);
    $(folder.domElement).addClass('planet-controls');

    function setGuiRowClass(guiEl, klass) {
        $(guiEl.domElement).parent().parent().addClass(klass);
    }

    folder = gui.addFolder('Relativistic effects');
    folder.add(p, 'aberration').onChange(updateShader);
    folder.add(p, 'beaming').onChange(updateShader);
    folder.add(p, 'doppler_shift').onChange(updateShader);
    setGuiRowClass(
        folder.add(p, 'gravitational_time_dilation').onChange(updateShader),
        'planet-controls indirect-planet-controls');
    setGuiRowClass(
        folder.add(p, 'lorentz_contraction').onChange(updateShader),
        'planet-controls indirect-planet-controls');

    folder.open();

    folder = gui.addFolder('Time');
    folder.add(p, 'light_travel_time').onChange(updateShader);
    folder.add(p, 'time_scale').min(0);

}

function onWindowResize( event ) {
    renderer.setSize( window.innerWidth, window.innerHeight );
}

function initializeCamera(camera) {

    var pitchAngle = 3.0, yawAngle = 0.0;

    camera.matrixWorldInverse.makeRotationX(degToRad(-pitchAngle));
    camera.matrixWorldInverse.multiply(new THREE.Matrix4().makeRotationY(degToRad(-yawAngle)));

    var m = camera.matrixWorldInverse.elements;

    camera.position.set(m[2], m[6], m[10]);
}

function updateCamera( event ) {

    var zoom_dist = camera.position.length();
    var m = camera.matrixWorldInverse.elements;
    var camera_matrix;

    if (shader.parameters.observer.motion) {
        camera_matrix = new THREE.Matrix3();
    }
    else {
        camera_matrix = observer.orientation;
    }

    camera_matrix.set(
        m[0], m[1], m[2],
        m[8], m[9], m[10],
        m[4], m[5], m[6]
    );

    if (shader.parameters.observer.motion) {

        observer.orientation = observer.orbitalFrame().multiply(camera_matrix);

    } else {

        var p = new THREE.Vector3(
            camera_matrix.elements[6],
            camera_matrix.elements[7],
            camera_matrix.elements[8]);

        var dist = shader.parameters.observer.distance;
        observer.position.set(-p.x*dist, -p.y*dist, -p.z*dist);
        observer.velocity.set(0,0,0);
    }
}

function frobeniusDistance(matrix1, matrix2) {
    var sum = 0.0;
    for (var i=0; i < matrix1.elements.length; i++) {
        var diff = matrix1.elements[i] - matrix2.elements[i];
        sum += diff*diff;
    }
    return Math.sqrt(sum);
}

function animate() {
    requestAnimationFrame( animate );

    camera.updateMatrixWorld();
    camera.matrixWorldInverse.getInverse( camera.matrixWorld );

    if (shader.needsUpdate || shader.hasMovingParts() ||
        frobeniusDistance(camera.matrixWorldInverse, lastCameraMat) > 1e-10) {

        shader.needsUpdate = false;
        render();
        lastCameraMat = camera.matrixWorldInverse.clone();
    }
    stats.update();
}

var lastCameraMat = new THREE.Matrix4().identity();

var getFrameDuration = (function() {
    var lastTimestamp = new Date().getTime();
    return function() {
        var timestamp = new Date().getTime();
        var diff = (timestamp - lastTimestamp) / 1000.0;
        lastTimestamp = timestamp;
        return diff;
    };
})();

function render() {
    var dt = getFrameDuration();
    observer.move(dt);
    if (shader.parameters.observer.motion) updateCamera();
    updateUniforms(dt); 
    renderer.render( scene, camera );
}
