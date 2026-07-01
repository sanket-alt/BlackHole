if ( !Detector.webgl ) Detector.addGetWebGLMessage();

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

    this.orientation.set(
        orbital_x.x, orbital_y.x, orbital_z.x,
        orbital_x.y, orbital_y.y, orbital_z.y,
        orbital_x.z, orbital_y.z, orbital_z.z
    );
};

var observer = new Observer();

var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 1000 );
var renderer = new THREE.WebGLRenderer();
renderer.setSize( window.innerWidth, window.innerHeight );
renderer.setClearColor( 0x000000, 1 );
document.body.appendChild( renderer.domElement );

var clock = new THREE.Clock();
var time = 0;
var camera_moved = false;
var camera_move_timeout;

var plane = new THREE.PlaneGeometry( 2, 2 );

var uniforms = {
    time: { value: 0.0 },
    resolution: { value: new THREE.Vector2( window.innerWidth, window.innerHeight ) },
    cam_pos: { value: observer.position },
    cam_x: { value: new THREE.Vector3( 1, 0, 0 ) },
    cam_y: { value: new THREE.Vector3( 0, 1, 0 ) },
    cam_z: { value: new THREE.Vector3( 0, 0, 1 ) },
    cam_vel: { value: observer.velocity },
    planet_distance: { value: 20.0 },
    planet_radius: { value: 6.371 },
    noise_scale: { value: 0.1 },
    noise_speed: { value: 1.0 },
    galaxy_texture: { value: new THREE.TextureLoader().load( 'galaxy.png' ) },
    star_texture: { value: new THREE.TextureLoader().load( 'stars.png' ) },
    planet_texture: { value: new THREE.TextureLoader().load( 'planet.jpg' ) },
    spectrum_texture: { value: new THREE.TextureLoader().load( 'spectrum.png' ) },
    magnetic_strength: { value: 1.0 },
    magnetic_dipole_tilt: { value: 0.0 },
    field_lines_density: { value: 1.0 },
    field_color: { value: new THREE.Vector3( 0.0, 1.0, 1.0 ) },
    show_plasma_effects: { value: 1 },
    hawking_radiation_enabled: { value: false },
    hawking_mode: { value: 0 },
    hawking_temperature: { value: 1e21 },
    bh_radius: { value: 1.0 }
};

var material = new THREE.ShaderMaterial( {
    uniforms: uniforms,
    vertexShader: `void main() { gl_Position = vec4(position, 1.0); }`,
    fragmentShader: `uniform vec2 resolution;
void main() { gl_FragColor = vec4(gl_FragCoord.xy / resolution, 0.5, 1.0); }`
} );

var mesh = new THREE.Mesh( plane, material );
scene.add( mesh );

window.addEventListener( 'resize', function() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
    uniforms.resolution.value = new THREE.Vector2( window.innerWidth, window.innerHeight );
});

function animate() {
    requestAnimationFrame( animate );
    
    time += clock.getDelta();
    uniforms.time.value = time;
    
    observer.orbitalFrame();
    uniforms.cam_pos.value = observer.position;
    uniforms.cam_x.value.setFromMatrixColumn(observer.orientation, 0);
    uniforms.cam_y.value.setFromMatrixColumn(observer.orientation, 1);
    uniforms.cam_z.value.setFromMatrixColumn(observer.orientation, 2);
    
    renderer.render( scene, camera );
}

animate();

function degToRad(a) { return Math.PI * a / 180.0; }
