/**
 * ASEView.js v0.1.0
 * JavaScript library for molecular visualization
 * https://github.com/kangmg/aseview
 * (c) 2024 MIT License
 */
/**
 * Scene manager for Three.js rendering
 */
class SceneManager {
    constructor(container, options = {}) {
        this.container = typeof container === 'string'
            ? document.querySelector(container)
            : container;

        if (!this.container) {
            throw new Error('Container element not found');
        }

        this.options = {
            backgroundColor: options.backgroundColor || '#1f2937',
            antialias: options.antialias !== false,
            ...options
        };

        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.directionalLight = null;
        this.moleculeGroup = null;

        // Axis helper
        this.axisScene = null;
        this.axisCamera = null;

        this._init();
    }

    _init() {
        // Create renderer
        this.renderer = new THREE.WebGLRenderer({
            antialias: this.options.antialias,
            alpha: true,
            preserveDrawingBuffer: true
        });
        this.renderer.setPixelRatio(window.devicePixelRatio || 1);
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.outputEncoding = THREE.sRGBEncoding;
        this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
        this.renderer.toneMappingExposure = 1.0;
        this.container.appendChild(this.renderer.domElement);

        // Create scene
        this.scene = new THREE.Scene();
        this._updateBackgroundColor();

        // Create camera
        this._setupCamera();

        // Create controls
        this._setupControls();

        // Setup lighting
        this._setupLighting();

        // Setup axis helper
        this._setupAxisHelper();

        // Handle resize
        window.addEventListener('resize', () => this._handleResize());
    }

    _setupCamera() {
        const width = this.container.clientWidth || 1;
        const height = this.container.clientHeight || 1;
        const aspect = width / height;

        if (this.options.viewMode === 'Orthographic') {
            const frustumSize = 20;
            this.camera = new THREE.OrthographicCamera(
                (frustumSize * aspect) / -2,
                (frustumSize * aspect) / 2,
                frustumSize / 2,
                frustumSize / -2,
                -1e3,
                1000
            );
        } else {
            this.camera = new THREE.PerspectiveCamera(55, aspect, 0.1, 2000);
        }
        this.camera.position.set(0, 0, 10);
    }

    _setupControls() {
        if (this.controls) {
            this.controls.dispose();
        }

        if (this.options.rotationMode === 'Orbit' && THREE.OrbitControls) {
            this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        } else if (THREE.TrackballControls) {
            this.controls = new THREE.TrackballControls(this.camera, this.renderer.domElement);
            this.controls.rotateSpeed = 5.0;
            this.controls.zoomSpeed = 1.2;
            this.controls.panSpeed = 0.8;
        }

        if (this.controls) {
            this.controls.update();
        }
    }

    _setupLighting() {
        this.directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        this.directionalLight.position.set(10, 10, 10);
        this.scene.add(this.directionalLight);

        const ambientLight = new THREE.AmbientLight(0xffffff, 0.2);
        this.scene.add(ambientLight);
    }

    _setupAxisHelper() {
        this.axisScene = new THREE.Scene();

        const axisSize = 100;
        this.axisCamera = new THREE.OrthographicCamera(-axisSize, axisSize, axisSize, -axisSize, -100, 100);
        this.axisCamera.position.set(0, 0, 50);

        const colors = {
            x: 0xff4444,
            y: 0x44cc44,
            z: 0x4488ff
        };

        const axisLength = 40;
        const arrowHeadLength = 12;
        const arrowHeadWidth = 6;

        const xDir = new THREE.Vector3(1, 0, 0);
        const yDir = new THREE.Vector3(0, 1, 0);
        const zDir = new THREE.Vector3(0, 0, 1);
        const origin = new THREE.Vector3(0, 0, 0);

        const xArrow = new THREE.ArrowHelper(xDir, origin, axisLength, colors.x, arrowHeadLength, arrowHeadWidth);
        const yArrow = new THREE.ArrowHelper(yDir, origin, axisLength, colors.y, arrowHeadLength, arrowHeadWidth);
        const zArrow = new THREE.ArrowHelper(zDir, origin, axisLength, colors.z, arrowHeadLength, arrowHeadWidth);

        this.axisScene.add(xArrow);
        this.axisScene.add(yArrow);
        this.axisScene.add(zArrow);

        // Create text labels
        const createTextSprite = (text, color) => {
            const canvas = document.createElement('canvas');
            canvas.width = 64;
            canvas.height = 64;
            const ctx = canvas.getContext('2d');
            ctx.font = 'bold 48px Arial';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillStyle = '#' + color.toString(16).padStart(6, '0');
            ctx.fillText(text, 32, 32);

            const texture = new THREE.CanvasTexture(canvas);
            const material = new THREE.SpriteMaterial({ map: texture, transparent: true });
            const sprite = new THREE.Sprite(material);
            sprite.scale.set(16, 16, 1);
            return sprite;
        };

        const xLabel = createTextSprite('X', colors.x);
        const yLabel = createTextSprite('Y', colors.y);
        const zLabel = createTextSprite('Z', colors.z);

        xLabel.position.set(axisLength + 12, 0, 0);
        yLabel.position.set(0, axisLength + 12, 0);
        zLabel.position.set(0, 0, axisLength + 12);

        this.axisScene.add(xLabel);
        this.axisScene.add(yLabel);
        this.axisScene.add(zLabel);

        const axisAmbient = new THREE.AmbientLight(0xffffff, 1.0);
        this.axisScene.add(axisAmbient);
    }

    _updateBackgroundColor() {
        const color = new THREE.Color(this.options.backgroundColor);
        this.scene.background = color;
        this.renderer.setClearColor(color, 1);
    }

    _handleResize() {
        const width = this.container.clientWidth || 1;
        const height = this.container.clientHeight || 1;

        this.renderer.setSize(width, height);

        if (this.camera.isPerspectiveCamera) {
            this.camera.aspect = width / height;
        } else {
            const frustumSize = 20;
            const aspect = width / height;
            this.camera.left = (frustumSize * aspect) / -2;
            this.camera.right = (frustumSize * aspect) / 2;
            this.camera.top = frustumSize / 2;
            this.camera.bottom = frustumSize / -2;
        }
        this.camera.updateProjectionMatrix();

        if (this.controls && typeof this.controls.handleResize === 'function') {
            this.controls.handleResize();
        }
    }

    render() {
        if (this.controls) {
            this.controls.update();
        }

        // Update light position
        if (this.directionalLight && this.camera) {
            const offset = new THREE.Vector3(5, 8, 5);
            const lightPos = this.camera.position.clone().add(
                offset.applyQuaternion(this.camera.quaternion)
            );
            this.directionalLight.position.copy(lightPos);
        }

        this.renderer.render(this.scene, this.camera);

        // Render axis helper
        if (this.axisScene && this.axisCamera) {
            this.axisCamera.quaternion.copy(this.camera.quaternion);

            const currentViewport = this.renderer.getViewport(new THREE.Vector4());
            const axisViewportSize = 100;

            this.renderer.setViewport(10, 10, axisViewportSize, axisViewportSize);
            this.renderer.setScissor(10, 10, axisViewportSize, axisViewportSize);
            this.renderer.setScissorTest(true);
            this.renderer.clearDepth();
            this.renderer.render(this.axisScene, this.axisCamera);
            this.renderer.setScissorTest(false);
            this.renderer.setViewport(currentViewport);
        }
    }

    animate() {
        const loop = () => {
            requestAnimationFrame(loop);
            this.render();
        };
        loop();
    }

    fitToMolecule(group) {
        const box = new THREE.Box3().setFromObject(group);
        if (!isFinite(box.max.x)) return;

        const center = box.getCenter(new THREE.Vector3());
        const size = box.getSize(new THREE.Vector3());
        const maxDim = Math.max(size.x, size.y, size.z, 1);

        if (this.camera.isPerspectiveCamera) {
            const halfFov = THREE.MathUtils.degToRad(this.camera.fov / 2);
            const distance = maxDim / (2 * Math.tan(halfFov)) + maxDim;
            this.camera.position.copy(center.clone().add(new THREE.Vector3(distance, distance, distance)));
        } else {
            const aspect = this.container.clientWidth / this.container.clientHeight;
            this.camera.left = -maxDim * aspect;
            this.camera.right = maxDim * aspect;
            this.camera.top = maxDim;
            this.camera.bottom = -maxDim;
        }

        if (this.controls) {
            this.controls.target.copy(center);
            this.controls.update();
        }
        this.camera.lookAt(center);
        this.camera.updateProjectionMatrix();
    }

    dispose() {
        if (this.controls) {
            this.controls.dispose();
        }
        if (this.renderer) {
            this.renderer.dispose();
            this.container.removeChild(this.renderer.domElement);
        }
    }
}

/**
 * Atom information: colors and radii for common elements
 */
const atomInfo = {
    'H':  { color: 0xFFFFFF, radius: 0.31 },
    'He': { color: 0xD9FFFF, radius: 0.28 },
    'Li': { color: 0xCC80FF, radius: 1.28 },
    'Be': { color: 0xC2FF00, radius: 0.96 },
    'B':  { color: 0xFFB5B5, radius: 0.84 },
    'C':  { color: 0x909090, radius: 0.76 },
    'N':  { color: 0x3050F8, radius: 0.71 },
    'O':  { color: 0xFF0D0D, radius: 0.66 },
    'F':  { color: 0x90E050, radius: 0.57 },
    'Ne': { color: 0xB3E3F5, radius: 0.58 },
    'Na': { color: 0xAB5CF2, radius: 1.66 },
    'Mg': { color: 0x8AFF00, radius: 1.41 },
    'Al': { color: 0xBFA6A6, radius: 1.21 },
    'Si': { color: 0xF0C8A0, radius: 1.11 },
    'P':  { color: 0xFF8000, radius: 1.07 },
    'S':  { color: 0xFFFF30, radius: 1.05 },
    'Cl': { color: 0x1FF01F, radius: 1.02 },
    'Ar': { color: 0x80D1E3, radius: 1.06 },
    'K':  { color: 0x8F40D4, radius: 2.03 },
    'Ca': { color: 0x3DFF00, radius: 1.76 },
    'Sc': { color: 0xE6E6E6, radius: 1.70 },
    'Ti': { color: 0xBFC2C7, radius: 1.60 },
    'V':  { color: 0xA6A6AB, radius: 1.53 },
    'Cr': { color: 0x8A99C7, radius: 1.39 },
    'Mn': { color: 0x9C7AC7, radius: 1.39 },
    'Fe': { color: 0xE06633, radius: 1.32 },
    'Co': { color: 0xF090A0, radius: 1.26 },
    'Ni': { color: 0x50D050, radius: 1.24 },
    'Cu': { color: 0xC88033, radius: 1.32 },
    'Zn': { color: 0x7D80B0, radius: 1.22 },
    'Ga': { color: 0xC28F8F, radius: 1.22 },
    'Ge': { color: 0x668F8F, radius: 1.20 },
    'As': { color: 0xBD80E3, radius: 1.19 },
    'Se': { color: 0xFFA100, radius: 1.20 },
    'Br': { color: 0xA62929, radius: 1.20 },
    'Kr': { color: 0x5CB8D1, radius: 1.16 },
    'Rb': { color: 0x702EB0, radius: 2.20 },
    'Sr': { color: 0x00FF00, radius: 1.95 },
    'Y':  { color: 0x94FFFF, radius: 1.90 },
    'Zr': { color: 0x94E0E0, radius: 1.75 },
    'Nb': { color: 0x73C2C9, radius: 1.64 },
    'Mo': { color: 0x54B5B5, radius: 1.54 },
    'Tc': { color: 0x3B9E9E, radius: 1.47 },
    'Ru': { color: 0x248F8F, radius: 1.46 },
    'Rh': { color: 0x0A7D8C, radius: 1.42 },
    'Pd': { color: 0x006985, radius: 1.39 },
    'Ag': { color: 0xC0C0C0, radius: 1.45 },
    'Cd': { color: 0xFFD98F, radius: 1.44 },
    'In': { color: 0xA67573, radius: 1.42 },
    'Sn': { color: 0x668080, radius: 1.39 },
    'Sb': { color: 0x9E63B5, radius: 1.39 },
    'Te': { color: 0xD47A00, radius: 1.38 },
    'I':  { color: 0x940094, radius: 1.39 },
    'Xe': { color: 0x429EB0, radius: 1.40 },
    'Cs': { color: 0x57178F, radius: 2.44 },
    'Ba': { color: 0x00C900, radius: 2.15 },
    'La': { color: 0x70D4FF, radius: 2.07 },
    'Ce': { color: 0xFFFFC7, radius: 2.04 },
    'Pt': { color: 0xD0D0E0, radius: 1.36 },
    'Au': { color: 0xFFD123, radius: 1.36 },
    'Hg': { color: 0xB8B8D0, radius: 1.32 },
    'Tl': { color: 0xA6544D, radius: 1.45 },
    'Pb': { color: 0x575961, radius: 1.46 },
    'Bi': { color: 0x9E4FB5, radius: 1.48 },
    'Po': { color: 0xAB5C00, radius: 1.40 },
    'At': { color: 0x754F45, radius: 1.50 },
    'Rn': { color: 0x428296, radius: 1.50 },
    'Fr': { color: 0x420066, radius: 2.60 },
    'Ra': { color: 0x007D00, radius: 2.21 },
    'Ac': { color: 0x70ABFA, radius: 2.15 },
    'Th': { color: 0x00BAFF, radius: 2.06 },
    'Pa': { color: 0x00A1FF, radius: 2.00 },
    'U':  { color: 0x008FFF, radius: 1.96 },
    'default': { color: 0xFF1493, radius: 1.50 }
};

/**
 * Get atom info for a given symbol
 */
function getAtomInfo(symbol) {
    return atomInfo[symbol] || atomInfo['default'];
}

/**
 * Detect bonds between atoms based on covalent radii
 */
function detectBonds(positions, symbols, threshold = 1.5) {
    const bonds = [];
    const n = positions.length;

    for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
            const info1 = getAtomInfo(symbols[i]);
            const info2 = getAtomInfo(symbols[j]);

            const dx = positions[i][0] - positions[j][0];
            const dy = positions[i][1] - positions[j][1];
            const dz = positions[i][2] - positions[j][2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

            const cutoff = (info1.radius + info2.radius) * threshold;
            if (dist <= cutoff) {
                bonds.push([i, j]);
            }
        }
    }

    return bonds;
}

/**
 * Atom rendering with different styles
 */

/**
 * Create atom mesh with specified style
 */
function createAtom(position, symbol, scale = 1.0, style = 'cartoon', options = {}) {
    const info = getAtomInfo(symbol);
    const radius = info.radius * scale;
    const color = options.color !== undefined ? options.color : info.color;

    const pos = Array.isArray(position)
        ? new THREE.Vector3(position[0], position[1], position[2])
        : position;

    switch (style) {
        case 'cartoon':
            return createAtomCartoon(pos, radius, color);
        case 'glossy':
            return createAtomGlossy(pos, radius, color);
        case 'metallic':
            return createAtomMetallic(pos, radius, color);
        case 'neon':
            return createAtomNeon(pos, radius, color);
        case 'bubble':
            return createAtomBubble(pos, radius, color);
        case '2d':
            return createAtom2D(pos, radius, color);
        default:
            return createAtomDefault(pos, radius, color);
    }
}

function createAtomDefault(position, radius, color) {
    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshPhongMaterial({
        color: color,
        shininess: 30
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(position);
    return mesh;
}

function createAtomCartoon(position, radius, color) {
    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshToonMaterial({
        color: color
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(position);
    return mesh;
}

function createAtomGlossy(position, radius, color) {
    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshPhongMaterial({
        color: color,
        shininess: 100,
        specular: 0x444444
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(position);
    return mesh;
}

function createAtomMetallic(position, radius, color) {
    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshStandardMaterial({
        color: color,
        metalness: 0.8,
        roughness: 0.2
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(position);
    return mesh;
}

function createAtomNeon(position, radius, color) {
    const group = new THREE.Group();

    // Core sphere
    const coreGeometry = new THREE.SphereGeometry(radius * 0.8, 32, 32);
    const coreMaterial = new THREE.MeshBasicMaterial({
        color: color
    });
    const core = new THREE.Mesh(coreGeometry, coreMaterial);
    group.add(core);

    // Glow effect
    const glowGeometry = new THREE.SphereGeometry(radius, 32, 32);
    const glowMaterial = new THREE.MeshBasicMaterial({
        color: color,
        transparent: true,
        opacity: 0.3
    });
    const glow = new THREE.Mesh(glowGeometry, glowMaterial);
    group.add(glow);

    group.position.copy(position);
    return group;
}

function createAtomBubble(position, radius, color) {
    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshPhongMaterial({
        color: color,
        transparent: true,
        opacity: 0.6,
        shininess: 100,
        specular: 0xffffff
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(position);
    return mesh;
}

function createAtom2D(position, radius, color, symbol) {
    const group = new THREE.Group();

    // Circle
    const circleGeometry = new THREE.CircleGeometry(radius, 32);
    const circleMaterial = new THREE.MeshBasicMaterial({
        color: color,
        side: THREE.DoubleSide
    });
    const circle = new THREE.Mesh(circleGeometry, circleMaterial);
    group.add(circle);

    // Border
    const edgeGeometry = new THREE.RingGeometry(radius * 0.9, radius, 32);
    const edgeMaterial = new THREE.MeshBasicMaterial({
        color: 0x000000,
        side: THREE.DoubleSide
    });
    const edge = new THREE.Mesh(edgeGeometry, edgeMaterial);
    edge.position.z = 0.01;
    group.add(edge);

    group.position.copy(position);
    return group;
}

/**
 * Get color for charge visualization (coolwarm colormap)
 */
function getChargeColor(charge, minCharge = -1, maxCharge = 1) {
    // Normalize charge to 0-1 range
    const range = maxCharge - minCharge;
    const normalized = range > 0 ? (charge - minCharge) / range : 0.5;

    // Coolwarm colormap: blue -> white -> red
    const negativeColor = new THREE.Color(0x3b4cc0);
    const neutralColor = new THREE.Color(0xf7f7f7);
    const positiveColor = new THREE.Color(0xb40426);

    const result = new THREE.Color();
    if (normalized < 0.5) {
        result.lerpColors(negativeColor, neutralColor, normalized * 2);
    } else {
        result.lerpColors(neutralColor, positiveColor, (normalized - 0.5) * 2);
    }

    return result.getHex();
}

/**
 * Bond rendering with different styles
 */

/**
 * Create bond mesh between two atoms
 */
function createBond(pos1, pos2, symbol1, symbol2, radius = 0.1, style = 'cartoon', options = {}) {
    const p1 = Array.isArray(pos1)
        ? new THREE.Vector3(pos1[0], pos1[1], pos1[2])
        : pos1;
    const p2 = Array.isArray(pos2)
        ? new THREE.Vector3(pos2[0], pos2[1], pos2[2])
        : pos2;

    switch (style) {
        case 'cartoon':
            return createBondCartoon(p1, p2, symbol1, symbol2, radius);
        case 'neon':
            return createBondNeon(p1, p2, symbol1, symbol2);
        case '2d':
            return createBond2D(p1, p2);
        default:
            return createBondDefault(p1, p2, symbol1, symbol2, radius);
    }
}

function createBondDefault(pos1, pos2, symbol1, symbol2, radius) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
    const direction = new THREE.Vector3().subVectors(pos2, pos1);
    const length = direction.length();
    const halfLength = length / 2;

    const color1 = getAtomInfo(symbol1).color;
    const color2 = getAtomInfo(symbol2).color;

    // First half (from atom 1 to midpoint)
    const geometry1 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material1 = new THREE.MeshPhongMaterial({ color: color1 });
    const cylinder1 = new THREE.Mesh(geometry1, material1);

    const half1Mid = new THREE.Vector3().addVectors(pos1, midpoint).multiplyScalar(0.5);
    cylinder1.position.copy(half1Mid);
    cylinder1.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    group.add(cylinder1);

    // Second half (from midpoint to atom 2)
    const geometry2 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material2 = new THREE.MeshPhongMaterial({ color: color2 });
    const cylinder2 = new THREE.Mesh(geometry2, material2);

    const half2Mid = new THREE.Vector3().addVectors(midpoint, pos2).multiplyScalar(0.5);
    cylinder2.position.copy(half2Mid);
    cylinder2.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    group.add(cylinder2);

    return group;
}

function createBondCartoon(pos1, pos2, symbol1, symbol2, radius) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
    const direction = new THREE.Vector3().subVectors(pos2, pos1);
    const length = direction.length();
    const halfLength = length / 2;

    // Black bonds for cartoon style
    const bondColor = 0x333333;

    const geometry1 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material1 = new THREE.MeshToonMaterial({ color: bondColor });
    const cylinder1 = new THREE.Mesh(geometry1, material1);

    const half1Mid = new THREE.Vector3().addVectors(pos1, midpoint).multiplyScalar(0.5);
    cylinder1.position.copy(half1Mid);
    cylinder1.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    group.add(cylinder1);

    const geometry2 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material2 = new THREE.MeshToonMaterial({ color: bondColor });
    const cylinder2 = new THREE.Mesh(geometry2, material2);

    const half2Mid = new THREE.Vector3().addVectors(midpoint, pos2).multiplyScalar(0.5);
    cylinder2.position.copy(half2Mid);
    cylinder2.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    group.add(cylinder2);

    return group;
}

function createBondNeon(pos1, pos2, symbol1, symbol2, radius) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
    const direction = new THREE.Vector3().subVectors(pos2, pos1);
    direction.length();

    const color1 = getAtomInfo(symbol1).color;
    const color2 = getAtomInfo(symbol2).color;

    // Use line for neon style
    const points = [pos1, midpoint];
    const geometry1 = new THREE.BufferGeometry().setFromPoints(points);
    const material1 = new THREE.LineBasicMaterial({
        color: color1,
        linewidth: 2
    });
    const line1 = new THREE.Line(geometry1, material1);
    group.add(line1);

    const points2 = [midpoint, pos2];
    const geometry2 = new THREE.BufferGeometry().setFromPoints(points2);
    const material2 = new THREE.LineBasicMaterial({
        color: color2,
        linewidth: 2
    });
    const line2 = new THREE.Line(geometry2, material2);
    group.add(line2);

    return group;
}

function createBond2D(pos1, pos2, symbol1, symbol2, radius) {
    const points = [pos1, pos2];
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({
        color: 0x000000,
        linewidth: 2
    });
    return new THREE.Line(geometry, material);
}

/**
 * MolecularViewer - Main viewer for molecular structures
 */

class MolecularViewer {
    /**
     * Create a molecular viewer
     * @param {string|HTMLElement} container - Container element or selector
     * @param {Object} options - Viewer options
     */
    constructor(container, options = {}) {
        this.options = {
            style: 'cartoon',
            backgroundColor: '#1f2937',
            showBond: true,
            showCell: false,
            atomSize: 1.0,
            bondThickness: 0.15,
            bondThreshold: 1.5,
            colorBy: 'Element',  // 'Element' or 'Charge'
            ...options
        };

        this.sceneManager = new SceneManager(container, {
            backgroundColor: this.options.backgroundColor,
            viewMode: this.options.viewMode,
            rotationMode: this.options.rotationMode
        });

        this.moleculeGroup = null;
        this.molecularData = null;

        // Start animation loop
        this.sceneManager.animate();
    }

    /**
     * Set molecular data
     * @param {Object} data - Molecular data object
     * @param {string[]} data.symbols - Atom symbols
     * @param {number[][]} data.positions - Atom positions [[x,y,z], ...]
     * @param {number[]} [data.charges] - Partial charges (optional)
     * @param {number[][]} [data.cell] - Unit cell vectors (optional)
     */
    setData(data) {
        this.molecularData = data;
        this._render();
    }

    /**
     * Update viewer options
     * @param {Object} options - Options to update
     */
    setOptions(options) {
        Object.assign(this.options, options);
        if (this.molecularData) {
            this._render();
        }
    }

    _render() {
        if (!this.molecularData) return;

        // Clear previous
        if (this.moleculeGroup) {
            this.sceneManager.scene.remove(this.moleculeGroup);
            this._disposeGroup(this.moleculeGroup);
        }

        this.moleculeGroup = new THREE.Group();

        const { symbols, positions, charges, cell } = this.molecularData;
        const n = positions.length;

        // Calculate charge range for coloring
        let minCharge = 0, maxCharge = 0;
        if (charges && this.options.colorBy === 'Charge') {
            minCharge = Math.min(...charges);
            maxCharge = Math.max(...charges);
            // Ensure symmetric range around 0
            const absMax = Math.max(Math.abs(minCharge), Math.abs(maxCharge));
            minCharge = -absMax;
            maxCharge = absMax;
        }

        // Detect bonds
        const bonds = this.options.showBond
            ? detectBonds(positions, symbols, this.options.bondThreshold)
            : [];

        // Render bonds
        bonds.forEach(([i, j]) => {
            const bond = createBond(
                positions[i],
                positions[j],
                symbols[i],
                symbols[j],
                this.options.bondThickness,
                this.options.style
            );
            if (bond) {
                this.moleculeGroup.add(bond);
            }
        });

        // Render atoms
        for (let i = 0; i < n; i++) {
            let atomColor;
            if (this.options.colorBy === 'Charge' && charges) {
                atomColor = getChargeColor(charges[i], minCharge, maxCharge);
            }

            const atom = createAtom(
                positions[i],
                symbols[i],
                this.options.atomSize,
                this.options.style,
                { color: atomColor }
            );
            if (atom) {
                this.moleculeGroup.add(atom);
            }
        }

        // Render unit cell
        if (this.options.showCell && cell) {
            this._renderCell(cell);
        }

        this.sceneManager.scene.add(this.moleculeGroup);
        this.sceneManager.fitToMolecule(this.moleculeGroup);
    }

    _renderCell(cell) {
        const v0 = new THREE.Vector3(0, 0, 0);
        const v1 = new THREE.Vector3(...cell[0]);
        const v2 = new THREE.Vector3(...cell[1]);
        const v3 = new THREE.Vector3(...cell[2]);

        const points = [
            v0, v1,
            v0, v2,
            v0, v3,
            v1, v1.clone().add(v2),
            v1, v1.clone().add(v3),
            v2, v2.clone().add(v1),
            v2, v2.clone().add(v3),
            v3, v3.clone().add(v1),
            v3, v3.clone().add(v2),
            v1.clone().add(v2), v1.clone().add(v2).add(v3),
            v1.clone().add(v3), v1.clone().add(v2).add(v3),
            v2.clone().add(v3), v1.clone().add(v2).add(v3)
        ];

        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({ color: 0xff0000 });
        const cellMesh = new THREE.LineSegments(geometry, material);
        this.moleculeGroup.add(cellMesh);
    }

    _disposeGroup(group) {
        group.traverse(child => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) {
                if (Array.isArray(child.material)) {
                    child.material.forEach(m => m.dispose());
                } else {
                    child.material.dispose();
                }
            }
        });
    }

    /**
     * Take a screenshot
     * @returns {string} Data URL of the image
     */
    screenshot() {
        this.sceneManager.render();
        return this.sceneManager.renderer.domElement.toDataURL('image/png');
    }

    /**
     * Dispose the viewer
     */
    dispose() {
        if (this.moleculeGroup) {
            this.sceneManager.scene.remove(this.moleculeGroup);
            this._disposeGroup(this.moleculeGroup);
        }
        this.sceneManager.dispose();
    }
}

/**
 * ASEView.js - JavaScript library for molecular visualization
 *
 * Usage:
 * <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
 * <script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
 * <script src="aseview.umd.js"></script>
 * <script>
 *   const viewer = new ASEView.MolecularViewer('#container', {
 *     style: 'cartoon',
 *     colorBy: 'Element'
 *   });
 *   viewer.setData({
 *     symbols: ['O', 'H', 'H'],
 *     positions: [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]
 *   });
 * </script>
 */


// Version
const version = '0.1.0';

export { MolecularViewer, SceneManager, atomInfo, createAtom, createBond, detectBonds, getAtomInfo, getChargeColor, version };
