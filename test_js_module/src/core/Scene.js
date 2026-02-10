/**
 * Scene manager for Three.js rendering
 */
export class SceneManager {
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
                -1000,
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
