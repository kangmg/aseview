/**
 * ASEView.js - High-level JavaScript API for molecular visualization
 *
 * Usage:
 *   const viewer = new ASEView.MolecularViewer('#container', { style: 'cartoon' });
 *   viewer.setData({ symbols: ['O','H','H'], positions: [[0,0,0],[0.96,0,0],[-0.24,0.93,0]] });
 */

(function(global) {
    'use strict';

    // Check dependencies
    if (typeof THREE === 'undefined') {
        console.error('ASEView requires Three.js. Please load it before ASEView.');
        return;
    }

    /**
     * MolecularViewer - Main viewer class for molecular structures
     */
    class MolecularViewer {
        /**
         * Create a molecular viewer
         * @param {string|HTMLElement} container - Container element or CSS selector
         * @param {Object} options - Viewer options
         */
        constructor(container, options = {}) {
            // Get container element
            if (typeof container === 'string') {
                this.container = document.querySelector(container);
            } else {
                this.container = container;
            }

            if (!this.container) {
                throw new Error('ASEView: Container element not found');
            }

            // Default options
            this.options = {
                style: 'cartoon',
                backgroundColor: '#1f2937',
                atomSize: 0.4,
                bondThickness: 0.1,
                bondThreshold: 1.2,
                showBond: true,
                showCell: false,
                rotationMode: 'trackball', // 'trackball' or 'orbit'
                ...options
            };

            this.moleculeGroup = null;
            this.cellGroup = null;
            this.data = null;

            this._initScene();
            this._initControls();
            this._initLighting();
            this._animate();
        }

        _initScene() {
            const width = this.container.clientWidth;
            const height = this.container.clientHeight;

            // Scene
            this.scene = new THREE.Scene();
            this.scene.background = new THREE.Color(this.options.backgroundColor);

            // Camera
            this.camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
            this.camera.position.set(0, 0, 10);

            // Renderer
            this.renderer = new THREE.WebGLRenderer({ antialias: true });
            this.renderer.setSize(width, height);
            this.renderer.setPixelRatio(window.devicePixelRatio);
            this.container.appendChild(this.renderer.domElement);

            // Handle resize
            this._resizeHandler = () => this._onResize();
            window.addEventListener('resize', this._resizeHandler);
        }

        _initControls() {
            if (this.options.rotationMode === 'orbit' && THREE.OrbitControls) {
                this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
                this.controls.enableDamping = true;
                this.controls.dampingFactor = 0.05;
            } else if (THREE.TrackballControls) {
                this.controls = new THREE.TrackballControls(this.camera, this.renderer.domElement);
                this.controls.rotateSpeed = 2.0;
                this.controls.zoomSpeed = 1.2;
                this.controls.panSpeed = 0.8;
            }
        }

        _initLighting() {
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
            this.scene.add(ambientLight);

            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(5, 5, 5);
            this.scene.add(directionalLight);

            const backLight = new THREE.DirectionalLight(0xffffff, 0.3);
            backLight.position.set(-5, -5, -5);
            this.scene.add(backLight);
        }

        _animate() {
            this._animationId = requestAnimationFrame(() => this._animate());
            if (this.controls) this.controls.update();
            this.renderer.render(this.scene, this.camera);
        }

        _onResize() {
            const width = this.container.clientWidth;
            const height = this.container.clientHeight;
            this.camera.aspect = width / height;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(width, height);
        }

        /**
         * Set molecular data
         * @param {Object} data - Molecular data
         * @param {string[]} data.symbols - Atom symbols ['C', 'H', 'O', ...]
         * @param {number[][]} data.positions - Atom positions [[x,y,z], ...]
         * @param {number[][]} [data.cell] - Unit cell vectors (optional)
         */
        setData(data) {
            this.data = data;
            this._render();
        }

        /**
         * Update viewer options
         * @param {Object} options - Options to update
         */
        setOptions(options) {
            Object.assign(this.options, options);

            // Update background if changed
            if (options.backgroundColor) {
                this.scene.background = new THREE.Color(options.backgroundColor);
            }

            // Re-render if we have data
            if (this.data) {
                this._render();
            }
        }

        /**
         * Set rendering style
         * @param {string} style - Style name: 'cartoon', 'glossy', 'metallic', 'neon', 'bubble', 'rowan', '2d', 'grey'
         */
        setStyle(style) {
            this.setOptions({ style });
        }

        _render() {
            if (!this.data) return;

            // Clear previous molecule
            if (this.moleculeGroup) {
                this.scene.remove(this.moleculeGroup);
                this._disposeGroup(this.moleculeGroup);
            }
            if (this.cellGroup) {
                this.scene.remove(this.cellGroup);
                this._disposeGroup(this.cellGroup);
            }

            this.moleculeGroup = new THREE.Group();

            const { symbols, positions, cell } = this.data;
            const style = this.options.style;
            const atomScale = this.options.atomSize;
            const bondThickness = this.options.bondThickness;

            // Create atoms
            for (let i = 0; i < positions.length; i++) {
                const pos = new THREE.Vector3(...positions[i]);
                const symbol = symbols[i];
                const atom = this._createAtom(pos, symbol, atomScale, style);
                if (atom) this.moleculeGroup.add(atom);
            }

            // Create bonds
            if (this.options.showBond) {
                const bonds = this._detectBonds(positions, symbols);
                for (const [i, j] of bonds) {
                    const p1 = new THREE.Vector3(...positions[i]);
                    const p2 = new THREE.Vector3(...positions[j]);
                    const bond = this._createBond(p1, p2, symbols[i], symbols[j], bondThickness, atomScale, style);
                    if (bond) this.moleculeGroup.add(bond);
                }
            }

            this.scene.add(this.moleculeGroup);

            // Create unit cell
            if (this.options.showCell && cell) {
                this._renderCell(cell);
            }

            // Fit camera to molecule
            this._fitCamera();
        }

        _createAtom(pos, symbol, atomScale, style) {
            // Use styles.js functions if available
            switch (style) {
                case 'cartoon':
                    if (typeof createAtomStyleCartoon !== 'undefined') {
                        return createAtomStyleCartoon(pos, symbol, atomScale);
                    }
                    break;
                case 'glossy':
                    if (typeof createAtomStyleGlossy !== 'undefined') {
                        return createAtomStyleGlossy(pos, symbol, atomScale);
                    }
                    break;
                case 'metallic':
                    if (typeof createAtomStyleMetallic !== 'undefined') {
                        return createAtomStyleMetallic(pos, symbol, atomScale);
                    }
                    break;
                case 'neon':
                    if (typeof createAtomStyleNeon !== 'undefined') {
                        return createAtomStyleNeon(pos, symbol, atomScale);
                    }
                    break;
                case 'bubble':
                    if (typeof createAtomStyleBubble !== 'undefined') {
                        return createAtomStyleBubble(pos, symbol, atomScale);
                    }
                    break;
                case 'rowan':
                    if (typeof createAtomStyleRowan !== 'undefined') {
                        return createAtomStyleRowan(pos, symbol, atomScale);
                    }
                    break;
                case '2d':
                    if (typeof createAtomStyle2D !== 'undefined') {
                        return createAtomStyle2D(pos, symbol, atomScale);
                    }
                    break;
                case 'grey':
                    if (typeof createAtomStyleGrey !== 'undefined') {
                        const info = atomInfo[symbol] || atomInfo.default;
                        return createAtomStyleGrey(pos, symbol, atomScale, info.color);
                    }
                    break;
            }

            // Fallback: default style
            if (typeof createAtomStyleDefault !== 'undefined') {
                return createAtomStyleDefault(pos, symbol, atomScale);
            }

            // Ultimate fallback: simple sphere
            const info = (typeof atomInfo !== 'undefined') ? (atomInfo[symbol] || atomInfo.default) : { radius: 0.5, color: 0x808080 };
            const geometry = new THREE.SphereGeometry(info.radius * atomScale, 32, 32);
            const material = new THREE.MeshLambertMaterial({ color: info.color });
            const mesh = new THREE.Mesh(geometry, material);
            mesh.position.copy(pos);
            return mesh;
        }

        _createBond(p1, p2, sym1, sym2, bondThickness, atomScale, style) {
            // Use styles.js functions if available
            switch (style) {
                case 'cartoon':
                    if (typeof createBondStyleCartoon !== 'undefined') {
                        return createBondStyleCartoon(p1, p2, sym1, sym2, bondThickness, atomScale);
                    }
                    break;
                case 'neon':
                    if (typeof createBondStyleNeon !== 'undefined') {
                        return createBondStyleNeon(p1, p2, sym1, sym2, bondThickness, atomScale);
                    }
                    break;
                case 'rowan':
                    if (typeof createBondStyleRowan !== 'undefined') {
                        return createBondStyleRowan(p1, p2, sym1, sym2, bondThickness, atomScale);
                    }
                    break;
                case '2d':
                    if (typeof createBondStyle2D !== 'undefined') {
                        return createBondStyle2D(p1, p2, sym1, sym2, bondThickness, atomScale);
                    }
                    break;
                case 'bubble':
                    if (typeof createBondStyleBubble !== 'undefined') {
                        return createBondStyleBubble(p1, p2, sym1, sym2, bondThickness, atomScale);
                    }
                    break;
            }

            // Default bond style
            if (typeof createBondStyleDefault !== 'undefined') {
                return createBondStyleDefault(p1, p2, sym1, sym2, bondThickness, atomScale, style);
            }

            // Ultimate fallback: simple cylinder
            const mid = p1.clone().add(p2).multiplyScalar(0.5);
            const dir = p2.clone().sub(p1);
            const length = dir.length();
            const geometry = new THREE.CylinderGeometry(bondThickness, bondThickness, length, 8);
            const material = new THREE.MeshLambertMaterial({ color: 0x808080 });
            const mesh = new THREE.Mesh(geometry, material);
            mesh.position.copy(mid);
            mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.normalize());
            return mesh;
        }

        _detectBonds(positions, symbols) {
            const bonds = [];
            const threshold = this.options.bondThreshold;

            for (let i = 0; i < positions.length; i++) {
                for (let j = i + 1; j < positions.length; j++) {
                    const p1 = new THREE.Vector3(...positions[i]);
                    const p2 = new THREE.Vector3(...positions[j]);

                    const info1 = (typeof atomInfo !== 'undefined') ? (atomInfo[symbols[i]] || atomInfo.default) : { radius: 0.5 };
                    const info2 = (typeof atomInfo !== 'undefined') ? (atomInfo[symbols[j]] || atomInfo.default) : { radius: 0.5 };

                    const maxDist = (info1.radius + info2.radius) * threshold;
                    if (p1.distanceTo(p2) < maxDist) {
                        bonds.push([i, j]);
                    }
                }
            }
            return bonds;
        }

        _renderCell(cell) {
            this.cellGroup = new THREE.Group();

            const v0 = new THREE.Vector3(0, 0, 0);
            const v1 = new THREE.Vector3(...cell[0]);
            const v2 = new THREE.Vector3(...cell[1]);
            const v3 = new THREE.Vector3(...cell[2]);

            const points = [
                v0, v1, v0, v2, v0, v3,
                v1, v1.clone().add(v2), v1, v1.clone().add(v3),
                v2, v2.clone().add(v1), v2, v2.clone().add(v3),
                v3, v3.clone().add(v1), v3, v3.clone().add(v2),
                v1.clone().add(v2), v1.clone().add(v2).add(v3),
                v1.clone().add(v3), v1.clone().add(v2).add(v3),
                v2.clone().add(v3), v1.clone().add(v2).add(v3)
            ];

            const geometry = new THREE.BufferGeometry().setFromPoints(points);
            const material = new THREE.LineBasicMaterial({ color: 0x808080 });
            const lines = new THREE.LineSegments(geometry, material);
            this.cellGroup.add(lines);
            this.scene.add(this.cellGroup);
        }

        _fitCamera() {
            const box = new THREE.Box3().setFromObject(this.moleculeGroup);
            const center = box.getCenter(new THREE.Vector3());
            const size = box.getSize(new THREE.Vector3()).length();

            this.camera.position.set(center.x, center.y, center.z + size * 1.5);
            if (this.controls) {
                this.controls.target.copy(center);
            }
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
         * @returns {string} Data URL of the screenshot
         */
        screenshot() {
            this.renderer.render(this.scene, this.camera);
            return this.renderer.domElement.toDataURL('image/png');
        }

        /**
         * Dispose the viewer and clean up resources
         */
        dispose() {
            // Stop animation
            if (this._animationId) {
                cancelAnimationFrame(this._animationId);
            }

            // Remove resize listener
            window.removeEventListener('resize', this._resizeHandler);

            // Dispose controls
            if (this.controls && this.controls.dispose) {
                this.controls.dispose();
            }

            // Clear scene
            if (this.moleculeGroup) {
                this.scene.remove(this.moleculeGroup);
                this._disposeGroup(this.moleculeGroup);
            }
            if (this.cellGroup) {
                this.scene.remove(this.cellGroup);
                this._disposeGroup(this.cellGroup);
            }

            // Dispose renderer
            this.renderer.dispose();
            this.container.removeChild(this.renderer.domElement);
        }
    }

    /**
     * InteractiveViewer - Full-featured viewer with UI controls
     */
    class InteractiveViewer {
        constructor(container, options = {}) {
            if (typeof container === 'string') {
                this.container = document.querySelector(container);
            } else {
                this.container = container;
            }

            if (!this.container) {
                throw new Error('ASEView: Container element not found');
            }

            this.options = {
                style: 'cartoon',
                backgroundColor: '#1f2937',
                atomSize: 0.4,
                bondThickness: 0.1,
                bondThreshold: 1.2,
                showBond: true,
                showCell: false,
                ...options
            };

            this.data = options.data || null;
            this._buildUI();
            this._initViewer();
            this._bindEvents();

            if (this.data) {
                this.viewer.setData(this.data);
            }
        }

        _buildUI() {
            this.container.innerHTML = '';
            this.container.style.cssText = 'display:flex;position:relative;overflow:hidden;background:#111827;';

            // CSS
            const style = document.createElement('style');
            style.textContent = `
                .asv-sidebar{width:260px;background:#1f2937;border-right:1px solid #374151;padding:1rem;overflow-y:auto;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;font-size:13px;color:#e5e7eb;flex-shrink:0;}
                .asv-viewer{flex:1;position:relative;min-width:0;}
                .asv-card{background:#374151;border:1px solid #4b5563;border-radius:6px;margin-bottom:12px;overflow:hidden;}
                .asv-card-header{padding:8px 12px;background:rgba(255,255,255,0.05);font-weight:600;font-size:12px;color:#9ca3af;text-transform:uppercase;letter-spacing:0.5px;}
                .asv-card-body{padding:12px;}
                .asv-row{display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;}
                .asv-row:last-child{margin-bottom:0;}
                .asv-label{color:#d1d5db;}
                .asv-select{background:#1f2937;border:1px solid #4b5563;color:#e5e7eb;padding:6px 8px;border-radius:4px;width:100%;margin-top:6px;}
                .asv-slider{width:100px;accent-color:#3b82f6;}
                .asv-toggle{position:relative;width:40px;height:20px;}
                .asv-toggle input{opacity:0;width:0;height:0;}
                .asv-toggle-slider{position:absolute;cursor:pointer;inset:0;background:#4b5563;border-radius:20px;transition:.2s;}
                .asv-toggle-slider:before{content:"";position:absolute;height:16px;width:16px;left:2px;bottom:2px;background:#fff;border-radius:50%;transition:.2s;}
                .asv-toggle input:checked+.asv-toggle-slider{background:#3b82f6;}
                .asv-toggle input:checked+.asv-toggle-slider:before{transform:translateX(20px);}
                .asv-value{color:#9ca3af;font-size:11px;min-width:35px;text-align:right;}
            `;
            this.container.appendChild(style);

            // Sidebar
            this.sidebar = document.createElement('div');
            this.sidebar.className = 'asv-sidebar';
            this.sidebar.innerHTML = `
                <div class="asv-card">
                    <div class="asv-card-header">Style</div>
                    <div class="asv-card-body">
                        <select class="asv-select" id="asv-style">
                            <option value="cartoon">Cartoon</option>
                            <option value="glossy">Glossy</option>
                            <option value="metallic">Metallic</option>
                            <option value="default">Default</option>
                            <option value="neon">Neon</option>
                            <option value="bubble">Bubble</option>
                            <option value="rowan">Rowan</option>
                            <option value="2d">2D</option>
                            <option value="grey">Grey</option>
                        </select>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Display</div>
                    <div class="asv-card-body">
                        <div class="asv-row">
                            <span class="asv-label">Bond</span>
                            <label class="asv-toggle">
                                <input type="checkbox" id="asv-bond" checked>
                                <span class="asv-toggle-slider"></span>
                            </label>
                        </div>
                        <div class="asv-row">
                            <span class="asv-label">Cell</span>
                            <label class="asv-toggle">
                                <input type="checkbox" id="asv-cell">
                                <span class="asv-toggle-slider"></span>
                            </label>
                        </div>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Parameters</div>
                    <div class="asv-card-body">
                        <div class="asv-row">
                            <span class="asv-label">Atom Size</span>
                            <input type="range" class="asv-slider" id="asv-atomsize" min="0.1" max="1.0" step="0.05" value="0.4">
                            <span class="asv-value" id="asv-atomsize-val">0.4</span>
                        </div>
                        <div class="asv-row">
                            <span class="asv-label">Bond Width</span>
                            <input type="range" class="asv-slider" id="asv-bondwidth" min="0.02" max="0.3" step="0.02" value="0.1">
                            <span class="asv-value" id="asv-bondwidth-val">0.1</span>
                        </div>
                    </div>
                </div>
            `;
            this.container.appendChild(this.sidebar);

            // Viewer container
            this.viewerContainer = document.createElement('div');
            this.viewerContainer.className = 'asv-viewer';
            this.container.appendChild(this.viewerContainer);
        }

        _initViewer() {
            this.viewer = new MolecularViewer(this.viewerContainer, this.options);
        }

        _bindEvents() {
            // Style
            this.sidebar.querySelector('#asv-style').value = this.options.style;
            this.sidebar.querySelector('#asv-style').addEventListener('change', (e) => {
                this.viewer.setStyle(e.target.value);
            });

            // Bond toggle
            this.sidebar.querySelector('#asv-bond').checked = this.options.showBond;
            this.sidebar.querySelector('#asv-bond').addEventListener('change', (e) => {
                this.viewer.setOptions({ showBond: e.target.checked });
            });

            // Cell toggle
            this.sidebar.querySelector('#asv-cell').checked = this.options.showCell;
            this.sidebar.querySelector('#asv-cell').addEventListener('change', (e) => {
                this.viewer.setOptions({ showCell: e.target.checked });
            });

            // Atom size slider
            this.sidebar.querySelector('#asv-atomsize').value = this.options.atomSize;
            this.sidebar.querySelector('#asv-atomsize').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-atomsize-val').textContent = val.toFixed(2);
                this.viewer.setOptions({ atomSize: val });
            });

            // Bond width slider
            this.sidebar.querySelector('#asv-bondwidth').value = this.options.bondThickness;
            this.sidebar.querySelector('#asv-bondwidth').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-bondwidth-val').textContent = val.toFixed(2);
                this.viewer.setOptions({ bondThickness: val });
            });
        }

        setData(data) {
            this.data = data;
            this.viewer.setData(data);
        }

        dispose() {
            this.viewer.dispose();
            this.container.innerHTML = '';
        }
    }

    /**
     * NormalModeViewer - Viewer for molecular vibration animations
     */
    class NormalModeViewer extends MolecularViewer {
        constructor(container, options = {}) {
            super(container, options);
            this.displacements = null;
            this.amplitude = options.amplitude || 1.0;
            this.animationSpeed = options.animationSpeed || 0.05;
            this.phase = 0;
            this.isAnimating = true;
        }

        /**
         * Set molecular data with displacements
         * @param {Object} data - Molecular data with displacements
         * @param {string[]} data.symbols - Atom symbols
         * @param {number[][]} data.positions - Equilibrium positions
         * @param {number[][]} data.displacements - Displacement vectors for each atom
         */
        setData(data) {
            this.data = {
                symbols: data.symbols,
                positions: data.positions,
                cell: data.cell
            };
            this.equilibriumPositions = data.positions.map(p => [...p]);
            this.displacements = data.displacements || data.positions.map(() => [0, 0, 0]);
            this._render();
        }

        setAmplitude(amplitude) {
            this.amplitude = amplitude;
        }

        setAnimationSpeed(speed) {
            this.animationSpeed = speed;
        }

        play() { this.isAnimating = true; }
        pause() { this.isAnimating = false; }
        toggleAnimation() { this.isAnimating = !this.isAnimating; }

        _animate() {
            this._animationId = requestAnimationFrame(() => this._animate());

            if (this.isAnimating && this.displacements && this.data) {
                this.phase += this.animationSpeed;
                const factor = Math.sin(this.phase) * this.amplitude;

                // Update positions
                this.data.positions = this.equilibriumPositions.map((pos, i) => {
                    const disp = this.displacements[i];
                    return [
                        pos[0] + disp[0] * factor,
                        pos[1] + disp[1] * factor,
                        pos[2] + disp[2] * factor
                    ];
                });
                this._render();
            }

            if (this.controls) this.controls.update();
            this.renderer.render(this.scene, this.camera);
        }
    }

    /**
     * OverlayViewer - Viewer for comparing two molecular structures
     */
    class OverlayViewer extends MolecularViewer {
        constructor(container, options = {}) {
            super(container, options);
            this.structure1 = null;
            this.structure2 = null;
            this.opacity1 = options.opacity1 || 1.0;
            this.opacity2 = options.opacity2 || 0.5;
            this.moleculeGroup2 = null;
        }

        /**
         * Set two structures for overlay comparison
         * @param {Object} data1 - First molecular structure
         * @param {Object} data2 - Second molecular structure
         */
        setStructures(data1, data2) {
            this.structure1 = data1;
            this.structure2 = data2;
            this.data = data1; // Primary structure
            this._renderOverlay();
        }

        setOpacity(opacity1, opacity2) {
            this.opacity1 = opacity1;
            this.opacity2 = opacity2;
            if (this.structure1 && this.structure2) {
                this._renderOverlay();
            }
        }

        _renderOverlay() {
            // Clear previous
            if (this.moleculeGroup) {
                this.scene.remove(this.moleculeGroup);
                this._disposeGroup(this.moleculeGroup);
            }
            if (this.moleculeGroup2) {
                this.scene.remove(this.moleculeGroup2);
                this._disposeGroup(this.moleculeGroup2);
            }

            // Render first structure
            this.data = this.structure1;
            this._render();
            this._setGroupOpacity(this.moleculeGroup, this.opacity1);

            // Render second structure
            const savedGroup = this.moleculeGroup;
            this.moleculeGroup = null;
            this.data = this.structure2;
            this._render();
            this.moleculeGroup2 = this.moleculeGroup;
            this.moleculeGroup = savedGroup;
            this._setGroupOpacity(this.moleculeGroup2, this.opacity2);
        }

        _setGroupOpacity(group, opacity) {
            if (!group) return;
            group.traverse(child => {
                if (child.material) {
                    if (Array.isArray(child.material)) {
                        child.material.forEach(m => {
                            m.transparent = opacity < 1;
                            m.opacity = opacity;
                        });
                    } else {
                        child.material.transparent = opacity < 1;
                        child.material.opacity = opacity;
                    }
                }
            });
        }

        dispose() {
            if (this.moleculeGroup2) {
                this.scene.remove(this.moleculeGroup2);
                this._disposeGroup(this.moleculeGroup2);
            }
            super.dispose();
        }
    }

    /**
     * InteractiveNormalModeViewer - NormalModeViewer with UI controls
     */
    class InteractiveNormalModeViewer {
        constructor(container, options = {}) {
            if (typeof container === 'string') {
                this.container = document.querySelector(container);
            } else {
                this.container = container;
            }

            this.options = {
                style: 'cartoon',
                backgroundColor: '#1f2937',
                atomSize: 0.4,
                bondThickness: 0.1,
                amplitude: 1.0,
                animationSpeed: 0.05,
                showBond: true,
                ...options
            };

            this.data = options.data || null;
            this._buildUI();
            this._initViewer();
            this._bindEvents();

            if (this.data) {
                this.viewer.setData(this.data);
            }
        }

        _buildUI() {
            this.container.innerHTML = '';
            this.container.style.cssText = 'display:flex;position:relative;overflow:hidden;background:#111827;';

            const style = document.createElement('style');
            style.textContent = `
                .asv-sidebar{width:260px;background:#1f2937;border-right:1px solid #374151;padding:1rem;overflow-y:auto;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;font-size:13px;color:#e5e7eb;flex-shrink:0;}
                .asv-viewer{flex:1;position:relative;min-width:0;}
                .asv-card{background:#374151;border:1px solid #4b5563;border-radius:6px;margin-bottom:12px;overflow:hidden;}
                .asv-card-header{padding:8px 12px;background:rgba(255,255,255,0.05);font-weight:600;font-size:12px;color:#9ca3af;text-transform:uppercase;letter-spacing:0.5px;}
                .asv-card-body{padding:12px;}
                .asv-row{display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;}
                .asv-row:last-child{margin-bottom:0;}
                .asv-label{color:#d1d5db;}
                .asv-select{background:#1f2937;border:1px solid #4b5563;color:#e5e7eb;padding:6px 8px;border-radius:4px;width:100%;margin-top:6px;}
                .asv-slider{width:100px;accent-color:#3b82f6;}
                .asv-toggle{position:relative;width:40px;height:20px;}
                .asv-toggle input{opacity:0;width:0;height:0;}
                .asv-toggle-slider{position:absolute;cursor:pointer;inset:0;background:#4b5563;border-radius:20px;transition:.2s;}
                .asv-toggle-slider:before{content:"";position:absolute;height:16px;width:16px;left:2px;bottom:2px;background:#fff;border-radius:50%;transition:.2s;}
                .asv-toggle input:checked+.asv-toggle-slider{background:#3b82f6;}
                .asv-toggle input:checked+.asv-toggle-slider:before{transform:translateX(20px);}
                .asv-value{color:#9ca3af;font-size:11px;min-width:35px;text-align:right;}
                .asv-btn{background:#3b82f6;color:#fff;border:none;padding:8px 16px;border-radius:4px;cursor:pointer;font-size:13px;width:100%;}
                .asv-btn:hover{background:#2563eb;}
            `;
            this.container.appendChild(style);

            this.sidebar = document.createElement('div');
            this.sidebar.className = 'asv-sidebar';
            this.sidebar.innerHTML = `
                <div class="asv-card">
                    <div class="asv-card-header">Animation</div>
                    <div class="asv-card-body">
                        <button class="asv-btn" id="asv-playpause">Pause</button>
                        <div class="asv-row" style="margin-top:12px;">
                            <span class="asv-label">Amplitude</span>
                            <input type="range" class="asv-slider" id="asv-amplitude" min="0.1" max="3.0" step="0.1" value="1.0">
                            <span class="asv-value" id="asv-amplitude-val">1.0</span>
                        </div>
                        <div class="asv-row">
                            <span class="asv-label">Speed</span>
                            <input type="range" class="asv-slider" id="asv-speed" min="0.01" max="0.15" step="0.01" value="0.05">
                            <span class="asv-value" id="asv-speed-val">0.05</span>
                        </div>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Style</div>
                    <div class="asv-card-body">
                        <select class="asv-select" id="asv-style">
                            <option value="cartoon">Cartoon</option>
                            <option value="glossy">Glossy</option>
                            <option value="metallic">Metallic</option>
                            <option value="default">Default</option>
                            <option value="neon">Neon</option>
                            <option value="bubble">Bubble</option>
                        </select>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Parameters</div>
                    <div class="asv-card-body">
                        <div class="asv-row">
                            <span class="asv-label">Atom Size</span>
                            <input type="range" class="asv-slider" id="asv-atomsize" min="0.1" max="1.0" step="0.05" value="0.4">
                            <span class="asv-value" id="asv-atomsize-val">0.4</span>
                        </div>
                    </div>
                </div>
            `;
            this.container.appendChild(this.sidebar);

            this.viewerContainer = document.createElement('div');
            this.viewerContainer.className = 'asv-viewer';
            this.container.appendChild(this.viewerContainer);
        }

        _initViewer() {
            this.viewer = new NormalModeViewer(this.viewerContainer, this.options);
        }

        _bindEvents() {
            this.sidebar.querySelector('#asv-playpause').addEventListener('click', (e) => {
                this.viewer.toggleAnimation();
                e.target.textContent = this.viewer.isAnimating ? 'Pause' : 'Play';
            });

            this.sidebar.querySelector('#asv-amplitude').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-amplitude-val').textContent = val.toFixed(1);
                this.viewer.setAmplitude(val);
            });

            this.sidebar.querySelector('#asv-speed').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-speed-val').textContent = val.toFixed(2);
                this.viewer.setAnimationSpeed(val);
            });

            this.sidebar.querySelector('#asv-style').addEventListener('change', (e) => {
                this.viewer.setStyle(e.target.value);
            });

            this.sidebar.querySelector('#asv-atomsize').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-atomsize-val').textContent = val.toFixed(2);
                this.viewer.setOptions({ atomSize: val });
            });
        }

        setData(data) {
            this.viewer.setData(data);
        }

        dispose() {
            this.viewer.dispose();
            this.container.innerHTML = '';
        }
    }

    /**
     * InteractiveOverlayViewer - OverlayViewer with UI controls
     */
    class InteractiveOverlayViewer {
        constructor(container, options = {}) {
            if (typeof container === 'string') {
                this.container = document.querySelector(container);
            } else {
                this.container = container;
            }

            this.options = {
                style: 'cartoon',
                backgroundColor: '#1f2937',
                atomSize: 0.4,
                bondThickness: 0.1,
                opacity1: 1.0,
                opacity2: 0.5,
                showBond: true,
                ...options
            };

            this._buildUI();
            this._initViewer();
            this._bindEvents();
        }

        _buildUI() {
            this.container.innerHTML = '';
            this.container.style.cssText = 'display:flex;position:relative;overflow:hidden;background:#111827;';

            const style = document.createElement('style');
            style.textContent = `
                .asv-sidebar{width:260px;background:#1f2937;border-right:1px solid #374151;padding:1rem;overflow-y:auto;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;font-size:13px;color:#e5e7eb;flex-shrink:0;}
                .asv-viewer{flex:1;position:relative;min-width:0;}
                .asv-card{background:#374151;border:1px solid #4b5563;border-radius:6px;margin-bottom:12px;overflow:hidden;}
                .asv-card-header{padding:8px 12px;background:rgba(255,255,255,0.05);font-weight:600;font-size:12px;color:#9ca3af;text-transform:uppercase;letter-spacing:0.5px;}
                .asv-card-body{padding:12px;}
                .asv-row{display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;}
                .asv-row:last-child{margin-bottom:0;}
                .asv-label{color:#d1d5db;}
                .asv-select{background:#1f2937;border:1px solid #4b5563;color:#e5e7eb;padding:6px 8px;border-radius:4px;width:100%;margin-top:6px;}
                .asv-slider{width:100px;accent-color:#3b82f6;}
                .asv-value{color:#9ca3af;font-size:11px;min-width:35px;text-align:right;}
            `;
            this.container.appendChild(style);

            this.sidebar = document.createElement('div');
            this.sidebar.className = 'asv-sidebar';
            this.sidebar.innerHTML = `
                <div class="asv-card">
                    <div class="asv-card-header">Overlay</div>
                    <div class="asv-card-body">
                        <div class="asv-row">
                            <span class="asv-label">Structure 1</span>
                            <input type="range" class="asv-slider" id="asv-opacity1" min="0" max="1" step="0.1" value="1.0">
                            <span class="asv-value" id="asv-opacity1-val">1.0</span>
                        </div>
                        <div class="asv-row">
                            <span class="asv-label">Structure 2</span>
                            <input type="range" class="asv-slider" id="asv-opacity2" min="0" max="1" step="0.1" value="0.5">
                            <span class="asv-value" id="asv-opacity2-val">0.5</span>
                        </div>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Style</div>
                    <div class="asv-card-body">
                        <select class="asv-select" id="asv-style">
                            <option value="cartoon">Cartoon</option>
                            <option value="glossy">Glossy</option>
                            <option value="metallic">Metallic</option>
                            <option value="default">Default</option>
                            <option value="neon">Neon</option>
                            <option value="bubble">Bubble</option>
                        </select>
                    </div>
                </div>
                <div class="asv-card">
                    <div class="asv-card-header">Parameters</div>
                    <div class="asv-card-body">
                        <div class="asv-row">
                            <span class="asv-label">Atom Size</span>
                            <input type="range" class="asv-slider" id="asv-atomsize" min="0.1" max="1.0" step="0.05" value="0.4">
                            <span class="asv-value" id="asv-atomsize-val">0.4</span>
                        </div>
                    </div>
                </div>
            `;
            this.container.appendChild(this.sidebar);

            this.viewerContainer = document.createElement('div');
            this.viewerContainer.className = 'asv-viewer';
            this.container.appendChild(this.viewerContainer);
        }

        _initViewer() {
            this.viewer = new OverlayViewer(this.viewerContainer, this.options);
        }

        _bindEvents() {
            this.sidebar.querySelector('#asv-opacity1').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-opacity1-val').textContent = val.toFixed(1);
                this.options.opacity1 = val;
                this.viewer.setOpacity(this.options.opacity1, this.options.opacity2);
            });

            this.sidebar.querySelector('#asv-opacity2').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-opacity2-val').textContent = val.toFixed(1);
                this.options.opacity2 = val;
                this.viewer.setOpacity(this.options.opacity1, this.options.opacity2);
            });

            this.sidebar.querySelector('#asv-style').addEventListener('change', (e) => {
                this.viewer.setStyle(e.target.value);
            });

            this.sidebar.querySelector('#asv-atomsize').addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                this.sidebar.querySelector('#asv-atomsize-val').textContent = val.toFixed(2);
                this.viewer.setOptions({ atomSize: val });
            });
        }

        setStructures(data1, data2) {
            this.viewer.setStructures(data1, data2);
        }

        dispose() {
            this.viewer.dispose();
            this.container.innerHTML = '';
        }
    }

    // Export to global namespace
    global.ASEView = global.ASEView || {};
    global.ASEView.MolecularViewer = MolecularViewer;
    global.ASEView.NormalModeViewer = NormalModeViewer;
    global.ASEView.OverlayViewer = OverlayViewer;
    global.ASEView.InteractiveViewer = InteractiveViewer;
    global.ASEView.InteractiveNormalModeViewer = InteractiveNormalModeViewer;
    global.ASEView.InteractiveOverlayViewer = InteractiveOverlayViewer;
    global.ASEView.version = '0.4.0';

})(typeof window !== 'undefined' ? window : this);
