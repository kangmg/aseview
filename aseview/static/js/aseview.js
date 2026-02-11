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

    // Export to global namespace
    global.ASEView = global.ASEView || {};
    global.ASEView.MolecularViewer = MolecularViewer;
    global.ASEView.version = '0.2.0';

})(typeof window !== 'undefined' ? window : this);
