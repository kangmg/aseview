/**
 * MolecularViewer - Main viewer for molecular structures
 */
import { SceneManager } from '../core/Scene.js';
import { createAtom, getChargeColor } from '../renderers/AtomRenderer.js';
import { createBond } from '../renderers/BondRenderer.js';
import { detectBonds, getAtomInfo } from '../utils/atomInfo.js';

export class MolecularViewer {
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
            bondThreshold: 1.2,
            colorBy: 'Element',  // 'Element' or 'Charge'
            normalizeCharge: true,  // true: use data range, false: use fixed range (-1 to 1)
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
        // minNegative: most negative charge (for blue), maxPositive: most positive charge (for red)
        let minNegative = -1, maxPositive = 1;
        const useChargeColor = charges && this.options.colorBy === 'Charge';

        if (useChargeColor) {
            if (this.options.normalizeCharge) {
                // Normalize positive and negative charges independently
                const positiveCharges = charges.filter(c => c > 0);
                const negativeCharges = charges.filter(c => c < 0);
                maxPositive = positiveCharges.length > 0 ? Math.max(...positiveCharges) : 0;
                minNegative = negativeCharges.length > 0 ? Math.min(...negativeCharges) : 0;
            }
            // else use fixed range (-1 to 1)
        }

        // Pre-calculate charge colors for all atoms
        const atomColors = [];
        if (useChargeColor) {
            for (let i = 0; i < n; i++) {
                atomColors[i] = getChargeColor(charges[i], minNegative, maxPositive);
            }
        }

        // Detect bonds
        const bonds = this.options.showBond
            ? detectBonds(positions, symbols, this.options.bondThreshold)
            : [];

        // Render bonds
        bonds.forEach(([i, j]) => {
            const bondOptions = {
                useChargeColor: useChargeColor
            };
            if (useChargeColor) {
                bondOptions.color1 = atomColors[i];
                bondOptions.color2 = atomColors[j];
            }

            const bond = createBond(
                positions[i],
                positions[j],
                symbols[i],
                symbols[j],
                this.options.bondThickness,
                this.options.style,
                bondOptions
            );
            if (bond) {
                this.moleculeGroup.add(bond);
            }
        });

        // Render atoms
        for (let i = 0; i < n; i++) {
            const atomColor = useChargeColor ? atomColors[i] : undefined;

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
