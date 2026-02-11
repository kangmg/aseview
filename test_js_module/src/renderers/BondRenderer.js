/**
 * Bond rendering with different styles
 */
import { getAtomInfo } from '../utils/atomInfo.js';

/**
 * Create bond mesh between two atoms
 * @param {Object} options - Options including color1, color2 for custom colors
 */
export function createBond(pos1, pos2, symbol1, symbol2, radius = 0.1, style = 'cartoon', options = {}) {
    const p1 = Array.isArray(pos1)
        ? new THREE.Vector3(pos1[0], pos1[1], pos1[2])
        : pos1;
    const p2 = Array.isArray(pos2)
        ? new THREE.Vector3(pos2[0], pos2[1], pos2[2])
        : pos2;

    // Use custom colors if provided, otherwise use element colors
    const color1 = options.color1 !== undefined ? options.color1 : getAtomInfo(symbol1).color;
    const color2 = options.color2 !== undefined ? options.color2 : getAtomInfo(symbol2).color;
    const useChargeColor = options.useChargeColor || false;

    switch (style) {
        case 'cartoon':
            return createBondCartoon(p1, p2, color1, color2, radius, useChargeColor);
        case 'neon':
            return createBondNeon(p1, p2, color1, color2, radius, useChargeColor);
        case '2d':
            return createBond2D(p1, p2, color1, color2, radius);
        case 'glossy':
        case 'metallic':
        case 'bubble':
        default:
            return createBondDefault(p1, p2, color1, color2, radius, useChargeColor);
    }
}

function createBondDefault(pos1, pos2, color1, color2, radius, useChargeColor = false) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
    const direction = new THREE.Vector3().subVectors(pos2, pos1);
    const length = direction.length();
    const halfLength = length / 2;

    // Use darkened colors for charge visualization, element colors otherwise
    const bondColor1 = useChargeColor ? darkenColor(color1, 0.7) : color1;
    const bondColor2 = useChargeColor ? darkenColor(color2, 0.7) : color2;

    // First half (from atom 1 to midpoint)
    const geometry1 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material1 = new THREE.MeshPhongMaterial({ color: bondColor1 });
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
    const material2 = new THREE.MeshPhongMaterial({ color: bondColor2 });
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

function createBondCartoon(pos1, pos2, color1, color2, radius, useChargeColor = false) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);
    const direction = new THREE.Vector3().subVectors(pos2, pos1);
    const length = direction.length();
    const halfLength = length / 2;

    // Use dark gray for element colors, charge colors (darkened) for charge visualization
    const bondColor1 = useChargeColor ? darkenColor(color1, 0.6) : 0x333333;
    const bondColor2 = useChargeColor ? darkenColor(color2, 0.6) : 0x333333;

    const geometry1 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material1 = new THREE.MeshToonMaterial({ color: bondColor1 });
    const cylinder1 = new THREE.Mesh(geometry1, material1);

    const half1Mid = new THREE.Vector3().addVectors(pos1, midpoint).multiplyScalar(0.5);
    cylinder1.position.copy(half1Mid);
    cylinder1.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    group.add(cylinder1);

    const geometry2 = new THREE.CylinderGeometry(radius, radius, halfLength, 16);
    const material2 = new THREE.MeshToonMaterial({ color: bondColor2 });
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

// Helper to darken a color
function darkenColor(color, factor) {
    const c = new THREE.Color(color);
    c.multiplyScalar(factor);
    return c.getHex();
}

function createBondNeon(pos1, pos2, color1, color2, radius, useChargeColor = false) {
    const group = new THREE.Group();

    const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5);

    // For neon style, use colors directly (glow effect)
    const bondColor1 = useChargeColor ? color1 : color1;
    const bondColor2 = useChargeColor ? color2 : color2;

    // Use line for neon style
    const points = [pos1, midpoint];
    const geometry1 = new THREE.BufferGeometry().setFromPoints(points);
    const material1 = new THREE.LineBasicMaterial({
        color: bondColor1,
        linewidth: 2
    });
    const line1 = new THREE.Line(geometry1, material1);
    group.add(line1);

    const points2 = [midpoint, pos2];
    const geometry2 = new THREE.BufferGeometry().setFromPoints(points2);
    const material2 = new THREE.LineBasicMaterial({
        color: bondColor2,
        linewidth: 2
    });
    const line2 = new THREE.Line(geometry2, material2);
    group.add(line2);

    return group;
}

function createBond2D(pos1, pos2, color1, color2, radius) {
    const points = [pos1, pos2];
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({
        color: 0x000000,
        linewidth: 2
    });
    return new THREE.Line(geometry, material);
}
