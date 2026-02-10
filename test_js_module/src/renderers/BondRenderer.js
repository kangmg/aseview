/**
 * Bond rendering with different styles
 */
import { getAtomInfo } from '../utils/atomInfo.js';

/**
 * Create bond mesh between two atoms
 */
export function createBond(pos1, pos2, symbol1, symbol2, radius = 0.1, style = 'cartoon', options = {}) {
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
            return createBondNeon(p1, p2, symbol1, symbol2, radius);
        case '2d':
            return createBond2D(p1, p2, symbol1, symbol2, radius);
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
    const length = direction.length();

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
