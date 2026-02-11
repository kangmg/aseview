/**
 * Atom rendering with different styles
 */
import { getAtomInfo } from '../utils/atomInfo.js';

/**
 * Create atom mesh with specified style
 */
export function createAtom(position, symbol, scale = 1.0, style = 'cartoon', options = {}) {
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
            return createAtom2D(pos, radius, color, symbol);
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
 * Positive and negative charges are normalized independently
 * 0 is always white, positive → red, negative → blue
 */
export function getChargeColor(charge, minNegative = -1, maxPositive = 1) {
    // Coolwarm colormap: blue -> white -> red
    const negativeColor = new THREE.Color(0x3b4cc0);  // blue for negative
    const neutralColor = new THREE.Color(0xf7f7f7);   // white for zero
    const positiveColor = new THREE.Color(0xb40426);  // red for positive

    const result = new THREE.Color();

    if (charge >= 0) {
        // Positive: white → red
        const normalized = maxPositive > 0 ? Math.min(charge / maxPositive, 1) : 0;
        result.lerpColors(neutralColor, positiveColor, normalized);
    } else {
        // Negative: blue → white
        const normalized = minNegative < 0 ? Math.min(Math.abs(charge) / Math.abs(minNegative), 1) : 0;
        result.lerpColors(neutralColor, negativeColor, normalized);
    }

    return result.getHex();
}
