
const atomInfo = {
    'H': { radius: 0.31, color: 0xFFFFFF }, 'He': { radius: 0.28, color: 0xD9FFFF },
    'Li': { radius: 1.28, color: 0xCC80FF }, 'Be': { radius: 0.96, color: 0xC2FF00 },
    'B': { radius: 0.84, color: 0xFFB5B5 }, 'C': { radius: 0.76, color: 0x909090 },
    'N': { radius: 0.71, color: 0x3050F8 }, 'O': { radius: 0.66, color: 0xFF0D0D },
    'F': { radius: 0.57, color: 0x90E050 }, 'Ne': { radius: 0.58, color: 0xB3E3F5 },
    'Na': { radius: 1.66, color: 0xAB5CF2 }, 'Mg': { radius: 1.41, color: 0x8AFF00 },
    'Al': { radius: 1.21, color: 0xBFA6A6 }, 'Si': { radius: 1.11, color: 0xF0C8A0 },
    'P': { radius: 1.07, color: 0xFF8000 }, 'S': { radius: 1.05, color: 0xFFFF30 },
    'Cl': { radius: 1.02, color: 0x1FF01F }, 'Ar': { radius: 1.06, color: 0x80D1E3 },
    'K': { radius: 2.03, color: 0x8F40D4 }, 'Ca': { radius: 1.76, color: 0x3DFF00 },
    'Sc': { radius: 1.70, color: 0xE6E6E6 }, 'Ti': { radius: 1.60, color: 0xBFC2C7 },
    'V': { radius: 1.53, color: 0xA6A6AB }, 'Cr': { radius: 1.39, color: 0x8A99C7 },
    'Mn': { radius: 1.61, color: 0x9C7AC7 }, 'Fe': { radius: 1.52, color: 0xE06633 },
    'Co': { radius: 1.50, color: 0xF090A0 }, 'Ni': { radius: 1.24, color: 0x50D050 },
    'Cu': { radius: 1.32, color: 0xC88033 }, 'Zn': { radius: 1.22, color: 0x7D80B0 },
    'Ga': { radius: 1.22, color: 0xC28F8F }, 'Ge': { radius: 1.20, color: 0x668F8F },
    'As': { radius: 1.19, color: 0xBD80E3 }, 'Se': { radius: 1.20, color: 0xFFA100 },
    'Br': { radius: 1.20, color: 0xA62929 }, 'Kr': { radius: 1.16, color: 0x5CB8D1 },
    'Rb': { radius: 2.20, color: 0x702EB0 }, 'Sr': { radius: 1.95, color: 0x00FF00 },
    'Y': { radius: 1.90, color: 0x94FFFF }, 'Zr': { radius: 1.75, color: 0x94E0E0 },
    'Nb': { radius: 1.64, color: 0x73C2C9 }, 'Mo': { radius: 1.54, color: 0x54B5B5 },
    'Tc': { radius: 1.47, color: 0x3B9E9E }, 'Ru': { radius: 1.46, color: 0x248F8F },
    'Rh': { radius: 1.42, color: 0x0A7D8C }, 'Pd': { radius: 1.39, color: 0x006985 },
    'Ag': { radius: 1.45, color: 0xC0C0C0 }, 'Cd': { radius: 1.44, color: 0xFFD98F },
    'In': { radius: 1.42, color: 0xA67573 }, 'Sn': { radius: 1.39, color: 0x668080 },
    'Sb': { radius: 1.39, color: 0x9E63B5 }, 'Te': { radius: 1.38, color: 0xD47A00 },
    'I': { radius: 1.39, color: 0x940094 }, 'Xe': { radius: 1.40, color: 0x429EB0 },
    'Cs': { radius: 2.44, color: 0x57178F }, 'Ba': { radius: 2.15, color: 0x00C900 },
    'La': { radius: 2.07, color: 0x70D4FF }, 'Ce': { radius: 2.04, color: 0xFFFFC7 },
    'Pr': { radius: 2.03, color: 0xD9FFC7 }, 'Nd': { radius: 2.01, color: 0xC7FFC7 },
    'Pm': { radius: 1.99, color: 0xA3FFC7 }, 'Sm': { radius: 1.98, color: 0x8FFFC7 },
    'Eu': { radius: 1.98, color: 0x61FFC7 }, 'Gd': { radius: 1.96, color: 0x45FFC7 },
    'Tb': { radius: 1.94, color: 0x30FFC7 }, 'Dy': { radius: 1.92, color: 0x1FFFC7 },
    'Ho': { radius: 1.92, color: 0x00FF9C }, 'Er': { radius: 1.89, color: 0x00E675 },
    'Tm': { radius: 1.90, color: 0x00D452 }, 'Yb': { radius: 1.87, color: 0x00BF38 },
    'Lu': { radius: 1.87, color: 0x00AB24 }, 'Hf': { radius: 1.75, color: 0x4DC2FF },
    'Ta': { radius: 1.70, color: 0x4DA6FF }, 'W': { radius: 1.62, color: 0x2194D6 },
    'Re': { radius: 1.51, color: 0x267DAB }, 'Os': { radius: 1.44, color: 0x266696 },
    'Ir': { radius: 1.41, color: 0x175487 }, 'Pt': { radius: 1.36, color: 0xD0D0E0 },
    'Au': { radius: 1.36, color: 0xFFD123 }, 'Hg': { radius: 1.32, color: 0xB8B8D0 },
    'Tl': { radius: 1.45, color: 0xA6544D }, 'Pb': { radius: 1.46, color: 0x575961 },
    'Bi': { radius: 1.48, color: 0x9E4FB5 }, 'Po': { radius: 1.40, color: 0xAB5C00 },
    'At': { radius: 1.50, color: 0x754F45 }, 'Rn': { radius: 1.50, color: 0x428296 },
    'Fr': { radius: 2.60, color: 0x420066 }, 'Ra': { radius: 2.21, color: 0x007D00 },
    'Ac': { radius: 2.15, color: 0x70ABFA }, 'Th': { radius: 2.06, color: 0x00BAFF },
    'Pa': { radius: 2.00, color: 0x00A1FF }, 'U': { radius: 1.96, color: 0x008FFF },
    'Np': { radius: 1.90, color: 0x0080FF }, 'Pu': { radius: 1.87, color: 0x006BFF },
    'Am': { radius: 1.80, color: 0x545CF2 }, 'Cm': { radius: 1.69, color: 0x785CE3 },
    // Transuranium elements (not in original list)
    'Bk': { radius: 1.70, color: 0x8A4FE3 }, 'Cf': { radius: 1.70, color: 0xA136D4 },
    'Es': { radius: 1.70, color: 0xB31FD4 }, 'Fm': { radius: 1.70, color: 0xB31FBA },
    'Md': { radius: 1.70, color: 0xB30DA6 }, 'No': { radius: 1.70, color: 0xBD0D87 },
    'Lr': { radius: 1.70, color: 0xC70066 }, 'Rf': { radius: 1.70, color: 0xCC0059 },
    'Db': { radius: 1.70, color: 0xD1004F }, 'Sg': { radius: 1.70, color: 0xD90045 },
    'Bh': { radius: 1.70, color: 0xE00038 }, 'Hs': { radius: 1.70, color: 0xE6002E },
    'Mt': { radius: 1.70, color: 0xEB0026 },
    'X': { radius: 0.7, color: 0xFFC0CB },
    'default': { radius: 0.7, color: 0xFFC0CB } // Pink as default for unknown elements
};

// CPK color scheme (classic Corey-Pauling-Koltun)
const cpkColors = {
    'H': 0xFFFFFF, 'He': 0xFFC0CB, 'Li': 0xB22222, 'Be': 0xFF1493,
    'B': 0x00FF00, 'C': 0xC8C8C8, 'N': 0x8F8FFF, 'O': 0xF00000,
    'F': 0xDAA520, 'Ne': 0xFF1493, 'Na': 0x0000FF, 'Mg': 0x228B22,
    'Al': 0x808090, 'Si': 0xDAA520, 'P': 0xFFA500, 'S': 0xFFC832,
    'Cl': 0x00FF00, 'Ar': 0xFF1493, 'K': 0xFF1493, 'Ca': 0x808090,
    'Sc': 0xFF1493, 'Ti': 0x808090, 'V': 0xFF1493, 'Cr': 0x808090,
    'Mn': 0x808090, 'Fe': 0xFFA500, 'Co': 0xFF1493, 'Ni': 0xA52A2A,
    'Cu': 0xA52A2A, 'Zn': 0xA52A2A, 'Ga': 0xFF1493, 'Ge': 0xFF1493,
    'As': 0xFF1493, 'Se': 0xFF1493, 'Br': 0xA52A2A, 'Kr': 0xFF1493,
    'Rb': 0xFF1493, 'Sr': 0xFF1493, 'Y': 0xFF1493, 'Zr': 0xFF1493,
    'Nb': 0xFF1493, 'Mo': 0xFF1493, 'Tc': 0xFF1493, 'Ru': 0xFF1493,
    'Rh': 0xFF1493, 'Pd': 0xFF1493, 'Ag': 0x808090, 'Cd': 0xFF1493,
    'In': 0xFF1493, 'Sn': 0xFF1493, 'Sb': 0xFF1493, 'Te': 0xFF1493,
    'I': 0xA020F0, 'Xe': 0xFF1493, 'Cs': 0xFF1493, 'Ba': 0xFFA500,
    'La': 0xFF1493, 'Ce': 0xFF1493, 'Pr': 0xFF1493, 'Nd': 0xFF1493,
    'Pm': 0xFF1493, 'Sm': 0xFF1493, 'Eu': 0xFF1493, 'Gd': 0xFF1493,
    'Tb': 0xFF1493, 'Dy': 0xFF1493, 'Ho': 0xFF1493, 'Er': 0xFF1493,
    'Tm': 0xFF1493, 'Yb': 0xFF1493, 'Lu': 0xFF1493, 'Hf': 0xFF1493,
    'Ta': 0xFF1493, 'W': 0xFF1493, 'Re': 0xFF1493, 'Os': 0xFF1493,
    'Ir': 0xFF1493, 'Pt': 0xFF1493, 'Au': 0xDAA520, 'Hg': 0xFF1493,
    'Tl': 0xFF1493, 'Pb': 0xFF1493, 'Bi': 0xFF1493, 'Po': 0xFF1493,
    'At': 0xFF1493, 'Rn': 0xFFFFFF, 'Fr': 0xFFFFFF, 'Ra': 0xFFFFFF,
    'Ac': 0xFFFFFF, 'Th': 0xFF1493, 'Pa': 0xFFFFFF, 'U': 0xFF1493,
    'Np': 0xFFFFFF, 'Pu': 0xFFFFFF, 'Am': 0xFFFFFF, 'Cm': 0xFFFFFF,
};

// Active color scheme: 'jmol' (default) or 'cpk'
let atomColorScheme = 'jmol';

/**
 * Return the element color for `symbol` under the current atomColorScheme.
 * Falls back to jmol (atomInfo.color) for unknown elements or schemes.
 */
function getAtomColorByScheme(symbol) {
    if (atomColorScheme === 'cpk') {
        return cpkColors[symbol] !== undefined
            ? cpkColors[symbol]
            : (atomInfo[symbol] || atomInfo['default']).color;
    }
    return (atomInfo[symbol] || atomInfo['default']).color; // 'jmol'
}

function getStandardMaterial(color, type) {
    switch (type) {
        case 'glossy': return new THREE.MeshPhongMaterial({ color: color, shininess: 100, specular: 0x222222 });
        case 'metallic': return new THREE.MeshStandardMaterial({ color: color, metalness: 0.3, roughness: 0.4 });
        default: return new THREE.MeshLambertMaterial({ color: color });
    }
}

function getBondPoints(p1, p2, sym1, sym2, atomScale, style) {
    const dir = p2.clone().sub(p1).normalize();

    let offsetFactor1 = 0.75;
    let offsetFactor2 = 0.75;

    if (style === '2d') {
        offsetFactor1 = 1.0;
        offsetFactor2 = 1.0;
    } else if (style === 'neon') {
        offsetFactor1 = 0.6;
        offsetFactor2 = 0.6;
    }

    const offset1 = (atomInfo[sym1] || atomInfo.default).radius * atomScale * offsetFactor1;
    const offset2 = (atomInfo[sym2] || atomInfo.default).radius * atomScale * offsetFactor2;

    const startPos = p1.clone().add(dir.clone().multiplyScalar(offset1));
    const endPos = p2.clone().sub(dir.clone().multiplyScalar(offset2));

    if (startPos.distanceTo(endPos) <= 0.01) return null;

    return { startPos, endPos, dir };
}

function createAtomStyle2D(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    canvas.width = 128;
    canvas.height = 128;
    context.beginPath();
    context.arc(64, 64, 60, 0, 2 * Math.PI, false);
    context.fillStyle = 'black';
    context.fill();
    context.beginPath();
    context.arc(64, 64, 58, 0, 2 * Math.PI, false);
    context.fillStyle = 'white';
    context.fill();
    context.font = 'bold 50px Arial';
    context.fillStyle = 'black';
    context.textAlign = 'center';
    context.textBaseline = 'middle';
    context.fillText(symbol, 64, 64);
    const texture = new THREE.CanvasTexture(canvas);
    const spriteMaterial = new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: true, alphaTest: 0.5 });
    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.position.copy(pos);
    sprite.scale.set(scaledRadius * 2.5, scaledRadius * 2.5, 1);
    return sprite;
}

function createAtomStyleCartoon(pos, symbol, atomScale, helpers) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const lod = (helpers && helpers._lod) || { sphere: 64, sphereHigh: 128 };
    const seg = lod.sphere;

    const geometry = getCachedSphereGeometry(scaledRadius, seg);
    const toonMaterial = getCachedMaterial('toon', getAtomColorByScheme(symbol), `cartoon_${seg}`);
    const sphere = new THREE.Mesh(geometry, toonMaterial);

    const outlineGeometry = getCachedSphereGeometry(scaledRadius * 1.05, seg);
    const outlineMaterial = getCachedMaterial('outline', 0x000000);
    const outline = new THREE.Mesh(outlineGeometry, outlineMaterial);

    const toonGroup = new THREE.Group();
    toonGroup.add(outline);
    toonGroup.add(sphere);
    toonGroup.position.copy(pos);
    return toonGroup;
}

function createAtomStyleNeon(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const atomColor = new THREE.Color(getAtomColorByScheme(symbol));

    // Use single sprite for glow effect only
    const glowCanvas = document.createElement('canvas');
    glowCanvas.width = 512;
    glowCanvas.height = 512;
    const glowContext = glowCanvas.getContext('2d');

    const center = 256;
    const maxRadius = 256;

    // Natural radial gradient - gradually increases from 40%, fades to 0 at edge
    const glowGradient = glowContext.createRadialGradient(center, center, 0, center, center, maxRadius);
    const glowColorStyle = atomColor.getStyle();

    glowGradient.addColorStop(0.0, 'rgba(0,0,0,0)');                            // Center: fully transparent
    glowGradient.addColorStop(0.25, 'rgba(0,0,0,0)');                           // Keep transparent until 25%
    glowGradient.addColorStop(0.35, `${glowColorStyle.slice(0, -1)}, 0.05)`); // Start gradually
    glowGradient.addColorStop(0.45, `${glowColorStyle.slice(0, -1)}, 0.15)`);
    glowGradient.addColorStop(0.60, `${glowColorStyle.slice(0, -1)}, 0.35)`);
    glowGradient.addColorStop(0.75, `${glowColorStyle.slice(0, -1)}, 0.6)`);
    glowGradient.addColorStop(0.88, `${glowColorStyle.slice(0, -1)}, 0.85)`);
    glowGradient.addColorStop(0.95, `${glowColorStyle.slice(0, -1)}, 1.0)`);  // Max intensity near edge
    glowGradient.addColorStop(1.0, `${glowColorStyle.slice(0, -1)}, 0.0)`);   // Fade out smoothly


    glowContext.fillStyle = glowGradient;
    glowContext.fillRect(0, 0, 512, 512);

    const glowTexture = new THREE.CanvasTexture(glowCanvas);
    glowTexture.needsUpdate = true;

    const glowMaterial = new THREE.SpriteMaterial({
        map: glowTexture,
        blending: THREE.AdditiveBlending,
        transparent: true,
        depthTest: false,
        alphaTest: 0.001
    });

    const glowSprite = new THREE.Sprite(glowMaterial);
    glowSprite.position.copy(pos);
    // Scale to match atom edge precisely
    glowSprite.scale.set(scaledRadius * 2.2, scaledRadius * 2.2, 1);

    // Add metadata
    glowSprite.userData = {
        symbol: symbol,
        radius: scaledRadius,
        originalColor: atomColor.clone()
    };

    return glowSprite;
}

function createAtomStyleGlossy(pos, symbol, atomScale, helpers) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const lod = (helpers && helpers._lod) || { sphere: 64 };
    const geometry = getCachedSphereGeometry(scaledRadius, lod.sphere);
    const material = getCachedMaterial('glossy', getAtomColorByScheme(symbol));
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleMetallic(pos, symbol, atomScale, helpers) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const lod = (helpers && helpers._lod) || { sphere: 64 };
    const geometry = getCachedSphereGeometry(scaledRadius, lod.sphere);
    const material = getCachedMaterial('metallic', getAtomColorByScheme(symbol));
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleDefault(pos, symbol, atomScale, helpers) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const lod = (helpers && helpers._lod) || { sphere: 64 };
    const geometry = getCachedSphereGeometry(scaledRadius, lod.sphere);
    const material = getCachedMaterial('lambert', getAtomColorByScheme(symbol));
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleRowan(pos, symbol, atomScale, helpers) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const lod = (helpers && helpers._lod) || { sphere: 64, sphereHigh: 128 };
    const group = new THREE.Group();

    const geometry = getCachedSphereGeometry(scaledRadius, lod.sphereHigh);
    const material = getCachedMaterial('lambertFlat', getAtomColorByScheme(symbol));
    const sphere = new THREE.Mesh(geometry, material);
    group.add(sphere);

    const outlineGeometry = getCachedSphereGeometry(scaledRadius * 1.05, lod.sphere);
    const outlineMaterial = getCachedMaterial('outline', 0x000000);
    const outline = new THREE.Mesh(outlineGeometry, outlineMaterial);
    group.add(outline);

    group.position.copy(pos);
    return group;
}

function createAtomStyleBubble(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;

    // Create canvas with radial gradient for bubble effect
    const canvas = document.createElement('canvas');
    canvas.width = 512;
    canvas.height = 512;
    const context = canvas.getContext('2d');

    const centerX = 256;
    const centerY = 256;
    const radius = 240;

    // Create radial gradient for bubble effect
    const gradient = context.createRadialGradient(
        centerX - radius * 0.3, centerY - radius * 0.3, 0,  // Highlight position
        centerX, centerY, radius
    );

    const atomColor = new THREE.Color(getAtomColorByScheme(symbol));
    const colorStyle = atomColor.getStyle();

    // Bubble gradient: bright highlight -> main color -> darker edge
    gradient.addColorStop(0, 'rgba(255, 255, 255, 0.9)');  // Bright highlight
    gradient.addColorStop(0.3, `${colorStyle.slice(0, -1)}, 0.7)`);  // Main color
    gradient.addColorStop(0.7, `${colorStyle.slice(0, -1)}, 0.5)`);  // Slightly darker
    gradient.addColorStop(1.0, `${colorStyle.slice(0, -1)}, 0.3)`);  // Transparent edge

    context.fillStyle = gradient;
    context.beginPath();
    context.arc(centerX, centerY, radius, 0, 2 * Math.PI);
    context.fill();

    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;

    const spriteMaterial = new THREE.SpriteMaterial({
        map: texture,
        transparent: true,
        opacity: 0.8,
        depthTest: true,
        depthWrite: false
    });

    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.position.copy(pos);
    // Match size with other styles (2x radius)
    sprite.scale.set(scaledRadius * 2, scaledRadius * 2, 1);

    return sprite;
}

function createAtomStyleGrey(pos, symbol, atomScale, color) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;

    // Remove 3D sphere, use canvas-drawn sprite only
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    canvas.width = 256;
    canvas.height = 256;

    const centerX = 128;
    const centerY = 128;
    const circleRadius = 120;

    // Draw black border
    context.beginPath();
    context.arc(centerX, centerY, circleRadius, 0, 2 * Math.PI, false);
    context.fillStyle = 'black';
    context.fill();

    // Draw circular background (sphere role) - slightly smaller than border
    context.beginPath();
    context.arc(centerX, centerY, circleRadius - 3, 0, 2 * Math.PI, false);
    context.fillStyle = `#${color.toString(16).padStart(6, '0')}`; // Convert color to hex
    context.fill();

    // Highlight effect with slightly brighter color
    context.beginPath();
    context.arc(centerX, centerY, circleRadius - 7, 0, 2 * Math.PI, false);
    const lighterColor = new THREE.Color(color).multiplyScalar(1.1);
    context.fillStyle = `#${lighterColor.getHexString()}`;
    context.fill();

    // Draw text
    context.font = 'bold 80px Arial';
    context.fillStyle = 'white'; // White text
    context.strokeStyle = 'black';
    context.lineWidth = 2;
    context.textAlign = 'center';
    context.textBaseline = 'middle';

    // Text outline effect
    context.strokeText(symbol, centerX, centerY);
    context.fillText(symbol, centerX, centerY);

    const texture = new THREE.CanvasTexture(canvas);
    const spriteMaterial = new THREE.SpriteMaterial({
        map: texture,
        transparent: true,
        depthTest: true,     // Enable depth test
        depthWrite: false,   // Disable depth write for transparent object
        alphaTest: 0.5
    });

    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.position.copy(pos); // Place directly at given position
    sprite.scale.set(scaledRadius * 2.5, scaledRadius * 2.5, 1);
    sprite.center.set(0.5, 0.5);

    return sprite; // Return sprite directly, not Group
}

function createBondStyle2D(p1, p2, sym1, sym2, bondThickness, atomScale) {
    const bondPoints = getBondPoints(p1, p2, sym1, sym2, atomScale, '2d');
    if (!bondPoints) return null;
    const { startPos, endPos, dir } = bondPoints;

    const bondGroup = new THREE.Group();
    const distance = startPos.distanceTo(endPos);

    const geometry = new THREE.CylinderGeometry(bondThickness, bondThickness, distance, 8, 1);
    const material = new THREE.MeshBasicMaterial({ color: 0xFFFFFF });
    const bond = new THREE.Mesh(geometry, material);
    bond.position.copy(startPos).add(endPos).multiplyScalar(0.5);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
    bondGroup.add(bond);

    const outline_geo = new THREE.CylinderGeometry(bondThickness + 0.01, bondThickness + 0.01, distance, 8, 1);
    const outline_mat = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.BackSide });
    const bond_outline = new THREE.Mesh(outline_geo, outline_mat);
    bond_outline.position.copy(startPos).add(endPos).multiplyScalar(0.5);
    bond_outline.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
    bondGroup.add(bond_outline);

    return bondGroup;
}

function createBondStyleCartoon(p1, p2, sym1, sym2, bondThickness, atomScale) {
    const bondPoints = getBondPoints(p1, p2, sym1, sym2, atomScale, 'cartoon');
    if (!bondPoints) return null;
    const { startPos, endPos, dir } = bondPoints;
    const distance = startPos.distanceTo(endPos);

    // Create gradient map for toon shading
    const colors = new Uint8Array(3);
    for (let c = 0; c <= 2; c++) {
        colors[c] = (c / 2) * 256;
    }
    const gradientMap = new THREE.DataTexture(colors, colors.length, 1, THREE.LuminanceFormat);
    gradientMap.needsUpdate = true;

    const geometry = new THREE.CylinderGeometry(bondThickness, bondThickness, distance, 8);
    const material = new THREE.MeshToonMaterial({
        color: 0x000000,
        gradientMap: gradientMap
    });
    const bond = new THREE.Mesh(geometry, material);
    bond.position.lerpVectors(startPos, endPos, 0.5);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
    return bond;
}

function createBondStyleDefault(p1, p2, sym1, sym2, bondThickness, atomScale, style) {
    const bondPoints = getBondPoints(p1, p2, sym1, sym2, atomScale, style);
    if (!bondPoints) return null;
    const { startPos, endPos } = bondPoints;

    const midPoint = startPos.clone().add(endPos).multiplyScalar(0.5);
    const color1 = getAtomColorByScheme(sym1);
    const color2 = getAtomColorByScheme(sym2);

    const bondGroup = new THREE.Group();
    const bond1 = createHalfBond(startPos, midPoint, color1, bondThickness, style);
    const bond2 = createHalfBond(midPoint, endPos, color2, bondThickness, style);
    if (bond1) bondGroup.add(bond1);
    if (bond2) bondGroup.add(bond2);
    return bondGroup;
}

function createHalfBond(start, end, color, bondThickness, style) {
    if (start.distanceTo(end) <= 0) return null;
    const path = new THREE.LineCurve3(start, end);
    const geometry = new THREE.TubeGeometry(path, 1, bondThickness, 8, false);
    const material = (style === 'neon') ? new THREE.MeshBasicMaterial({ color: 0xffffff }) : getStandardMaterial(color, style);
    const bond = new THREE.Mesh(geometry, material);
    return bond;
}

function createBondStyleNeon(p1, p2, sym1, sym2, bondThickness, atomScale) {
    // Get atom info
    const info1 = atomInfo[sym1] || atomInfo['default'];
    const info2 = atomInfo[sym2] || atomInfo['default'];
    const radius1 = info1.radius * atomScale;
    const radius2 = info2.radius * atomScale;

    // Calculate bond vector
    const bondVector = new THREE.Vector3().subVectors(p2, p1);
    const bondLength = bondVector.length();
    const bondDirection = bondVector.clone().normalize();

    // Calculate start/end points from atom surfaces
    const startPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(radius1));
    const endPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(bondLength - radius2));

    // Return null if bond is too short (atoms overlap)
    const adjustedBondLength = bondLength - radius1 - radius2;
    if (adjustedBondLength <= 0) return null;

    const bondGroup = new THREE.Group();

    // Create cylindrical radial gradient effect with multiple layers

    // 1. Outermost layer - very transparent glow
    const path1 = new THREE.LineCurve3(startPos, endPos);
    const outerGeometry = new THREE.TubeGeometry(path1, 2, bondThickness * 1.2, 8, false);
    const outerMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.1, // Very transparent
        blending: THREE.AdditiveBlending,
        depthWrite: false
    });
    const outerTube = new THREE.Mesh(outerGeometry, outerMaterial);
    outerTube.renderOrder = -2;
    bondGroup.add(outerTube);

    // 2. Middle layer
    const path2 = new THREE.LineCurve3(startPos, endPos);
    const midGeometry = new THREE.TubeGeometry(path2, 2, bondThickness, 12, false);
    const midMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.3,
        blending: THREE.AdditiveBlending,
        depthWrite: false
    });
    const midTube = new THREE.Mesh(midGeometry, midMaterial);
    midTube.renderOrder = -1;
    bondGroup.add(midTube);

    // 3. Core layer - center line is brightest
    const path3 = new THREE.LineCurve3(startPos, endPos);
    const coreGeometry = new THREE.TubeGeometry(path3, 2, bondThickness * 0.8, 16, false);
    const coreMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.6, // Center is most opaque
        blending: THREE.AdditiveBlending,
        depthWrite: false
    });
    const coreTube = new THREE.Mesh(coreGeometry, coreMaterial);
    coreTube.renderOrder = 0;
    bondGroup.add(coreTube);

    return bondGroup;
}

function createBondStyleRowan(p1, p2, sym1, sym2, bondThickness, atomScale) {
    const bondPoints = getBondPoints(p1, p2, sym1, sym2, atomScale, 'rowan');
    if (!bondPoints) return null;
    const { startPos, endPos, dir } = bondPoints;
    const distance = startPos.distanceTo(endPos);

    const geometry = new THREE.CylinderGeometry(bondThickness, bondThickness, distance, 8);
    const material = new THREE.MeshBasicMaterial({ color: 0x000000 });
    const bond = new THREE.Mesh(geometry, material);
    bond.position.lerpVectors(startPos, endPos, 0.5);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
    return bond;
}

function createBondStyleBubble(p1, p2, sym1, sym2, bondThickness, atomScale) {
    const bondPoints = getBondPoints(p1, p2, sym1, sym2, atomScale, 'default');
    if (!bondPoints) return null;
    const { startPos, endPos } = bondPoints;

    const midPoint = startPos.clone().add(endPos).multiplyScalar(0.5);
    const color1 = getAtomColorByScheme(sym1);
    const color2 = getAtomColorByScheme(sym2);

    const bondGroup = new THREE.Group();

    // Create semi-transparent tubes with gradient-like appearance
    const bond1 = createHalfBondBubble(startPos, midPoint, color1, bondThickness);
    const bond2 = createHalfBondBubble(midPoint, endPos, color2, bondThickness);

    if (bond1) bondGroup.add(bond1);
    if (bond2) bondGroup.add(bond2);

    return bondGroup;
}

function createHalfBondBubble(start, end, color, bondThickness) {
    if (start.distanceTo(end) <= 0) return null;

    const path = new THREE.LineCurve3(start, end);
    const geometry = new THREE.TubeGeometry(path, 2, bondThickness, 12, false);
    const material = new THREE.MeshPhongMaterial({
        color: color,
        transparent: true,
        opacity: 0.4,  // More transparent for lighter appearance
        shininess: 80,
        specular: 0xffffff
    });
    const bond = new THREE.Mesh(geometry, material);

    return bond;
}

function createBondStyleGrey(p1, p2, sym1, sym2, bondThickness, atomScale, colorMap) {
    // Get atom info
    const info1 = atomInfo[sym1] || atomInfo['default'];
    const info2 = atomInfo[sym2] || atomInfo['default'];
    const radius1 = info1.radius * atomScale;
    const radius2 = info2.radius * atomScale;

    // Calculate bond vector
    const bondVector = new THREE.Vector3().subVectors(p2, p1);
    const bondLength = bondVector.length();
    const bondDirection = bondVector.clone().normalize();

    // Calculate start/end points from atom surfaces
    const startPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(radius1));
    const endPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(bondLength - radius2));

    // Return null if bond is too short (atoms overlap)
    const adjustedBondLength = bondLength - radius1 - radius2;
    if (adjustedBondLength <= 0) return null;

    // Calculate midpoint (based on atom surfaces)
    const midPoint = startPos.clone().add(endPos).multiplyScalar(0.5);

    const final_color1 = colorMap[sym1];
    const final_color2 = colorMap[sym2];

    const bondGroup = new THREE.Group();
    const bond1 = createHalfBondGrey(startPos, midPoint, final_color1, bondThickness);
    const bond2 = createHalfBondGrey(midPoint, endPos, final_color2, bondThickness);

    if (bond1) bondGroup.add(bond1);
    if (bond2) bondGroup.add(bond2);

    return bondGroup;
}

function createHalfBondGrey(start, end, color, bondThickness) {
    if (start.distanceTo(end) <= 0) return null;

    const distance = start.distanceTo(end);
    const direction = new THREE.Vector3().subVectors(end, start).normalize();

    // Use CylinderGeometry for cleaner appearance
    const geometry = new THREE.CylinderGeometry(bondThickness, bondThickness, distance, 16, 1);
    const material = new THREE.MeshBasicMaterial({
        color: color,
        depthWrite: true,
        depthTest: true
    });
    const bond = new THREE.Mesh(geometry, material);

    // Position and orient the cylinder
    bond.position.copy(start).add(end).multiplyScalar(0.5);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);

    return bond;
}

// ─────────────────────────────────────────────────────────────
// Convex Hull helper (incremental / gift-wrapping in 3D)
// Returns { vertices: Float32Array, indices: Uint32Array } for the hull.
// ─────────────────────────────────────────────────────────────
function _computeConvexHull(points) {
    // points: array of THREE.Vector3
    // Simplified incremental convex hull that works for common polyhedra (4-12 vertices)
    const n = points.length;
    if (n < 3) return null;
    if (n === 3) {
        // Single triangle
        return {
            faceVertexIndices: [[0, 1, 2]],
            faceColors: null
        };
    }

    // Gift-wrapping (Jarvis March) in 3D — find initial face, then expand
    // Step 1: Find two extreme points (min X, then min Y)
    let a = 0;
    for (let i = 1; i < n; i++) {
        if (points[i].x < points[a].x ||
            (points[i].x === points[a].x && points[i].y < points[a].y)) a = i;
    }

    // Step 2: Build faces using a BFS-based gift wrapping
    const faces = [];      // each face = [i, j, k] indices
    const edgeQueue = [];  // edges to process: { i, j, oppositeFace }
    const processedEdges = new Set();

    const edgeKey = (i, j) => i < j ? `${i}_${j}` : `${j}_${i}`;

    // Find initial triangle: pick point a, find b and c forming first face
    // b = point that makes smallest positive angle from (a, down-direction)
    let b = -1;
    for (let i = 0; i < n; i++) {
        if (i === a) continue;
        if (b === -1) { b = i; continue; }
        const ab = points[b].clone().sub(points[a]);
        const ac = points[i].clone().sub(points[a]);
        // prefer c over b if the cross product points "outward" (downward Z or smallest angle)
        if (ab.cross(ac).z < 0) b = i;
    }

    // c = point that makes a valid triangle with normals pointing "outward"
    let c = -1;
    for (let i = 0; i < n; i++) {
        if (i === a || i === b) continue;
        if (c === -1) { c = i; continue; }
        const ab = points[b].clone().sub(points[a]);
        const ac = points[c].clone().sub(points[a]);
        const ai = points[i].clone().sub(points[a]);
        const normal = ab.clone().cross(ac);
        // prefer i if it's more "extreme"
        if (normal.dot(ai) < 0) c = i;
    }

    if (c === -1) return null;

    // Ensure consistent winding (normal points away from centroid)
    const centroid = new THREE.Vector3();
    points.forEach(p => centroid.add(p));
    centroid.divideScalar(n);

    const faceNormal = (i, j, k) => {
        const ab = points[j].clone().sub(points[i]);
        const ac = points[k].clone().sub(points[i]);
        return ab.cross(ac);
    };

    let normal0 = faceNormal(a, b, c);
    const toCenter = centroid.clone().sub(points[a]);
    if (normal0.dot(toCenter) > 0) {
        // flip winding
        [b, c] = [c, b];
    }

    faces.push([a, b, c]);
    edgeQueue.push({ i: b, j: a, faceIdx: 0 });
    edgeQueue.push({ i: c, j: b, faceIdx: 0 });
    edgeQueue.push({ i: a, j: c, faceIdx: 0 });

    const maxIter = n * n * 3;
    let iter = 0;

    while (edgeQueue.length > 0 && iter++ < maxIter) {
        const { i, j } = edgeQueue.shift();
        const key = edgeKey(i, j);
        if (processedEdges.has(key)) continue;
        processedEdges.add(key);

        // Find the best k (visible from outside, minimum angle)
        const edge = points[j].clone().sub(points[i]);
        let bestK = -1;
        let bestAngle = Infinity;
        // Reference normal: (i -> j) x (i -> any existing face point on this edge's face)
        let refDir = null;
        // Find a point already in a face containing edge i-j (reversed: j-i)
        for (const f of faces) {
            const hasI = f.includes(i), hasJ = f.includes(j);
            if (hasI && hasJ) {
                const other = f.find(v => v !== i && v !== j);
                if (other !== undefined) {
                    refDir = points[other].clone().sub(points[i]);
                    break;
                }
            }
        }
        if (!refDir) continue;

        const refNormal = edge.clone().cross(refDir);

        for (let k = 0; k < n; k++) {
            if (k === i || k === j) continue;
            // skip if face i,j,k or j,i,k already exists
            const alreadyExists = faces.some(f => {
                const s = new Set(f);
                return s.has(i) && s.has(j) && s.has(k);
            });
            if (alreadyExists) continue;

            const toK = points[k].clone().sub(points[i]);
            const candidateNormal = edge.clone().cross(toK);
            // Angle from refNormal
            const dot = refNormal.dot(candidateNormal) / (refNormal.length() * candidateNormal.length() + 1e-10);
            const angle = Math.acos(Math.max(-1, Math.min(1, dot)));
            // We want the face that "wraps" — smallest positive angle (most "outward")
            if (angle < bestAngle) {
                // Check that centroid is on the inside
                const n_ = faceNormal(i, j, k);
                const toCentroid = centroid.clone().sub(points[i]);
                // Allow if normal points away from centroid
                if (n_.dot(toCentroid) < 0) {
                    bestAngle = angle;
                    bestK = k;
                }
            }
        }

        if (bestK === -1) continue;

        // Determine winding for new face
        let ni = i, nj = j, nk = bestK;
        const n_ = faceNormal(ni, nj, nk);
        const toCentroid = centroid.clone().sub(points[ni]);
        if (n_.dot(toCentroid) > 0) {
            [nj, nk] = [nk, nj];
        }

        const alreadyExists = faces.some(f => {
            const s = new Set(f);
            return s.has(ni) && s.has(nj) && s.has(nk);
        });
        if (!alreadyExists) {
            faces.push([ni, nj, nk]);
            edgeQueue.push({ i: nj, j: ni, faceIdx: faces.length - 1 });
            edgeQueue.push({ i: nk, j: nj, faceIdx: faces.length - 1 });
            edgeQueue.push({ i: ni, j: nk, faceIdx: faces.length - 1 });
        }
    }

    return faces.length >= 1 ? { faceVertexIndices: faces } : null;
}

/**
 * Create a convex polyhedron from neighbor atom positions around a center atom.
 * @param {THREE.Vector3} centerPos - Position of center atom (not a vertex)
 * @param {THREE.Vector3[]} neighborPositions - Positions of neighbor atoms (vertices)
 * @param {number[]} neighborColors - Colors (hex int) of neighbor atoms
 * @param {number} opacity - Face transparency (0-1)
 * @param {boolean} showEdge - Whether to render wireframe edges
 * @param {number} edgeOpacity - Edge line opacity (0-1)
 * @param {number} edgeColor - Edge color (hex int, default 0x111111)
 * @returns {THREE.Group|null}
 */
function createPolyhedronFaces(centerPos, neighborPositions, neighborColors, opacity = 0.25,
    showEdge = true, edgeOpacity = 0.7, edgeColor = 0x111111) {
    if (!neighborPositions || neighborPositions.length < 4) return null;

    const hull = _computeConvexHull(neighborPositions);
    if (!hull || !hull.faceVertexIndices || hull.faceVertexIndices.length === 0) return null;

    const faces = hull.faceVertexIndices;

    // Build BufferGeometry with vertex colors (each face gets averaged color of its 3 vertices)
    const vertexCount = faces.length * 3;
    const positions = new Float32Array(vertexCount * 3);
    const colors = new Float32Array(vertexCount * 3);

    let idx = 0;
    for (const [i, j, k] of faces) {
        const verts = [neighborPositions[i], neighborPositions[j], neighborPositions[k]];
        const cols = [neighborColors[i], neighborColors[j], neighborColors[k]];

        // Average color of 3 vertices
        const avgColor = new THREE.Color(
            ((cols[0] >> 16 & 0xff) + (cols[1] >> 16 & 0xff) + (cols[2] >> 16 & 0xff)) / (3 * 255),
            ((cols[0] >> 8 & 0xff) + (cols[1] >> 8 & 0xff) + (cols[2] >> 8 & 0xff)) / (3 * 255),
            ((cols[0] & 0xff) + (cols[1] & 0xff) + (cols[2] & 0xff)) / (3 * 255)
        );

        for (const v of verts) {
            positions[idx * 3] = v.x;
            positions[idx * 3 + 1] = v.y;
            positions[idx * 3 + 2] = v.z;
            colors[idx * 3] = avgColor.r;
            colors[idx * 3 + 1] = avgColor.g;
            colors[idx * 3 + 2] = avgColor.b;
            idx++;
        }
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: opacity,
        side: THREE.DoubleSide,
        depthWrite: false,
    });

    const group = new THREE.Group();
    group.add(new THREE.Mesh(geometry, material));

    // ── Wireframe edges ──────────────────────────────────────────
    if (showEdge) {
        const edgesGeo = new THREE.EdgesGeometry(geometry);
        const edgeMat = new THREE.LineBasicMaterial({
            color: edgeColor,
            transparent: edgeOpacity < 1.0,
            opacity: edgeOpacity,
            depthWrite: false,
        });
        group.add(new THREE.LineSegments(edgesGeo, edgeMat));
    }

    return group;
}

/**
 * Create a filled ring face (fan triangulation from centroid).
 * @param {THREE.Vector3[]} ringPositions - Ordered ring atom positions
 * @param {number} faceColor - Hex color int
 * @param {number} opacity - Transparency (0-1)
 * @returns {THREE.Mesh|null}
 */
function createRingFace(ringPositions, faceColor, opacity = 0.3) {
    if (!ringPositions || ringPositions.length < 3) return null;

    const n = ringPositions.length;
    // Centroid of ring
    const centroid = new THREE.Vector3();
    ringPositions.forEach(p => centroid.add(p));
    centroid.divideScalar(n);

    // Fan triangulation: centroid -> (p[i], p[i+1])
    const vertCount = n * 3;
    const positions = new Float32Array(vertCount * 3);
    let idx = 0;

    for (let i = 0; i < n; i++) {
        const p1 = ringPositions[i];
        const p2 = ringPositions[(i + 1) % n];
        positions[idx++] = centroid.x; positions[idx++] = centroid.y; positions[idx++] = centroid.z;
        positions[idx++] = p1.x; positions[idx++] = p1.y; positions[idx++] = p1.z;
        positions[idx++] = p2.x; positions[idx++] = p2.y; positions[idx++] = p2.z;
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshBasicMaterial({
        color: faceColor,
        transparent: true,
        opacity: opacity,
        side: THREE.DoubleSide,
        depthWrite: false,
    });

    return new THREE.Mesh(geometry, material);
}

// ═══════════════════════════════════════════════════════════════
// Performance Optimization Utilities
// ═══════════════════════════════════════════════════════════════

/**
 * Determine LOD (Level of Detail) sphere segments based on atom count
 * and performance mode setting.
 *
 * @param {number} atomCount - Total number of atoms to render
 * @param {string} performanceMode - "auto", "high", or "normal"
 * @returns {{ sphere: number, sphereHigh: number, bond: number }}
 */
function getLODSegments(atomCount, performanceMode) {
    if (performanceMode === 'normal') {
        return { sphere: 64, sphereHigh: 128, bond: 8 };
    }

    const isHigh = performanceMode === 'high' ||
                   (performanceMode === 'auto' && atomCount > 100);

    if (!isHigh) {
        return { sphere: 64, sphereHigh: 128, bond: 8 };
    }

    // Adaptive: reduce segments as atom count grows
    if (atomCount > 2000) {
        return { sphere: 8, sphereHigh: 16, bond: 4 };
    } else if (atomCount > 500) {
        return { sphere: 12, sphereHigh: 24, bond: 6 };
    } else {
        // 100 < atomCount <= 500
        return { sphere: 16, sphereHigh: 32, bond: 6 };
    }
}

/**
 * Geometry & material cache — avoids creating duplicate SphereGeometry
 * and material instances for atoms of the same element and style.
 */
const _geoCache = {};
const _matCache = {};

function clearPerformanceCache() {
    // Dispose cached geometries
    for (const key in _geoCache) {
        if (_geoCache[key] && _geoCache[key].dispose) _geoCache[key].dispose();
    }
    for (const key in _matCache) {
        if (_matCache[key] && _matCache[key].dispose) _matCache[key].dispose();
    }
    for (const key in _geoCache) delete _geoCache[key];
    for (const key in _matCache) delete _matCache[key];
}

function getCachedSphereGeometry(radius, segments) {
    const key = `sphere_${radius.toFixed(4)}_${segments}`;
    if (!_geoCache[key]) {
        _geoCache[key] = new THREE.SphereGeometry(radius, segments, segments);
        _geoCache[key]._isCached = true;
    }
    return _geoCache[key];
}

function getCachedMaterial(type, colorHex, extraKey) {
    const key = `mat_${type}_${colorHex}_${extraKey || ''}`;
    if (!_matCache[key]) {
        switch (type) {
            case 'toon': {
                const colors = new Uint8Array(2);
                colors[0] = 0; colors[1] = 255;
                const gradientMap = new THREE.DataTexture(colors, 2, 1, THREE.LuminanceFormat);
                gradientMap.minFilter = THREE.NearestFilter;
                gradientMap.magFilter = THREE.NearestFilter;
                gradientMap.needsUpdate = true;
                _matCache[key] = new THREE.MeshToonMaterial({
                    color: new THREE.Color(colorHex),
                    gradientMap: gradientMap
                });
                break;
            }
            case 'outline':
                _matCache[key] = new THREE.MeshBasicMaterial({
                    color: colorHex,
                    side: THREE.BackSide
                });
                break;
            case 'glossy':
                _matCache[key] = new THREE.MeshPhongMaterial({
                    color: colorHex, shininess: 100, specular: 0x222222
                });
                break;
            case 'metallic':
                _matCache[key] = new THREE.MeshStandardMaterial({
                    color: colorHex, metalness: 0.3, roughness: 0.4
                });
                break;
            case 'lambert':
                _matCache[key] = new THREE.MeshLambertMaterial({ color: colorHex });
                break;
            case 'lambertFlat':
                _matCache[key] = new THREE.MeshLambertMaterial({
                    color: new THREE.Color(colorHex), flatShading: false
                });
                break;
            default:
                _matCache[key] = new THREE.MeshLambertMaterial({ color: colorHex });
        }
        _matCache[key]._isCached = true;
    }
    return _matCache[key];
}

/**
 * Spatial grid for O(n) bond detection.
 * Instead of O(n²) all-pairs, partition atoms into grid cells
 * and only check atoms in neighboring cells.
 *
 * @param {THREE.Vector3[]} positions - Atom positions
 * @param {string[]} symbols - Atom symbols
 * @param {number} bondCutoff - Bond threshold multiplier for covalent radii
 * @returns {Array<{i: number, j: number}>} - Bond pairs
 */
function findBondsSpatialGrid(positions, symbols, bondCutoff) {
    const n = positions.length;
    if (n === 0) return [];

    // For small systems, brute force is fine and avoids grid overhead
    if (n < 80) {
        return _findBondsBruteForce(positions, symbols, bondCutoff);
    }

    // Determine max possible bond distance
    let maxCutoff = 0;
    for (let i = 0; i < n; i++) {
        const r = (atomInfo[symbols[i]] || atomInfo['default']).radius;
        const d = r * 2 * bondCutoff;
        if (d > maxCutoff) maxCutoff = d;
    }

    // Build spatial grid
    const cellSize = maxCutoff + 0.001;
    const grid = {};
    const bonds = [];

    function cellKey(x, y, z) {
        return `${x},${y},${z}`;
    }

    // Insert atoms into grid cells
    for (let i = 0; i < n; i++) {
        const p = positions[i];
        const cx = Math.floor(p.x / cellSize);
        const cy = Math.floor(p.y / cellSize);
        const cz = Math.floor(p.z / cellSize);
        const key = cellKey(cx, cy, cz);
        if (!grid[key]) grid[key] = [];
        grid[key].push(i);
    }

    // Check neighboring cells (3x3x3 neighborhood)
    const checked = new Set();
    for (let i = 0; i < n; i++) {
        const p = positions[i];
        const cx = Math.floor(p.x / cellSize);
        const cy = Math.floor(p.y / cellSize);
        const cz = Math.floor(p.z / cellSize);

        for (let dx = -1; dx <= 1; dx++) {
            for (let dy = -1; dy <= 1; dy++) {
                for (let dz = -1; dz <= 1; dz++) {
                    const key = cellKey(cx + dx, cy + dy, cz + dz);
                    const cell = grid[key];
                    if (!cell) continue;

                    for (let k = 0; k < cell.length; k++) {
                        const j = cell[k];
                        if (j <= i) continue;

                        const pairKey = i * n + j;
                        if (checked.has(pairKey)) continue;
                        checked.add(pairKey);

                        const dist = positions[i].distanceTo(positions[j]);
                        const r1 = (atomInfo[symbols[i]] || atomInfo['default']).radius;
                        const r2 = (atomInfo[symbols[j]] || atomInfo['default']).radius;
                        const adjustedCutoff = (r1 + r2) * bondCutoff;

                        if (dist <= adjustedCutoff) {
                            bonds.push({ i: i, j: j });
                        }
                    }
                }
            }
        }
    }

    return bonds;
}

function _findBondsBruteForce(positions, symbols, bondCutoff) {
    const bonds = [];
    const n = positions.length;
    for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
            const dist = positions[i].distanceTo(positions[j]);
            const r1 = (atomInfo[symbols[i]] || atomInfo['default']).radius;
            const r2 = (atomInfo[symbols[j]] || atomInfo['default']).radius;
            if (dist <= (r1 + r2) * bondCutoff) {
                bonds.push({ i: i, j: j });
            }
        }
    }
    return bonds;
}

// ======== Hydrogen Bond Functions ========

/**
 * Detect hydrogen bonds using the D-H...A geometric criteria.
 *
 * Donors:    H bonded to N, O, or F
 * Acceptors: N, O, F
 * Criteria:
 *   - H...A distance  ≤ 2.5 Å
 *   - D...A distance  ≤ 3.5 Å
 *   - D-H...A angle   ≥ 120°
 * Excluded:
 *   - D and A directly bonded (1,2)
 *   - D and A share a common bonded neighbor (1,3: H-D-X-A)
 *
 * @param {THREE.Vector3[]} positions   - Atom positions
 * @param {string[]}        symbols     - Element symbols
 * @param {Array<{i,j}>}    covalentBonds - Existing covalent bond pairs
 * @returns {Array<{h, a}>} H-bond pairs: index of H and acceptor atom
 */
function detectHydrogenBonds(positions, symbols, covalentBonds, options = {}) {
    const donorElements   = new Set(['N', 'O', 'F']);
    const acceptorElements = new Set(['N', 'O', 'F']);

    const HA_CUTOFF = Number.isFinite(options.hBondThreshold) ? options.hBondThreshold : 2.5;   // Å  H...A distance
    const DA_CUTOFF = 3.5;   // Å  D...A distance
    const ANGLE_MIN = 120.0; // °  D-H...A angle

    // Build adjacency list and direct-bond set from covalent bonds
    const adjacency = {};
    const bondSet   = new Set();
    if (covalentBonds) {
        covalentBonds.forEach(({ i, j }) => {
            if (!adjacency[i]) adjacency[i] = [];
            if (!adjacency[j]) adjacency[j] = [];
            adjacency[i].push(j);
            adjacency[j].push(i);
            bondSet.add(`${Math.min(i, j)}_${Math.max(i, j)}`);
        });
    }

    const hbonds = [];

    for (let h = 0; h < symbols.length; h++) {
        if (symbols[h] !== 'H') continue;

        const hNeighbors = adjacency[h] || [];
        // H must be covalently bonded to a strong donor (N/O/F); excludes C-H, S-H
        const donorIdx = hNeighbors.find(d => donorElements.has(symbols[d]));
        if (donorIdx === undefined) continue;

        const pH = positions[h];
        const pD = positions[donorIdx];
        const dNeighbors = adjacency[donorIdx] || [];

        for (let a = 0; a < symbols.length; a++) {
            if (a === h || a === donorIdx) continue;

            if (!acceptorElements.has(symbols[a])) continue;

            const haCutoff = HA_CUTOFF;
            const daCutoff = DA_CUTOFF;
            const angleMin = ANGLE_MIN;

            // --- Topology filters (cheap, no sqrt) ---

            // Exclude 1,2: D directly bonded to A
            if (bondSet.has(`${Math.min(donorIdx, a)}_${Math.max(donorIdx, a)}`)) continue;

            // Exclude 1,3: D-X-A (H-D-X-A = 3-bond path), e.g. H-N-C=O
            const is13 = dNeighbors.some(
                x => x !== h && bondSet.has(`${Math.min(x, a)}_${Math.max(x, a)}`)
            );
            if (is13) continue;

            // --- Distance filters ---
            const pA = positions[a];

            const haDist = pH.distanceTo(pA);
            if (haDist > haCutoff) continue;

            const daDist = pD.distanceTo(pA);
            if (daDist > daCutoff) continue;

            // --- Angle filter: D-H...A angle at H ---
            const vecHD = new THREE.Vector3().subVectors(pD, pH).normalize();
            const vecHA = new THREE.Vector3().subVectors(pA, pH).normalize();
            const cosAngle = Math.max(-1, Math.min(1, vecHD.dot(vecHA)));
            const angleDeg = Math.acos(cosAngle) * (180 / Math.PI);
            if (angleDeg < angleMin) continue;

            hbonds.push({ h, a });
        }
    }

    return hbonds;
}

/**
 * Create a dashed-cylinder representation of a hydrogen bond.
 *
 * @param {THREE.Vector3} p1
 * @param {THREE.Vector3} p2
 * @returns {THREE.Group}
 */
function createHydrogenBond(p1, p2) {
    const distance  = p1.distanceTo(p2);
    const direction = new THREE.Vector3().subVectors(p2, p1).normalize();

    const hBondRadius = 0.05;  // thinner than covalent bonds (0.1)
    const dashLength  = 0.20;
    const gapLength   = 0.12;

    const group = new THREE.Group();
    const up    = new THREE.Vector3(0, 1, 0);

    let traveled = 0;
    let isDash   = true;

    while (traveled < distance) {
        const segLen   = isDash ? dashLength : gapLength;
        const actual   = Math.min(segLen, distance - traveled);

        if (isDash && actual > 0.01) {
            const start   = p1.clone().addScaledVector(direction, traveled);
            const end     = p1.clone().addScaledVector(direction, traveled + actual);
            const midPt   = start.clone().lerp(end, 0.5);

            const geometry = new THREE.CylinderGeometry(hBondRadius, hBondRadius, actual, 6, 1);
            const material = new THREE.MeshBasicMaterial({
                color:       0x222222,
                transparent: true,
                opacity:     0.65,
                depthTest:   true,
                depthWrite:  false
            });

            const cylinder = new THREE.Mesh(geometry, material);
            cylinder.position.copy(midPt);
            cylinder.quaternion.setFromUnitVectors(up, direction);
            group.add(cylinder);
        }

        traveled += actual;
        isDash = !isDash;
    }

    return group;
}
