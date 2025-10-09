
const atomInfo = {
    'H': { radius: 0.37, color: 0xFFFFFF }, 'He': { radius: 0.32, color: 0xD9FFFF },
    'Li': { radius: 1.34, color: 0xCC80FF }, 'Be': { radius: 0.90, color: 0xC2FF00 },
    'B': { radius: 0.82, color: 0xFFB5B5 }, 'C': { radius: 0.77, color: 0x909090 },
    'N': { radius: 0.75, color: 0x3050F8 }, 'O': { radius: 0.73, color: 0xFF0D0D },
    'F': { radius: 0.71, color: 0x90E050 }, 'Ne': { radius: 0.69, color: 0xB3E3F5 },
    'Na': { radius: 1.54, color: 0xAB5CF2 }, 'Mg': { radius: 1.30, color: 0x8AFF00 },
    'Al': { radius: 1.18, color: 0xBFA6A6 }, 'Si': { radius: 1.11, color: 0xF0C8A0 },
    'P': { radius: 1.06, color: 0xFF8000 }, 'S': { radius: 1.02, color: 0xFFFF30 },
    'Cl': { radius: 0.99, color: 0x1FF01F }, 'Ar': { radius: 0.97, color: 0x80D1E3 },
    'K': { radius: 1.96, color: 0x8F40D4 }, 'Ca': { radius: 1.74, color: 0x3DFF00 },
    'Sc': { radius: 1.44, color: 0xE6E6E6 }, 'Ti': { radius: 1.36, color: 0xBFC2C7 },
    'V': { radius: 1.25, color: 0xA6A6AB }, 'Cr': { radius: 1.27, color: 0x8A99C7 },
    'Mn': { radius: 1.39, color: 0x9C7AC7 }, 'Fe': { radius: 1.25, color: 0xE06633 },
    'Co': { radius: 1.26, color: 0xF090A0 }, 'Ni': { radius: 1.21, color: 0x50D050 },
    'Cu': { radius: 1.38, color: 0xC88033 }, 'Zn': { radius: 1.31, color: 0x7D80B0 },
    'Ga': { radius: 1.26, color: 0xC28F8F }, 'Ge': { radius: 1.22, color: 0x668F8F },
    'As': { radius: 1.19, color: 0xBD80E3 }, 'Se': { radius: 1.16, color: 0xFFA100 },
    'Br': { radius: 1.14, color: 0xA62929 }, 'Kr': { radius: 1.10, color: 0x5CB8D1 },
    'Rb': { radius: 2.11, color: 0x702EB0 }, 'Sr': { radius: 1.92, color: 0x00FF00 },
    'Y': { radius: 1.62, color: 0x94FFFF }, 'Zr': { radius: 1.48, color: 0x94E0E0 },
    'Nb': { radius: 1.37, color: 0x73C2C9 }, 'Mo': { radius: 1.45, color: 0x54B5B5 },
    'Tc': { radius: 1.56, color: 0x3B9E9E }, 'Ru': { radius: 1.26, color: 0x248F8F },
    'Rh': { radius: 1.35, color: 0x0A7D8C }, 'Pd': { radius: 1.31, color: 0x006985 },
    'Ag': { radius: 1.53, color: 0xC0C0C0 }, 'Cd': { radius: 1.48, color: 0xFFD98F },
    'In': { radius: 1.44, color: 0xA67575 }, 'Sn': { radius: 1.41, color: 0x668080 },
    'Sb': { radius: 1.38, color: 0x9E63B5 }, 'Te': { radius: 1.35, color: 0xD47A00 },
    'I': { radius: 1.33, color: 0x940094 }, 'Xe': { radius: 1.30, color: 0x429EB0 },
    'Cs': { radius: 2.25, color: 0x57178F }, 'Ba': { radius: 1.98, color: 0x00C900 },
    'La': { radius: 1.69, color: 0x70D4FF }, 'Ce': { radius: 1.65, color: 0xFFFFC7 },
    'Pr': { radius: 1.65, color: 0xD9FFC7 }, 'Nd': { radius: 1.64, color: 0xC7FFC7 },
    'Pm': { radius: 1.63, color: 0xA3FFC7 }, 'Sm': { radius: 1.62, color: 0x8FFFC7 },
    'Eu': { radius: 1.61, color: 0x61FFC7 }, 'Gd': { radius: 1.59, color: 0x45FFC7 },
    'Tb': { radius: 1.58, color: 0x30FFC7 }, 'Dy': { radius: 1.59, color: 0x1FFFC7 },
    'Ho': { radius: 1.58, color: 0x00FF9C }, 'Er': { radius: 1.57, color: 0x00E675 },
    'Tm': { radius: 1.56, color: 0x00D452 }, 'Yb': { radius: 1.56, color: 0x00BF38 },
    'Lu': { radius: 1.57, color: 0x00AB24 }, 'Hf': { radius: 1.50, color: 0x4DC2FF },
    'Ta': { radius: 1.38, color: 0x4DA6FF }, 'W': { radius: 1.46, color: 0x2194D6 },
    'Re': { radius: 1.59, color: 0x267DAB }, 'Os': { radius: 1.28, color: 0x266696 },
    'Ir': { radius: 1.37, color: 0x175487 }, 'Pt': { radius: 1.28, color: 0xD0D0E0 },
    'Au': { radius: 1.44, color: 0xFFD123 }, 'Hg': { radius: 1.49, color: 0xB8B8D0 },
    'Tl': { radius: 1.48, color: 0xA6544D }, 'Pb': { radius: 1.47, color: 0x575961 },
    'Bi': { radius: 1.46, color: 0x9E4FB5 }, 'Po': { radius: 1.45, color: 0xAB5C00 },
    'At': { radius: 1.47, color: 0x754F45 }, 'Rn': { radius: 1.50, color: 0x428296 },
    'Fr': { radius: 2.23, color: 0x420066 }, 'Ra': { radius: 2.01, color: 0x007D00 },
    'Ac': { radius: 1.86, color: 0x70ABFA }, 'Th': { radius: 1.75, color: 0x00BAFF },
    'Pa': { radius: 1.69, color: 0x00A1FF }, 'U': { radius: 1.70, color: 0x008FFF },
    'Np': { radius: 1.71, color: 0x0080FF }, 'Pu': { radius: 1.72, color: 0x006BFF },
    'Am': { radius: 1.66, color: 0x545CF2 }, 'Cm': { radius: 1.66, color: 0x785CE3 },
    'Bk': { radius: 1.68, color: 0x8A4FE3 }, 'Cf': { radius: 1.68, color: 0xA136D4 },
    'Es': { radius: 1.65, color: 0xB31FD4 }, 'Fm': { radius: 1.67, color: 0xB31FBA },
    'Md': { radius: 1.73, color: 0xB30DA6 }, 'No': { radius: 1.76, color: 0xBD0D87 },
    'Lr': { radius: 1.61, color: 0xC70066 }, 'Rf': { radius: 1.57, color: 0xCC0059 },
    'Db': { radius: 1.49, color: 0xD1004F }, 'Sg': { radius: 1.43, color: 0xD90045 },
    'Bh': { radius: 1.41, color: 0xE00038 }, 'Hs': { radius: 1.34, color: 0xE6002E },
    'Mt': { radius: 1.29, color: 0xEB0026 },
    'default': { radius: 0.7, color: 0xFF00FF } // Default for unknown elements
};

function getStandardMaterial(color, type) {
    switch(type) {
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

function createAtomStyleCartoon(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    
    // Create simple 2-step gradient map for clean toon shading
    const colors = new Uint8Array(2);
    colors[0] = 0;      // Dark
    colors[1] = 255;    // Light
    const gradientMap = new THREE.DataTexture(colors, 2, 1, THREE.LuminanceFormat);
    gradientMap.minFilter = THREE.NearestFilter;
    gradientMap.magFilter = THREE.NearestFilter;
    gradientMap.needsUpdate = true;
    
    const geometry = new THREE.SphereGeometry(scaledRadius, 64, 64);  // Increased from 32 to 64
    const toonMaterial = new THREE.MeshToonMaterial({
        color: new THREE.Color(info.color),
        gradientMap: gradientMap
    });
    const sphere = new THREE.Mesh(geometry, toonMaterial);

    const outlineGeometry = new THREE.SphereGeometry(scaledRadius * 1.05, 64, 64);  // Increased from 32 to 64
    const outlineMaterial = new THREE.MeshBasicMaterial({
        color: 0x000000,
        side: THREE.BackSide
    });
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
    const atomColor = new THREE.Color(info.color);
    
    // 단일 sprite로 glow 효과만 사용
    const glowCanvas = document.createElement('canvas');
    glowCanvas.width = 512;
    glowCanvas.height = 512;
    const glowContext = glowCanvas.getContext('2d');
    
    const center = 256;
    const maxRadius = 256;
    
    // 자연스러운 radial gradient - 40% 지점부터 서서히 증가, 테두리에서 0
    const glowGradient = glowContext.createRadialGradient(center, center, 0, center, center, maxRadius);
    const glowColorStyle = atomColor.getStyle();
    
	glowGradient.addColorStop(0.0, 'rgba(0,0,0,0)');                            // 중심: 완전 투명
	glowGradient.addColorStop(0.25, 'rgba(0,0,0,0)');                           // 25%까지 투명 유지
	glowGradient.addColorStop(0.35, `${glowColorStyle.slice(0, -1)}, 0.05)`); // 서서히 시작
	glowGradient.addColorStop(0.45, `${glowColorStyle.slice(0, -1)}, 0.15)`); 
	glowGradient.addColorStop(0.60, `${glowColorStyle.slice(0, -1)}, 0.35)`); 
	glowGradient.addColorStop(0.75, `${glowColorStyle.slice(0, -1)}, 0.6)`);  
	glowGradient.addColorStop(0.88, `${glowColorStyle.slice(0, -1)}, 0.85)`); 
	glowGradient.addColorStop(0.95, `${glowColorStyle.slice(0, -1)}, 1.0)`);  // 가장자리 근처 최대치
	glowGradient.addColorStop(1.0, `${glowColorStyle.slice(0, -1)}, 0.0)`);   // 부드럽게 사라지기

    
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
    // 원자 테두리까지 정확히 맞춤
    glowSprite.scale.set(scaledRadius * 2.2, scaledRadius * 2.2, 1);
    
    // 메타데이터 추가
    glowSprite.userData = {
        symbol: symbol,
        radius: scaledRadius,
        originalColor: atomColor.clone()
    };
    
    return glowSprite;
}

function createAtomStyleGlossy(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const geometry = new THREE.SphereGeometry(scaledRadius, 64, 64);
    const material = getStandardMaterial(info.color, 'glossy');
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleMetallic(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || info['default'];
    const scaledRadius = info.radius * atomScale;
    const geometry = new THREE.SphereGeometry(scaledRadius, 64, 64);
    const material = getStandardMaterial(info.color, 'metallic');
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleDefault(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const geometry = new THREE.SphereGeometry(scaledRadius, 64, 64);
    const material = getStandardMaterial(info.color, 'default');
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(pos);
    return sphere;
}

function createAtomStyleRowan(pos, symbol, atomScale) {
    const info = atomInfo[symbol] || atomInfo['default'];
    const scaledRadius = info.radius * atomScale;
    const group = new THREE.Group();

    // Use higher resolution to avoid banding artifacts
    const geometry = new THREE.SphereGeometry(scaledRadius, 128, 128);
    const material = new THREE.MeshLambertMaterial({
        color: new THREE.Color(info.color),
        flatShading: false
    });
    const sphere = new THREE.Mesh(geometry, material);
    group.add(sphere);

    const outlineGeometry = new THREE.SphereGeometry(scaledRadius * 1.05, 64, 64);
    const outlineMaterial = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.BackSide });
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
    
    const atomColor = new THREE.Color(info.color);
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
    
    // 3D 구체는 제거하고, 캔버스로 그린 스프라이트만 사용
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    canvas.width = 256;
    canvas.height = 256;
    
    const centerX = 128;
    const centerY = 128;
    const circleRadius = 120;
    
    // 검은 테두리 그리기
    context.beginPath();
    context.arc(centerX, centerY, circleRadius, 0, 2 * Math.PI, false);
    context.fillStyle = 'black';
    context.fill();
    
    // 원형 배경 그리기 (구체 역할) - 테두리보다 약간 작게
    context.beginPath();
    context.arc(centerX, centerY, circleRadius - 3, 0, 2 * Math.PI, false);
    context.fillStyle = `#${color.toString(16).padStart(6, '0')}`; // color를 hex로 변환
    context.fill();
    
    // 약간 더 밝은 색으로 하이라이트 효과
    context.beginPath();
    context.arc(centerX, centerY, circleRadius - 7, 0, 2 * Math.PI, false);
    const lighterColor = new THREE.Color(color).multiplyScalar(1.1);
    context.fillStyle = `#${lighterColor.getHexString()}`;
    context.fill();
    
    // 텍스트 그리기
    context.font = 'bold 80px Arial';
    context.fillStyle = 'white'; // 텍스트를 흰색으로
    context.strokeStyle = 'black';
    context.lineWidth = 2;
    context.textAlign = 'center';
    context.textBaseline = 'middle';
    
    // 텍스트에 외곽선 효과
    context.strokeText(symbol, centerX, centerY);
    context.fillText(symbol, centerX, centerY);
    
    const texture = new THREE.CanvasTexture(canvas);
    const spriteMaterial = new THREE.SpriteMaterial({
        map: texture,
        transparent: true,
        depthTest: true,     // 깊이 테스트 활성화
        depthWrite: false,    // 투명 객체이므로 깊이 쓰기는 비활성화
        alphaTest: 0.5
    });
    
    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.position.copy(pos); // 주어진 위치에 직접 배치
    sprite.scale.set(scaledRadius * 2.5, scaledRadius * 2.5, 1);
    sprite.center.set(0.5, 0.5);
    
    return sprite; // Group이 아닌 sprite 직접 반환
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
    const color1 = (atomInfo[sym1] || atomInfo.default).color;
    const color2 = (atomInfo[sym2] || atomInfo.default).color;

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
    // 원자 정보 가져오기
    const info1 = atomInfo[sym1] || atomInfo['default'];
    const info2 = atomInfo[sym2] || atomInfo['default'];
    const radius1 = info1.radius * atomScale;
    const radius2 = info2.radius * atomScale;
    
    // 결합 벡터 계산
    const bondVector = new THREE.Vector3().subVectors(p2, p1);
    const bondLength = bondVector.length();
    const bondDirection = bondVector.clone().normalize();
    
    // 원자 표면에서 시작하고 끝나는 점 계산
    const startPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(radius1));
    const endPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(bondLength - radius2));
    
    // 결합이 너무 짧으면 (원자가 겹치면) null 반환
    const adjustedBondLength = bondLength - radius1 - radius2;
    if (adjustedBondLength <= 0) return null;
    
    const bondGroup = new THREE.Group();
    
    // 다중 레이어로 원통 radial gradient 효과 생성
    
    // 1. 가장 바깥쪽 레이어 - 매우 투명한 glow
    const path1 = new THREE.LineCurve3(startPos, endPos);
    const outerGeometry = new THREE.TubeGeometry(path1, 2, bondThickness * 1.2, 8, false);
    const outerMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.1, // 매우 투명
        blending: THREE.AdditiveBlending,
        depthWrite: false
    });
    const outerTube = new THREE.Mesh(outerGeometry, outerMaterial);
    outerTube.renderOrder = -2;
    bondGroup.add(outerTube);
    
    // 2. 중간 레이어
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
    
    // 3. 코어 레이어 - 중심선이 가장 밝음 (원래 코드의 오타 수정)
    const path3 = new THREE.LineCurve3(startPos, endPos);
    const coreGeometry = new THREE.TubeGeometry(path3, 2, bondThickness * 0.8, 16, false);
    const coreMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.6, // 중심이 가장 불투명
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
    const color1 = (atomInfo[sym1] || atomInfo.default).color;
    const color2 = (atomInfo[sym2] || atomInfo.default).color;

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
    // 원자 정보 가져오기
    const info1 = atomInfo[sym1] || atomInfo['default'];
    const info2 = atomInfo[sym2] || atomInfo['default'];
    const radius1 = info1.radius * atomScale;
    const radius2 = info2.radius * atomScale;
    
    // 결합 벡터 계산
    const bondVector = new THREE.Vector3().subVectors(p2, p1);
    const bondLength = bondVector.length();
    const bondDirection = bondVector.clone().normalize();
    
    // 원자 표면에서 시작하고 끝나는 점 계산
    const startPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(radius1));
    const endPos = new THREE.Vector3().addVectors(p1, bondDirection.clone().multiplyScalar(bondLength - radius2));
    
    // 결합이 너무 짧으면 (원자가 겹치면) null 반환
    const adjustedBondLength = bondLength - radius1 - radius2;
    if (adjustedBondLength <= 0) return null;
    
    // 중점 계산 (원자 표면 기준)
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