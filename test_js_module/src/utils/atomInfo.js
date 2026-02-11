/**
 * Atom information: colors and radii for common elements
 */
export const atomInfo = {
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
export function getAtomInfo(symbol) {
    return atomInfo[symbol] || atomInfo['default'];
}

/**
 * Detect bonds between atoms based on covalent radii
 */
export function detectBonds(positions, symbols, threshold = 1.2) {
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
