// Overlay Viewer - renderMolecule function
// This renders all molecules simultaneously with different colors

function renderMolecule() {
    if (!scene || !camera || !renderer) return;
    if (!molecularData || molecularData.length === 0) return;

    if (moleculeGroup) {
        scene.remove(moleculeGroup);
        disposeObject3D(moleculeGroup);
    }

    moleculeGroup = new THREE.Group();
    
    // Color palette for different molecules
    const moleculeColors = [
        0xff6b6b, // Red
        0x4ecdc4, // Cyan  
        0xffe66d, // Yellow
        0x95e1d3, // Mint
        0xf38181, // Pink
        0xaa96da, // Purple
        0xfcbad3, // Light pink
        0xa8e6cf, // Light green
    ];
    
    const styleName = (settings.style || 'cartoon').toLowerCase();
    const atomScale = Math.max(settings.atomSize, 0.05);
    const bondRadius = Math.max(settings.bondThickness, 0.005);
    const bondCutoff = settings.bondThreshold;
    
    // Render each molecule
    molecularData.forEach((molecule, molIndex) => {
        if (!molecule || !molecule.positions || !molecule.symbols) return;
        
        const positions = molecule.positions.map(pos => new THREE.Vector3(...pos));
        const symbols = molecule.symbols;
        const moleculeColor = moleculeColors[molIndex % moleculeColors.length];
        const opacity = 0.8; // Can be controlled by UI later
        
        // Create molecule group
        const molGroup = new THREE.Group();
        
        // Render bonds
        if (settings.showBond) {
            for (let i = 0; i < positions.length; i++) {
                for (let j = i + 1; j < positions.length; j++) {
                    const dist = positions[i].distanceTo(positions[j]);
                    if (dist <= bondCutoff) {
                        const bondGeometry = new THREE.CylinderGeometry(bondRadius, bondRadius, dist, 8);
                        const bondMaterial = new THREE.MeshPhongMaterial({
                            color: moleculeColor,
                            transparent: true,
                            opacity: opacity
                        });
                        const bond = new THREE.Mesh(bondGeometry, bondMaterial);
                        
                        // Position and orient bond
                        const midPoint = positions[i].clone().add(positions[j]).multiplyScalar(0.5);
                        bond.position.copy(midPoint);
                        const direction = positions[j].clone().sub(positions[i]).normalize();
                        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
                        
                        molGroup.add(bond);
                    }
                }
            }
        }
        
        // Render atoms
        for (let i = 0; i < positions.length; i++) {
            const atomGeometry = new THREE.SphereGeometry(atomScale * 0.5, 32, 32);
            const atomMaterial = new THREE.MeshPhongMaterial({
                color: moleculeColor,
                transparent: true,
                opacity: opacity
            });
            const atom = new THREE.Mesh(atomGeometry, atomMaterial);
            atom.position.copy(positions[i]);
            molGroup.add(atom);
        }
        
        moleculeGroup.add(molGroup);
    });
    
    scene.add(moleculeGroup);
    
    if (cameraNeedsFit) {
        fitCameraToMolecule(moleculeGroup);
    }
}
