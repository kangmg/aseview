import json
import numpy as np
import os
import pkg_resources
from typing import Dict, Any, List, Optional, Union
from ase import Atoms


class MolecularData:
    """Class to handle molecular data conversion between ASE Atoms and JSON format."""
    
    @staticmethod
    def from_atoms(atoms: Atoms) -> Dict[str, Any]:
        """Convert ASE Atoms object to JSON-serializable dictionary."""
        positions = atoms.get_positions().tolist()
        symbols = atoms.get_chemical_symbols()
        
        data = {
            "positions": positions,
            "symbols": symbols
        }
        
        # Add cell information if present
        if atoms.pbc.any():
            data["cell"] = atoms.get_cell().tolist()
            
        # Add forces if present
        if hasattr(atoms, 'arrays') and 'forces' in atoms.arrays:
            data["forces"] = atoms.arrays['forces'].tolist()
            
        return data
    
    @staticmethod
    def to_atoms(data: Dict[str, Any]) -> Atoms:
        """Convert JSON-serializable dictionary to ASE Atoms object."""
        atoms = Atoms(symbols=data["symbols"], positions=data["positions"])
        
        if "cell" in data:
            atoms.set_cell(data["cell"])
            atoms.pbc = [True, True, True]
            
        return atoms


class BaseViewer:
    """Base class for all molecular viewers."""
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str]):
        """
        Initialize the viewer with molecular data.
        
        Args:
            data: Can be an ASE Atoms object, a dictionary with molecular data,
                  or a path to a JSON file containing molecular data.
        """
        if isinstance(data, str):
            # Assume it's a path to a JSON file
            with open(data, 'r') as f:
                self.data = json.load(f)
        elif isinstance(data, Atoms):
            self.data = MolecularData.from_atoms(data)
        else:
            self.data = data
            
        # Convert to list if it's a single molecule
        if isinstance(self.data, dict):
            self.data = [self.data]
    
    def _get_js_data(self) -> str:
        """Convert molecular data to JavaScript-compatible JSON string."""
        return json.dumps(self.data)
    
    def show(self):
        """Display the viewer in a Jupyter notebook."""
        try:
            from IPython.display import display, HTML
            html = self._generate_html()
            display(HTML(html))
        except ImportError:
            print("IPython is not available. Use get_html() to get the HTML content.")
    
    def get_html(self) -> str:
        """Get the HTML content for the viewer."""
        return self._generate_html()
    
    def save_html(self, filename: str) -> None:
        """Save the viewer as an HTML file."""
        html = self.get_html()
        with open(filename, 'w') as f:
            f.write(html)
        print(f"Saved viewer to {filename}")
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the viewer."""
        # This will be implemented by subclasses
        raise NotImplementedError


class MolecularViewer(BaseViewer):
    """Molecular viewer with advanced settings and controls."""
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str], name: str = "molecule", **kwargs):
        super().__init__(data)
        self.name = name  # Base name for downloads
        self.settings = {
            "bondThreshold": 1.5,
            "bondThickness": 0.2,
            "atomSize": 1.0,
            "animationSpeed": 30,
            "forceScale": 1.0,
            "backgroundColor": "#1f2937",
            "style": "Cartoon",
            "showCell": True,
            "showBond": True,
            "showEnergyPlot": False,
            "showForces": False,
            "viewMode": "Perspective",
            "rotationMode": "TrackBall",
            "selectionMode": "Lasso",
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the molecular viewer."""
        js_data = self._get_js_data()
        
        # Try to load the static HTML file from the package
        try:
            template_path = pkg_resources.resource_filename('aseview', 'templates/molecular_viewer.html')
            with open(template_path, 'r') as f:
                html = f.read()
            
            # Inject download name and molecular data
            script = f"""
<script>
    document.addEventListener('DOMContentLoaded', function() {{
        if (typeof setDownloadName === 'function') {{
            setDownloadName('{self.name}');
        }}
        if (typeof setMolecularData === 'function') {{
            setMolecularData({js_data});
        }}
    }});
</script>
</body>
"""
            html = html.replace('</body>', script)
            return html
        except FileNotFoundError:
            # Fallback to inline HTML generation
            html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Molecular Viewer</title>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
    <script src="/static/js/styles.js"></script>
    <style>
        body {
            margin: 0;
            font-family: sans-serif;
            background-color: #111827;
            color: #f9fafb;
            overflow: hidden;
        }
        #container {
            width: 100vw;
            height: 100vh;
            position: relative;
        }
        #controls {
            position: absolute;
            top: 10px;
            left: 10px;
            z-index: 100;
            background: rgba(255, 255, 255, 0.8);
            padding: 10px;
            border-radius: 5px;
            color: black;
        }
        .control-group {
            margin-bottom: 10px;
        }
        .control-group label {
            display: block;
            margin-bottom: 5px;
        }
    </style>
</head>
<body>
    <div id="container"></div>
    <div id="controls">
        <div class="control-group">
            <label for="style-selector">Style</label>
            <select id="style-selector">
                <option value="default">Default</option>
                <option value="2d">2D</option>
                <option value="cartoon">Cartoon</option>
                <option value="neon">Neon</option>
                <option value="glossy">Glossy</option>
                <option value="metallic">Metallic</option>
                <option value="rowan">Rowan</option>
                <option value="grey">Grey</option>
            </select>
        </div>
        <div class="control-group">
            <label for="control-type">Controls</label>
            <select id="control-type">
                <option value="orbit">Orbit</option>
                <option value="trackball">Trackball</option>
            </select>
        </div>
        <div class="control-group">
            <label for="camera-mode">Camera</label>
            <select id="camera-mode">
                <option value="perspective">Perspective</option>
                <option value="orthographic">Orthographic</option>
            </select>
        </div>
        <div class="control-group">
            <label for="atom-scale">Atom Scale: <span id="atom-scale-value">0.5</span></label>
            <input type="range" id="atom-scale" min="0.1" max="2" step="0.1" value="0.5" />
        </div>
        <div class="control-group">
            <label for="bond-scale">Bond Scale: <span id="bond-scale-value">0.1</span></label>
            <input type="range" id="bond-scale" min="0.05" max="0.5" step="0.01" value="0.1" />
        </div>
        <div class="control-group">
            <label for="bond-cutoff-scale">Bond Cutoff Scale: <span id="bond-cutoff-scale-value">1.0</span></label>
            <input type="range" id="bond-cutoff-scale" min="0.5" max="1.5" step="0.05" value="1.0" />
        </div>
        <div class="control-group">
            <label for="force-scale">Force Vector Scale: <span id="force-scale-value">0.5</span></label>
            <input type="range" id="force-scale" min="0.1" max="2" step="0.1" value="0.5" />
        </div>
    </div>
    <script>
        const molecularData = {js_data};
        let scene, camera, renderer, controls, light;
        let molecule;

        let currentStyle = 'default';
        let controlType = 'orbit';
        let cameraMode = 'perspective';
        let atomScale = 0.5;
        let bondScale = 0.1;
        let bondCutoffScale = 1.0;
        let forceVectorScale = 0.5;

        function init() {
            scene = new THREE.Scene();
            renderer = new THREE.WebGLRenderer({ antialias: true });
            renderer.setSize(window.innerWidth, window.innerHeight);
            document.getElementById('container').appendChild(renderer.domElement);
            renderer.sortObjects = true;

            setupCamera();
            setupControls();

            light = new THREE.DirectionalLight(0xffffff, 0.8);
            scene.add(light);
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.2);
            scene.add(ambientLight);

            scene.background = new THREE.Color(0x111827);

            drawMolecule();
            animate();
        }

        function setupCamera() {
            const aspect = window.innerWidth / window.innerHeight;
            if (cameraMode === 'perspective') {
                camera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
            } else {
                const frustumSize = 10;
                camera = new THREE.OrthographicCamera(
                    (frustumSize * aspect) / -2,
                    (frustumSize * aspect) / 2,
                    frustumSize / 2,
                    frustumSize / -2,
                    0.1,
                    1000
                );
            }
            camera.position.z = 5;
        }

        function setupControls() {
            if (controls) {
                controls.dispose();
            }
            if (controlType === 'orbit') {
                controls = new THREE.OrbitControls(camera, renderer.domElement);
            } else {
                controls = new THREE.TrackballControls(camera, renderer.domElement);
                controls.rotateSpeed = 5.0;
                controls.zoomSpeed = 1.2;
                controls.panSpeed = 0.8;
            }
        }

        function drawMolecule() {
            if (molecule) {
                scene.remove(molecule);
            }
            molecule = new THREE.Group();

            const positions = molecularData.map((m) => m.positions.map((p) => new THREE.Vector3(...p)))[0];
            const symbols = molecularData[0].symbols;

            const uniqueSymbols = [...new Set(symbols)];
            const greyColorMap = {};
            for (let i = 0; i < uniqueSymbols.length; i++) {
                const lightness = (i + 1) / (uniqueSymbols.length + 1);
                greyColorMap[uniqueSymbols[i]] = new THREE.Color().setHSL(0, 0, lightness);
            }

            if (molecularData[0].forces) {
                const forceColor = 0xff0000;
                const headLength = 0.2;
                const headWidth = 0.1;
                for (let i = 0; i < positions.length; i++) {
                    const forceVec = new THREE.Vector3(...molecularData[0].forces[i]);
                    if (forceVec.lengthSq() === 0) continue;
                    const origin = positions[i];
                    const direction = forceVec.clone().normalize();
                    const length = forceVec.length() * forceVectorScale;
                    const arrowHelper = new THREE.ArrowHelper(direction, origin, length, forceColor, headLength, headWidth);
                    arrowHelper.renderOrder = 3;
                    molecule.add(arrowHelper);
                }
            }

            if (molecularData[0].cell) {
                const cellMaterial = new THREE.LineBasicMaterial({ color: 0x808080 });
                const cellVectors = molecularData[0].cell.map((c) => new THREE.Vector3(...c));
                const points = [
                    new THREE.Vector3(0, 0, 0), cellVectors[0],
                    new THREE.Vector3(0, 0, 0), cellVectors[1],
                    new THREE.Vector3(0, 0, 0), cellVectors[2],
                    cellVectors[0], cellVectors[0].clone().add(cellVectors[1]),
                    cellVectors[0], cellVectors[0].clone().add(cellVectors[2]),
                    cellVectors[1], cellVectors[1].clone().add(cellVectors[0]),
                    cellVectors[1], cellVectors[1].clone().add(cellVectors[2]),
                    cellVectors[2], cellVectors[2].clone().add(cellVectors[0]),
                    cellVectors[2], cellVectors[2].clone().add(cellVectors[1]),
                    cellVectors[0].clone().add(cellVectors[1]), cellVectors[0].clone().add(cellVectors[1]).add(cellVectors[2]),
                    cellVectors[0].clone().add(cellVectors[2]), cellVectors[0].clone().add(cellVectors[1]).add(cellVectors[2]),
                    cellVectors[1].clone().add(cellVectors[2]), cellVectors[0].clone().add(cellVectors[1]).add(cellVectors[2]),
                ];
                const geometry = new THREE.BufferGeometry().setFromPoints(points);
                const lineMesh = new THREE.LineSegments(geometry, cellMaterial);
                lineMesh.renderOrder = 2;
                molecule.add(lineMesh);
            }

            for (let i = 0; i < positions.length; i++) {
                for (let j = i + 1; j < positions.length; j++) {
                    const dist = positions[i].distanceTo(positions[j]);
                    if (dist < 1.6) {
                        let bond;
                        switch (currentStyle) {
                            case '2d':
                                bond = createBondStyle2D(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale);
                                break;
                            case 'cartoon':
                                bond = createBondStyleCartoon(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale);
                                break;
                            case 'neon':
                                bond = createBondStyleNeon(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale);
                                break;
                            case 'rowan':
                                bond = createBondStyleRowan(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale);
                                break;
                            case 'grey':
                                bond = createBondStyleGrey(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale, greyColorMap);
                                break;
                            default:
                                bond = createBondStyleDefault(positions[i], positions[j], symbols[i], symbols[j], bondScale, atomScale, currentStyle);
                                break;
                        }
                        if (bond) {
                            bond.renderOrder = 1;
                            molecule.add(bond);
                        }
                    }
                }
            }

            for (let i = 0; i < positions.length; i++) {
                let atom;
                switch (currentStyle) {
                    case '2d':
                        atom = createAtomStyle2D(positions[i], symbols[i], atomScale);
                        break;
                    case 'cartoon':
                        atom = createAtomStyleCartoon(positions[i], symbols[i], atomScale);
                        break;
                    case 'neon':
                        atom = createAtomStyleNeon(positions[i], symbols[i], atomScale);
                        break;
                    case 'glossy':
                        atom = createAtomStyleGlossy(positions[i], symbols[i], atomScale);
                        break;
                    case 'metallic':
                        atom = createAtomStyleMetallic(positions[i], symbols[i], atomScale);
                        break;
                    case 'rowan':
                        atom = createAtomStyleRowan(positions[i], symbols[i], atomScale);
                        break;
                    case 'grey':
                        atom = createAtomStyleGrey(positions[i], symbols[i], atomScale, greyColorMap[symbols[i]]);
                        break;
                    default:
                        atom = createAtomStyleDefault(positions[i], symbols[i], atomScale);
                        break;
                }
                atom.renderOrder = 0;
                molecule.add(atom);
            }
            scene.add(molecule);
        }

        function animate() {
            requestAnimationFrame(animate);
            controls.update();
            light.position.copy(camera.position);
            renderer.render(scene, camera);
        }

        document.getElementById('style-selector').addEventListener('change', (event) => {
            currentStyle = event.target.value;
            drawMolecule();
        });

        document.getElementById('control-type').addEventListener('change', (event) => {
            controlType = event.target.value;
            setupControls();
        });

        document.getElementById('camera-mode').addEventListener('change', (event) => {
            cameraMode = event.target.value;
            const oldPos = camera.position.clone();
            const oldZoom = camera.zoom;
            setupCamera();
            camera.position.copy(oldPos);
            camera.zoom = oldZoom;
            camera.updateProjectionMatrix();
            setupControls();
        });

        document.getElementById('atom-scale').addEventListener('input', (event) => {
            atomScale = parseFloat(event.target.value);
            document.getElementById('atom-scale-value').textContent = atomScale.toFixed(2);
            drawMolecule();
        });

        document.getElementById('bond-scale').addEventListener('input', (event) => {
            bondScale = parseFloat(event.target.value);
            document.getElementById('bond-scale-value').textContent = bondScale.toFixed(2);
            drawMolecule();
        });

        document.getElementById('bond-cutoff-scale').addEventListener('input', (event) => {
            bondCutoffScale = parseFloat(event.target.value);
            document.getElementById('bond-cutoff-scale-value').textContent = bondCutoffScale.toFixed(2);
            drawMolecule();
        });

        document.getElementById('force-scale').addEventListener('input', (event) => {
            forceVectorScale = parseFloat(event.target.value);
            document.getElementById('force-scale-value').textContent = forceVectorScale.toFixed(2);
            drawMolecule();
        });

        window.addEventListener('resize', () => {
            const aspect = window.innerWidth / window.innerHeight;
            if (camera.isPerspectiveCamera) {
                camera.aspect = aspect;
            } else {
                const frustumSize = 10;
                camera.left = (frustumSize * aspect) / -2;
                camera.right = (frustumSize * aspect) / 2;
                camera.top = frustumSize / 2;
                camera.bottom = frustumSize / -2;
            }
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
            if (controls.handleResize) controls.handleResize();
        }, false);

        init();
    </script>
</body>
</html>


class NormalViewer(BaseViewer):
    """Normal mode viewer for molecular vibrations."""
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str], **kwargs):
        super().__init__(data)
        self.settings = {
            "bondThreshold": 1.5,
            "bondThickness": 0.2,
            "atomSize": 1.0,
            "animationSpeed": 30,
            "displacementAmplitude": 1.0,
            "backgroundColor": "#1f2937",
            "style": "Cartoon",
            "showCell": True,
            "showBond": True,
            "showModeVector": False,
            "viewMode": "Perspective",
            "rotationMode": "TrackBall",
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the normal mode viewer."""
        js_data = self._get_js_data()
        
        # Try to load the static HTML file from the package
        try:
            template_path = pkg_resources.resource_filename('aseview', 'static/normal_viewer.html')
            with open(template_path, 'r') as f:
                html = f.read()
        except FileNotFoundError:
            # Fallback to inline HTML generation
            html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Normal Mode Viewer</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ font-family: Arial, sans-serif; background: #111827; color: #f9fafb; display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        .container {{ text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Normal Mode Viewer</h1>
        <p>Normal mode data loaded successfully</p>
        <div>
            <p>Number of configurations: {len(self.data)}</p>
            <p>Atoms: {len(self.data[0]['symbols']) if self.data else 0}</p>
        </div>
    </div>
    <script>
        function setMolecularData(data) {{
            console.log('Molecular data set:', data);
        }}
        setMolecularData({js_data});
    </script>
</body>
</html>
            """
            return html
        
        # Inject the molecular data into the HTML
        # We'll add a script tag at the end of the body to initialize the viewer
        script = f"""
<script>
    // Set the molecular data
    document.addEventListener('DOMContentLoaded', function() {{
        if (typeof setMolecularData === 'function') {{
            setMolecularData({js_data});
        }}
    }});
</script>
</body>
"""
        
        # Replace the closing body tag with our script and closing body tag
        html = html.replace('</body>', script)
        return html


class OverlayViewer(BaseViewer):
    """Overlay viewer for comparing multiple molecules."""
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str], **kwargs):
        super().__init__(data)
        self.settings = {
            "bondThreshold": 1.5,
            "bondThickness": 0.2,
            "atomSize": 1.0,
            "backgroundColor": "#1f2937",
            "style": "Cartoon",
            "colorBy": "Atom",
            "showCell": True,
            "showBond": True,
            "viewMode": "Perspective",
            "rotationMode": "TrackBall",
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the overlay viewer."""
        js_data = self._get_js_data()
        
        # Try to load the static HTML file from the package
        try:
            template_path = pkg_resources.resource_filename('aseview', 'static/overlay_viewer.html')
            with open(template_path, 'r') as f:
                html = f.read()
        except FileNotFoundError:
            # Fallback to inline HTML generation
            html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Overlay Viewer</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ font-family: Arial, sans-serif; background: #111827; color: #f9fafb; display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        .container {{ text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Overlay Viewer</h1>
        <p>Overlay data loaded successfully</p>
        <div>
            <p>Number of molecules: {len(self.data)}</p>
            <p>Atoms in first molecule: {len(self.data[0]['symbols']) if self.data else 0}</p>
        </div>
    </div>
    <script>
        function setMolecularData(data) {{
            console.log('Molecular data set:', data);
        }}
        setMolecularData({js_data});
    </script>
</body>
</html>
            """
            return html
        
        # Inject the molecular data into the HTML
        # We'll add a script tag at the end of the body to initialize the viewer
        script = f"""
<script>
    // Set the molecular data
    document.addEventListener('DOMContentLoaded', function() {{
        if (typeof setMolecularData === 'function') {{
            setMolecularData({js_data});
        }}
    }});
</script>
</body>
"""
        
        # Replace the closing body tag with our script and closing body tag
        html = html.replace('</body>', script)
        return html
