import json
from IPython.display import HTML
from typing import List, Union, Optional, Dict
import numpy as np
from scipy.spatial.distance import pdist, squareform
import ase
from ase import Atoms
from ase.build import molecule
from ase.io.trajectory import TrajectoryReader
from ase.calculators.calculator import PropertyNotImplementedError
from rdkit import Chem
from rdkit.Chem import AllChem

atomic_symbols2hex = {
    'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00', 'B': '#FFB5B5',
    'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'F': '#90E050', 'Ne': '#B3E3F5',
    'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000',
    'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00',
    'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7', 'Mn': '#9C7AC7',
    'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050', 'Cu': '#C88033', 'Zn': '#7D80B0',
    'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929',
    'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0',
    'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C',
    'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080',
    'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F',
    'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7',
    'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7',
    'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38',
    'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB',
    'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0',
    'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45',
    'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF',
    'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2',
    'Cm': '#785CE3', 'default': '#FFC0CB' # Pink as default
}

covalent_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71,
    'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.70,
    'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24,
    'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20,
    'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42,
    'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 'Cs': 2.44, 'Ba': 2.15,
    'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98,
    'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
    'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41,
    'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
    'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00,
    'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69, 'default': 0.7
}



def viewer(
    atoms_or_traj: Union[Atoms, List[Atoms], TrajectoryReader],
    rotation_mode: str = 'orbit',
    bg_color: str = '#2D3748',
    show_bonds: bool = True,
    show_forces: bool = False,
    show_cell: bool = True,
    show_energy: bool = True,
    atom_scale: float = 1.0,
    bond_cutoff_scale: float = 1.2,
    bond_thickness: float = 0.08,
    force_scale: float = 0.3,
    animation_fps: int = 10,
    width: int = 800,
    height: int = 600,
    write_html: Optional[str] = None,
    ) -> HTML:
    """
    Visualizes ASE Atoms objects or TrajectoryReader in 3D within a Jupyter Notebook.
    (Full version including GIF saving functionality)
    """
    if isinstance(atoms_or_traj, Atoms):
        frames = [atoms_or_traj]
    else:
        frames = [atoms for atoms in atoms_or_traj]

    trajectory_data = []
    has_energy = False
    for atoms in frames:
        frame_data = {
            "symbols": atoms.get_chemical_symbols(),
            "positions": atoms.get_positions().tolist(),
            "cell": atoms.get_cell().tolist(),
        }
        try:
            forces = atoms.get_forces(apply_constraint=False)
            frame_data["forces"] = forces.tolist()
        except (RuntimeError, PropertyNotImplementedError): pass
        try:
            energy = atoms.get_potential_energy()
            frame_data["energy"] = energy
            has_energy = True
        except (RuntimeError, PropertyNotImplementedError): frame_data["energy"] = None
        trajectory_data.append(frame_data)

    json_data = json.dumps(trajectory_data)
    effective_show_energy = show_energy and has_energy

    atom_info_parts = []
    symbols = sorted([s for s in covalent_radii if s != 'default'])

    for symbol in symbols:
        radius = covalent_radii.get(symbol)
        hex_color = atomic_symbols2hex.get(symbol, atomic_symbols2hex['default'])
        color_int_str = f"0x{hex_color.lstrip('#')}"
        atom_info_parts.append(f"'{symbol}':{{color:{color_int_str},radius:{radius}}}")

    default_radius = covalent_radii['default']
    default_hex_color = atomic_symbols2hex['default']
    default_color_int_str = f"0x{default_hex_color.lstrip('#')}"
    atom_info_parts.append(f"'default':{{color:{default_color_int_str},radius:{default_radius}}}")

    atom_info_js = ", ".join(atom_info_parts)

    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ASE 3D Viewer</title>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css">
        <style>
            body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; overflow: hidden; user-select: none; }}
            .main-container {{ display: flex; flex-direction: column; width: {width}px; background-color: #1A202C; border-radius: 8px; border: 1px solid #4A5568; }}
            .viewer-wrapper {{ position: relative; width: 100%; height: {height}px; overflow: hidden; }}
            .viewer-container {{ width: 100%; height: 100%; display: flex; }}
            .sidebar {{
                width: 240px; min-width: 240px; padding: 20px;
                background-color: #1A202C; color: #E2E8F0;
                display: flex; flex-direction: column; gap: 20px;
                overflow-y: auto; overflow-x: hidden;
                transition: width 0.3s ease, min-width 0.3s ease, padding 0.3s ease;
                flex-shrink: 0; border-right: 1px solid #4A5568;
                resize: none;
            }}
            .viewer-wrapper.sidebar-collapsed .sidebar {{ width: 0; min-width: 0; padding: 20px 0; border-right-color: transparent; }}
            #sidebar-toggle-btn {{
                position: absolute; top: 10px; left: 10px; z-index: 1000;
                background-color: #2D3748; color: white; border: 1px solid #4A5568;
                border-radius: 6px; width: 30px; height: 30px;
                display: flex; align-items: center; justify-content: center;
                cursor: pointer; font-size: 16px; font-weight: bold;
                font-family: "Courier New", Courier, monospace;
                transition: background-color 0.2s;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2);
            }}
            #sidebar-toggle-btn:hover {{ background-color: #4A5568; }}
            
            #top-right-controls {{
                position: absolute;
                top: 15px;
                right: 15px;
                z-index: 100;
            }}

            #save-gif-btn {{
                background-color: rgba(74, 85, 104, 0.7);
                border: 1px solid rgba(255, 255, 255, 0.1);
                border-radius: 8px;
                color: #E2E8F0;
                padding: 10px 16px;
                cursor: pointer;
                font-size: 14px;
                font-weight: 600;
                display: flex;
                align-items: center;
                gap: 8px;
                transition: all 0.3s ease;
                backdrop-filter: blur(5px);
            }}
            #save-gif-btn:hover {{
                background-color: rgba(99, 110, 130, 0.85);
                transform: translateY(-2px);
            }}
            #save-gif-btn.saving {{
                background: #6B7280;
                cursor: not-allowed;
                transform: none;
            }}
            #save-gif-btn .fa-spinner {{
                animation: fa-spin 1s linear infinite;
            }}
            @keyframes fa-spin {{
                0% {{ transform: rotate(0deg); }}
                100% {{ transform: rotate(360deg); }}
            }}


            .resize-handle {{ position: absolute; bottom: 0; left: 0; right: 0; height: 4px; background-color: #4A5568; cursor: ns-resize; transition: background-color 0.2s; }}
            .resize-handle:hover {{ background-color: #63B3ED; }}
            .resize-handle::before {{ content: ''; position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 30px; height: 2px; background-color: #A0AEC0; border-radius: 1px; }}
            .sidebar h2 {{ font-size: 18px; margin: 0 0 10px 0; padding-bottom: 10px; border-bottom: 1px solid #4A5568; color: #A0AEC0; white-space: nowrap; }}
            details {{ color: #E2E8F0; border-bottom: 1px solid #2d3748; white-space: nowrap; }}
            details[open] {{ padding-bottom: 15px; }}
            summary {{ cursor: pointer; font-size: 1rem; font-weight: bold; padding: 10px 0; list-style-position: inside; }}
            summary::marker {{ color: #A0AEC0; }}
            .details-content {{ padding-top: 10px; display: flex; flex-direction: column; gap: 12px; }}
            .toggle-switch {{ display: flex; justify-content: space-between; align-items: center; width: 100%; font-size: 14px; }}
            .toggle-switch input {{ display: none; }}
            .toggle-switch .slider {{ position: relative; cursor: pointer; width: 40px; height: 22px; background-color: #ccc; border-radius: 22px; transition: background-color 0.2s; flex-shrink: 0; }}
            .toggle-switch .slider:before {{ position: absolute; content: ""; height: 16px; width: 16px; left: 3px; bottom: 3px; background-color: white; border-radius: 50%; transition: transform 0.2s; }}
            .toggle-switch input:checked + .slider {{ background-color: #4CAF50; }}
            .toggle-switch input:checked + .slider:before {{ transform: translateX(18px); }}
            select {{ background-color: #2D3748; color: #E2E8F0; border: 1px solid #4A5568; border-radius: 4px; padding: 4px; width: 100%; }}
            input[type="color"] {{ width: 100%; height: 25px; border: 1px solid #4A5568; border-radius: 4px; padding: 2px; background-color: #2D3748;}}
            .animation-controls {{ display: flex; align-items: center; gap: 10px; }}
            .control-btn {{ background-color: #4A5568; border: none; border-radius: 50%; width: 32px; height: 32px; cursor: pointer; display: flex; justify-content: center; align-items: center; flex-shrink: 0; transition: background-color 0.2s; color: white; }}
            .control-btn:hover {{ background-color: #636e82; }}
            .control-btn.saving {{ background: #6B7280; cursor: not-allowed; }}
            .control-btn svg {{ fill: white; }}
            input[type="range"] {{ width: 100%; }}
            .slider-group {{ display: flex; flex-direction: column; gap: 4px;}}
            .slider-label, .slider-value {{ font-size: 12px; color: #A0AEC0; }}
            #canvas-container {{ flex-grow: 1; height: 100%; cursor: grab; background-color: {bg_color}; min-width: 0; position: relative; }}
            #plot-container {{ width: 100%; padding: 10px; box-sizing: border-box; background-color: #1A202C; display: { 'block' if effective_show_energy else 'none' }; }}
            .resizing {{ user-select: none; }}
            .resizing * {{ pointer-events: none; }}
        </style>
    </head>
    <body>
        <div class="main-container">
            <div class="viewer-wrapper" id="viewer-wrapper">
                <div class="viewer-container">
                    <div class="sidebar" id="sidebar">
                        <h2 style="text-align: right;">Molecule Viewer</h2>
                        <details open><summary>Display Options</summary><div class="details-content">
                            <label class="toggle-switch"><span>Bonds</span><input type="checkbox" id="bonds-toggle"><span class="slider"></span></label>
                            <label class="toggle-switch"><span>Forces</span><input type="checkbox" id="forces-toggle"><span class="slider"></span></label>
                            <label class="toggle-switch"><span>Unit Cell</span><input type="checkbox" id="cell-toggle"><span class="slider"></span></label>
                            <label class="toggle-switch"><span>Shadows</span><input type="checkbox" id="shadows-toggle"><span class="slider"></span></label>
                            <label class="toggle-switch" id="energy-toggle-label" style="display: {'flex' if has_energy else 'none'};"><span>Energy Plot</span><input type="checkbox" id="energy-plot-toggle"><span class="slider"></span></label>
                        </div></details>
                        <details open><summary>Advanced Settings</summary><div class="details-content">
                            <div class="slider-group"><label for="rotation-mode-select" class="slider-label">Rotation Mode</label><select id="rotation-mode-select"><option value="orbit">Orbit</option><option value="trackball">Trackball</option></select></div>
                            <div class="slider-group"><label for="style-select" class="slider-label">Style</label><select id="style-select"><option value="matte">Matte</option><option value="glossy">Glossy</option><option value="metallic">Metallic</option><option value="2d">2D Structure</option><option value="toon">Toon</option><option value="neon">Neon</option></select></div>
                            <div class="slider-group"><label for="atom-scale-slider" class="slider-label">Atom Size</label><input type="range" id="atom-scale-slider" min="0.1" max="2.0" step="0.05" value="{atom_scale}"><span id="atom-scale-label" class="slider-value">Scale: {atom_scale:.2f}</span></div>
                            <div class="slider-group"><label for="bond-thickness-slider" class="slider-label">Bond Thickness</label><input type="range" id="bond-thickness-slider" min="0.02" max="0.2" step="0.01" value="{bond_thickness}"><span id="bond-thickness-label" class="slider-value">Value: {bond_thickness:.2f}</span></div>
                            <div class="slider-group"><label for="bond-cutoff-slider" class="slider-label">Bond Cutoff Scale</label><input type="range" id="bond-cutoff-slider" min="0.8" max="1.5" step="0.01" value="{bond_cutoff_scale}"><span id="bond-cutoff-label" class="slider-value">Scale: {bond_cutoff_scale:.2f}</span></div>
                            <div class="slider-group"><label for="force-scale-slider" class="slider-label">Force Vector Scale</label><input type="range" id="force-scale-slider" min="0.1" max="1.0" step="0.05" value="{force_scale}"><span id="force-scale-label" class="slider-value">Scale: {force_scale:.2f}</span></div>
                            <div class="slider-group"><label for="animation-speed-slider" class="slider-label">Animation Speed (FPS)</label><input type="range" id="animation-speed-slider" min="1" max="180" step="1" value="{animation_fps}"><span id="animation-speed-label" class="slider-value">FPS: {animation_fps}</span></div>
                            <div class="slider-group"><label for="bg-color-picker" class="slider-label">Background Color</label><input type="color" id="bg-color-picker" value="{bg_color}"></div>
                        </div></details>
                        <div id="animation-section" style="padding-top: 10px; display: none;"><h3 style="margin:0 0 10px 0;">Animation</h3><div class="animation-controls">
                            <button id="play-pause-btn" class="control-btn" title="Play/Pause"><svg id="play-icon" width="14" height="14" viewBox="0 0 24 24"><path d="M8 5v14l11-7z"/></svg><svg id="pause-icon" width="14" height="14" viewBox="0 0 24 24" style="display:none;"><path d="M6 19h4V5H6v14zm8-14v14h4V5h-4z"/></svg></button>
                            <input type="range" id="frame-slider" min="0" value="0">
                            <button id="copy-xyz-btn" class="control-btn" title="Copy XYZ to Clipboard"><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path></svg></button>
                        </div><span id="frame-label" class="slider-value">Frame: 0 / 0</span></div>
                        <div class="resize-handle" id="resize-handle"></div>
                    </div>
                    <div id="canvas-container"></div>
                </div>
                <div id="sidebar-toggle-btn">⟨</div>
                <div id="top-right-controls">
                     <button id="save-gif-btn" title="Save as GIF">
                        <i class="fas fa-file-download"></i>
                        <span>Save as GIF</span>
                    </button>
                </div>
            </div>
            <div id="plot-container">
                <canvas id="energy-plot"></canvas>
            </div>
        </div>

        <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/chart.js/2.9.4/Chart.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/gif.js/0.2.0/gif.js"></script>

        <script>
            // Data and Initial State
            const trajectoryData = {json_data};
            const hasEnergy = { 'true' if has_energy else 'false' };
            const initialShowBonds = {str(show_bonds).lower()};
            const initialShowForces = {str(show_forces).lower()};
            const initialShowCell = {str(show_cell).lower()};
            const initialShowEnergy = {str(effective_show_energy).lower()};
            const atomInfo = {{{atom_info_js}}};
            
            // Global Variables
            let scene, camera, renderer, controls, energyPlot, ambientLight, directionalLight;
            let atomGroup, bondGroup, forceGroup, cellGroup;
            let animationInterval = null, isPlaying = false;
            let currentFrameIndex = 0;
            let isResizing = false;

            // UI Elements
            const viewerWrapper = document.getElementById('viewer-wrapper');
            const sidebar = document.getElementById('sidebar');
            const sidebarToggleBtn = document.getElementById('sidebar-toggle-btn');
            const resizeHandle = document.getElementById('resize-handle');
            const rotationModeSelect = document.getElementById('rotation-mode-select');
            const styleSelect = document.getElementById('style-select');
            const bgColorPicker = document.getElementById('bg-color-picker');
            const atomScaleSlider = document.getElementById('atom-scale-slider'), atomScaleLabel = document.getElementById('atom-scale-label');
            const bondCutoffSlider = document.getElementById('bond-cutoff-slider'), bondCutoffLabel = document.getElementById('bond-cutoff-label');
            const bondThicknessSlider = document.getElementById('bond-thickness-slider'), bondThicknessLabel = document.getElementById('bond-thickness-label');
            const forceScaleSlider = document.getElementById('force-scale-slider'), forceScaleLabel = document.getElementById('force-scale-label');
            const animationSpeedSlider = document.getElementById('animation-speed-slider'), animationSpeedLabel = document.getElementById('animation-speed-label');
            const animationSection = document.getElementById('animation-section'), slider = document.getElementById('frame-slider'), frameLabel = document.getElementById('frame-label');
            const playPauseBtn = document.getElementById('play-pause-btn'), playIcon = document.getElementById('play-icon'), pauseIcon = document.getElementById('pause-icon');
            const copyXyzBtn = document.getElementById('copy-xyz-btn');
            const saveGifBtn = document.getElementById('save-gif-btn');
            const bondsToggle = document.getElementById('bonds-toggle'), forcesToggle = document.getElementById('forces-toggle'), cellToggle = document.getElementById('cell-toggle'), shadowsToggle = document.getElementById('shadows-toggle'), energyPlotToggle = document.getElementById('energy-plot-toggle');

            // Config Variables
            let shadowsEnabled = false;
            let currentStyle = 'matte';
            let currentAtomScale = {atom_scale};
            let currentBondCutoffFactor = {bond_cutoff_scale};
            let currentBondThickness = {bond_thickness};
            let currentForceScale = {force_scale};
            let currentAnimationSpeed = 1000 / {animation_fps};

            // Main execution starts here
            init();

            function init() {{
                const container = document.getElementById('canvas-container');
                scene = new THREE.Scene();
                scene.background = new THREE.Color(bgColorPicker.value);
                camera = new THREE.PerspectiveCamera(75, container.clientWidth / container.clientHeight, 0.1, 1000);
                camera.position.z = 20;
                scene.add(camera);

                renderer = new THREE.WebGLRenderer({{ antialias: true, preserveDrawingBuffer: true }});
                renderer.setSize(container.clientWidth, container.clientHeight);
                renderer.setPixelRatio(window.devicePixelRatio);
                renderer.shadowMap.enabled = shadowsEnabled;
                container.appendChild(renderer.domElement);

                ambientLight = new THREE.AmbientLight(0xcccccc, 0.8);
                scene.add(ambientLight);
                directionalLight = new THREE.DirectionalLight(0xffffff, 0.6);
                directionalLight.position.set(1, 1, 2);
                directionalLight.castShadow = true;
                camera.add(directionalLight);

                atomGroup = new THREE.Group(); bondGroup = new THREE.Group();
                forceGroup = new THREE.Group(); cellGroup = new THREE.Group();
                scene.add(atomGroup, bondGroup, forceGroup, cellGroup);

                setupControls(rotationModeSelect.value);
                setupUI();
                setupResizeHandle();
                if (hasEnergy) setupEnergyPlot();
                
                window.addEventListener('resize', onWindowResize);
                updateScene(0);
                animate();
            }}

            function animate() {{
                requestAnimationFrame(animate);
                controls.update();
                renderer.render(scene, camera);
            }}

            function onWindowResize() {{
                const container = document.getElementById('canvas-container');
                if (!container || container.clientWidth === 0 || container.clientHeight === 0) return;
                camera.aspect = container.clientWidth / container.clientHeight;
                camera.updateProjectionMatrix();
                renderer.setSize(container.clientWidth, container.clientHeight);
                if(energyPlot) energyPlot.resize();
            }}

            function togglePlay() {{
                isPlaying = !isPlaying;
                if (isPlaying) {{
                    playIcon.style.display = 'none';
                    pauseIcon.style.display = 'block';
                    const speed = 1000 / parseInt(document.getElementById('animation-speed-slider').value);
                    animationInterval = setInterval(() => {{
                        let nextFrame = (parseInt(slider.value) + 1) % trajectoryData.length;
                        updateScene(nextFrame);
                    }}, speed);
                }} else {{
                    playIcon.style.display = 'block';
                    pauseIcon.style.display = 'none';
                    clearInterval(animationInterval);
                }}
            }}

            function setupControls(mode) {{
                if (controls) controls.dispose();
                const container = document.getElementById('canvas-container');
                if (mode === 'trackball') {{
                    controls = new THREE.TrackballControls(camera, container);
                }} else {{
                    controls = new THREE.OrbitControls(camera, container);
                    controls.enableDamping = true;
                }}
            }}

            function setupUI() {{
                sidebarToggleBtn.addEventListener('click', () => {{
                    const isCollapsed = viewerWrapper.classList.toggle('sidebar-collapsed');
                    sidebarToggleBtn.innerHTML = isCollapsed ? '⟩' : '⟨';
                    setTimeout(onWindowResize, 310);
                }});

                rotationModeSelect.value = "{rotation_mode}";
                rotationModeSelect.addEventListener('change', (e) => setupControls(e.target.value));
                bgColorPicker.addEventListener('input', (e) => {{ scene.background.set(e.target.value); }});

                // Only show animation controls if there is a trajectory
                if (trajectoryData.length > 1) {{
                    animationSection.style.display = 'block';
                    slider.max = trajectoryData.length - 1;
                    slider.addEventListener('input', (e) => updateScene(parseInt(e.target.value)));
                    playPauseBtn.addEventListener('click', togglePlay);
                    copyXyzBtn.addEventListener('click', copyXyzToClipboard);
                }}
                
                // The GIF save button should always be listening
                saveGifBtn.addEventListener('click', saveAsGif);

                bondsToggle.checked = initialShowBonds; bondGroup.visible = initialShowBonds;
                forcesToggle.checked = initialShowForces; forceGroup.visible = initialShowForces;
                cellToggle.checked = initialShowCell; cellGroup.visible = initialShowCell;
                shadowsToggle.checked = shadowsEnabled;

                bondsToggle.addEventListener('change', (e) => {{ bondGroup.visible = e.target.checked; }});
                forcesToggle.addEventListener('change', (e) => {{ forceGroup.visible = e.target.checked; }});
                cellToggle.addEventListener('change', (e) => {{ cellGroup.visible = e.target.checked; }});

                shadowsToggle.addEventListener('change', (e) => {{
                    shadowsEnabled = e.target.checked;
                    updateScene(currentFrameIndex, false);
                }});

                if (hasEnergy) {{
                    energyPlotToggle.checked = initialShowEnergy;
                    document.getElementById('plot-container').style.display = initialShowEnergy ? 'block' : 'none';
                    energyPlotToggle.addEventListener('change', (e) => {{
                        document.getElementById('plot-container').style.display = e.target.checked ? 'block' : 'none';
                    }});
                }}

                atomScaleSlider.addEventListener('input', (e) => {{ currentAtomScale = parseFloat(e.target.value); atomScaleLabel.textContent = `Scale: ${{currentAtomScale.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                bondCutoffSlider.addEventListener('input', (e) => {{ currentBondCutoffFactor = parseFloat(e.target.value); bondCutoffLabel.textContent = `Scale: ${{currentBondCutoffFactor.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                bondThicknessSlider.addEventListener('input', (e) => {{ currentBondThickness = parseFloat(e.target.value); bondThicknessLabel.textContent = `Value: ${{currentBondThickness.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                forceScaleSlider.addEventListener('input', (e) => {{ currentForceScale = parseFloat(e.target.value); forceScaleLabel.textContent = `Scale: ${{currentForceScale.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                
                animationSpeedSlider.addEventListener('input', (e) => {{
                    const fps = parseInt(e.target.value);
                    animationSpeedLabel.textContent = `FPS: ${{fps}}`;
                    currentAnimationSpeed = 1000 / fps;
                    if (isPlaying) {{
                        clearInterval(animationInterval);
                        animationInterval = setInterval(() => {{ let nextFrame = (parseInt(slider.value) + 1) % trajectoryData.length; updateScene(nextFrame); }}, currentAnimationSpeed);
                    }}
                }});
                styleSelect.addEventListener('change', (e) => {{ currentStyle = e.target.value; updateScene(currentFrameIndex, false); }});
            }}

            function updateScene(frameIdx, updateSlider=true) {{
                currentFrameIndex = frameIdx;
                clearGroup(atomGroup); clearGroup(bondGroup); clearGroup(forceGroup); clearGroup(cellGroup);
                const frameData = trajectoryData[currentFrameIndex];
                if (!frameData) return;

                renderer.shadowMap.enabled = shadowsEnabled;
                const isLightNeeded = (currentStyle !== '2d' && currentStyle !== 'neon');
                ambientLight.visible = isLightNeeded;
                directionalLight.visible = isLightNeeded;

                const positions = frameData.positions.map(p => new THREE.Vector3(...p));
                const symbols = frameData.symbols;

                positions.forEach((pos, i) => {{ const atom = createAtom(pos, symbols[i]); if (atom) atomGroup.add(atom); }});

                if (bondsToggle.checked) {{
                    for (let i = 0; i < positions.length; i++) {{
                        for (let j = i + 1; j < positions.length; j++) {{
                            const r_i = (atomInfo[symbols[i]] || atomInfo['default']).radius;
                            const r_j = (atomInfo[symbols[j]] || atomInfo['default']).radius;
                            const cutoff = (r_i + r_j) * currentBondCutoffFactor;
                            if (positions[i].distanceTo(positions[j]) < cutoff) createBond(positions[i], positions[j], symbols[i], symbols[j]);
                        }}
                    }}
                }}

                if (forcesToggle.checked && frameData.forces) {{
                    frameData.forces.forEach((forceArr, i) => {{
                        const forceVec = new THREE.Vector3(...forceArr);
                        const length = forceVec.length() * currentForceScale;
                        if (length > 0.01) {{
                            const arrow = new THREE.ArrowHelper(forceVec.normalize(), positions[i], length, 0xE53E3E, length * 0.25, length * 0.2);
                            arrow.line.material.depthTest = false; arrow.cone.material.depthTest = false;
                            forceGroup.add(arrow);
                        }}
                    }});
                }}

                if (cellToggle.checked && frameData.cell && frameData.cell.flat().some(v => v !== 0)) drawCell(frameData.cell);

                if (updateSlider && trajectoryData.length > 1) {{
                    slider.value = currentFrameIndex;
                    frameLabel.textContent = `Frame: ${{currentFrameIndex}} / ${{trajectoryData.length - 1}}`;
                }}
                if (hasEnergy) updatePlotHighlight(currentFrameIndex);
            }}

            function createAtom(pos, symbol) {{
                const info = atomInfo[symbol] || atomInfo['default'];
                const scaledRadius = info.radius * currentAtomScale;
                let geometry, material, sphere;
                switch (currentStyle) {{
                    case '2d':
                        const canvas = document.createElement('canvas'); const context = canvas.getContext('2d');
                        canvas.width = 128; canvas.height = 128;
                        context.beginPath(); context.arc(64, 64, 60, 0, 2 * Math.PI, false); context.fillStyle = 'black'; context.fill();
                        context.beginPath(); context.arc(64, 64, 56, 0, 2 * Math.PI, false); context.fillStyle = 'white'; context.fill();
                        context.font = 'bold 50px Arial'; context.fillStyle = 'black'; context.textAlign = 'center'; context.textBaseline = 'middle';
                        context.fillText(symbol, 64, 64);
                        const texture = new THREE.CanvasTexture(canvas);
                        const spriteMaterial = new THREE.SpriteMaterial({{ map: texture, transparent: true, depthTest: false }});
                        const sprite = new THREE.Sprite(spriteMaterial);
                        sprite.position.copy(pos);
                        sprite.scale.set(scaledRadius * 2.5, scaledRadius * 2.5, 1);
                        return sprite;
                    case 'toon':
                        geometry = new THREE.SphereGeometry(scaledRadius, 32, 32);
                        const toonMaterial = new THREE.MeshToonMaterial({{ color: new THREE.Color(info.color), transparent: true, opacity: info.opacity || 1.0 }});
                        sphere = new THREE.Mesh(geometry, toonMaterial);
                        sphere.position.copy(pos);
                        const outlineGeometry = new THREE.SphereGeometry(scaledRadius * 1.1, 32, 32);
                        const outlineMaterial = new THREE.MeshBasicMaterial({{ color: 0x000000, transparent: true, opacity: info.opacity || 1.0, side: THREE.BackSide }});
                        const outline = new THREE.Mesh(outlineGeometry, outlineMaterial);
                        outline.position.copy(pos);
                        const toonGroup = new THREE.Group();
                        toonGroup.add(outline); toonGroup.add(sphere);
                        if (shadowsEnabled) {{ sphere.castShadow = true; sphere.receiveShadow = true; outline.castShadow = true; }}
                        return toonGroup;
                    case 'neon':
                        geometry = new THREE.SphereGeometry(scaledRadius, 32, 32); material = new THREE.MeshBasicMaterial({{ color: new THREE.Color(info.color).multiplyScalar(1.5) }});
                        sphere = new THREE.Mesh(geometry, material); sphere.position.copy(pos);
                        return sphere;
                    default:
                        geometry = new THREE.SphereGeometry(scaledRadius, 32, 32); material = getStandardMaterial(info.color, currentStyle);
                        sphere = new THREE.Mesh(geometry, material); sphere.position.copy(pos);
                        if (shadowsEnabled) {{ sphere.castShadow = true; sphere.receiveShadow = true; }}
                        return sphere;
                }}
            }}

            function createBond(p1, p2, sym1, sym2) {{
                if (currentStyle === '2d') {{
                    const dir = p2.clone().sub(p1).normalize();
                    const offset1 = (atomInfo[sym1] || atomInfo.default).radius*currentAtomScale*0.4; const offset2 = (atomInfo[sym2] || atomInfo.default).radius*currentAtomScale*0.4;
                    const startPos = p1.clone().add(dir.clone().multiplyScalar(offset1)); const endPos = p2.clone().sub(dir.clone().multiplyScalar(offset2));
                    if (startPos.distanceTo(endPos) <= 0) return;
                    const geometry = new THREE.CylinderGeometry(currentBondThickness, currentBondThickness, startPos.distanceTo(endPos), 8, 1);
                    const material = new THREE.MeshBasicMaterial({{ color: 0x666666 }}); const bond = new THREE.Mesh(geometry, material);
                    bond.position.copy(startPos).add(endPos).multiplyScalar(0.5); bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
                    bondGroup.add(bond);
                    return;
                }}
                if (currentStyle === 'toon') {{
                    const path = new THREE.Line3(p1, p2); const distance = path.distance();
                    if (distance <= 0) return;
                    const geometry = new THREE.CylinderGeometry(currentBondThickness, currentBondThickness, distance, 8);
                    const material = new THREE.MeshToonMaterial({{ color: 0x000000 }}); const bond = new THREE.Mesh(geometry, material);
                    bond.position.lerpVectors(p1, p2, 0.5); bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), path.delta(new THREE.Vector3()).normalize());
                    if (shadowsEnabled) {{ bond.castShadow = true; }}
                    bondGroup.add(bond);
                    return;
                }}
                const midPoint = p1.clone().add(p2).multiplyScalar(0.5);
                const bond1 = createHalfBond(p1, midPoint, (atomInfo[sym1] || atomInfo.default).color); const bond2 = createHalfBond(midPoint, p2, (atomInfo[sym2] || atomInfo.default).color);
                if (bond1) bondGroup.add(bond1); if (bond2) bondGroup.add(bond2);
            }}

            function createHalfBond(start, end, color) {{
                if (start.distanceTo(end) <= 0) return null;
                const path = new THREE.LineCurve3(start, end); const geometry = new THREE.TubeGeometry(path, 1, currentBondThickness, 8, false);
                const material = (currentStyle === 'neon') ? new THREE.MeshBasicMaterial({{ color: 0xffffff }}) : getStandardMaterial(color, currentStyle);
                const bond = new THREE.Mesh(geometry, material);
                if (shadowsEnabled) {{ bond.castShadow = true; }}
                return bond;
            }}

            function getStandardMaterial(color, type) {{
                switch(type) {{
                    case 'glossy': return new THREE.MeshPhongMaterial({{ color: color, shininess: 100, specular: 0x222222 }});
                    case 'metallic': return new THREE.MeshStandardMaterial({{ color: color, metalness: 0.3, roughness: 0.4 }});
                    default: return new THREE.MeshLambertMaterial({{ color: color }});
                }}
            }}
            
            function clearGroup(group) {{
                while (group.children.length > 0) {{
                    const child = group.children[0]; group.remove(child);
                    if(child.geometry) child.geometry.dispose();
                    if(child.material) {{
                        if (Array.isArray(child.material)) {{ child.material.forEach(m => m.dispose()); }}
                        else if (child.material.map) {{ child.material.map.dispose(); }}
                        if (child.material.dispose) {{ child.material.dispose(); }}
                    }}
                    if(child.children && child.children.length > 0) clearGroup(child);
                }}
            }}

            function drawCell(cell) {{
                const [v1, v2, v3] = cell.map(c => new THREE.Vector3(...c));
                const points = [ new THREE.Vector3(0,0,0), v1, v2, v3, v1.clone().add(v2), v1.clone().add(v3), v2.clone().add(v3), v1.clone().add(v2).add(v3) ];
                const edges = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]];
                const lineMaterial = new THREE.LineBasicMaterial({{ color: 0x666666 }});
                edges.forEach(edge => {{
                    const geometry = new THREE.BufferGeometry().setFromPoints([points[edge[0]], points[edge[1]]]);
                    cellGroup.add(new THREE.Line(geometry, lineMaterial));
                }});
            }}

            function setupResizeHandle() {{
                let startY = 0, startHeight = 0;
                resizeHandle.addEventListener('mousedown', (e) => {{
                    e.preventDefault(); isResizing = true;
                    startY = e.clientY;
                    startHeight = viewerWrapper.offsetHeight;
                    document.body.classList.add('resizing');
                    document.addEventListener('mousemove', handleMouseMove);
                    document.addEventListener('mouseup', handleMouseUp);
                }});
                function handleMouseMove(e) {{
                    if (!isResizing) return;
                    const deltaY = e.clientY - startY;
                    const newHeight = Math.max(300, Math.min(800, startHeight + deltaY));
                    viewerWrapper.style.height = newHeight + 'px';
                    onWindowResize();
                }}
                function handleMouseUp() {{
                    isResizing = false;
                    document.body.classList.remove('resizing');
                    document.removeEventListener('mousemove', handleMouseMove);
                    document.removeEventListener('mouseup', handleMouseUp);
                }}
            }}

            function setupEnergyPlot() {{
                const plotCtx = document.getElementById('energy-plot').getContext('2d');
                const energyData = trajectoryData.map((frame, i) => ({{x: i, y: frame.energy}}));
                energyPlot = new Chart(plotCtx, {{
                    type: 'scatter',
                    data: {{
                        datasets: [
                            {{ label: 'Energy', data: energyData, backgroundColor: 'rgba(54, 162, 235, 0.6)', borderColor: 'rgba(54, 162, 235, 1)', showLine: true, pointRadius: 3, order: 1 }},
                            {{ label: 'Current Frame', data: [], backgroundColor: 'rgba(255, 99, 132, 1)', pointRadius: 6, pointStyle: 'rectRot', order: 0 }}
                        ]
                    }},
                    options: {{
                        scales: {{
                            x: {{ title: {{ display: true, text: 'Frame Index', color: '#A0AEC0' }}, ticks: {{ color: '#A0AEC0' }}, grid: {{ color: '#4A5568' }} }},
                            y: {{ title: {{ display: true, text: 'Energy [eV]', color: '#A0AEC0' }}, ticks: {{ color: '#A0AEC0' }}, grid: {{ color: '#4A5568' }} }}
                        }},
                        plugins: {{ legend: {{ display: false }} }},
                        onClick: (event, elements) => {{ if (elements.length > 0 && elements[0].datasetIndex === 0) updateScene(elements[0].index); }}
                    }}
                }});
            }}

            function updatePlotHighlight(frameIndex) {{
                if (!energyPlot || !trajectoryData[frameIndex] || trajectoryData[frameIndex].energy === null) return;
                const point = trajectoryData[frameIndex];
                energyPlot.data.datasets[1].data = [{{x: frameIndex, y: point.energy}}];
                energyPlot.update('none');
            }}

            function copyXyzToClipboard() {{
                const frameData = trajectoryData[currentFrameIndex];
                if (!frameData) return;
                const nAtoms = frameData.symbols.length;
                const comment = `Frame ${{currentFrameIndex}}, Energy = ${{frameData.energy ? frameData.energy.toFixed(6) : "N/A"}} eV`;
                let xyzString = `${{nAtoms}}\n${{comment}}\n`;
                for (let i = 0; i < nAtoms; i++) {{
                    const symbol = frameData.symbols[i];
                    const [x, y, z] = frameData.positions[i];
                    xyzString += `${{symbol.padEnd(4)}} ${{x.toFixed(8).padStart(12)}} ${{y.toFixed(8).padStart(12)}} ${{z.toFixed(8).padStart(12)}}\n`;
                }}
                navigator.clipboard.writeText(xyzString).then(() => {{
                    const originalContent = copyXyzBtn.innerHTML;
                    const checkIcon = `<svg width="16" height="16" viewBox="0 0 24 24"><path fill="none" stroke="currentColor" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" d="M20 6L9 17l-5-5"></path></svg>`;
                    copyXyzBtn.innerHTML = checkIcon;
                    setTimeout(() => {{ copyXyzBtn.innerHTML = originalContent; }}, 2000);
                }});
            }}
            
            async function saveAsGif() {{
                const btn = document.getElementById('save-gif-btn');
                btn.classList.add('saving');
                btn.querySelector('i').className = "fas fa-spinner";
                btn.querySelector('span').textContent = 'Recording...';
                btn.disabled = true;

                let workerUrl;
                try {{
                    const workerScriptResponse = await fetch('https://cdnjs.cloudflare.com/ajax/libs/gif.js/0.2.0/gif.worker.js');
                    if (!workerScriptResponse.ok) throw new Error('Failed to fetch worker script');
                    const workerScriptText = await workerScriptResponse.text();
                    const workerBlob = new Blob([workerScriptText], {{ type: 'application/javascript' }});
                    workerUrl = URL.createObjectURL(workerBlob);

                    const wasPlaying = isPlaying;
                    if (wasPlaying) togglePlay();

                    const gif = new GIF({{
                        workers: 2, quality: 10,
                        width: renderer.domElement.width, height: renderer.domElement.height,
                        workerScript: workerUrl, background: scene.background.getStyle(),
                    }});

                    gif.on('finished', function(blob) {{
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = `trajectory.gif`;
                        document.body.appendChild(a);
                        a.click();
                        document.body.removeChild(a);
                        URL.revokeObjectURL(url);
                        URL.revokeObjectURL(workerUrl);

                        btn.classList.remove('saving');
                        btn.querySelector('i').className = "fas fa-file-download";
                        btn.querySelector('span').textContent = 'Save as GIF';
                        btn.disabled = false;
                        if (wasPlaying) togglePlay();
                    }});

                    const cameraState = camera.clone();
                    let frame = 0;
                    const nFrames = (trajectoryData.length > 1) ? trajectoryData.length : 1;
                    const animationSpeed = 1000 / parseInt(document.getElementById('animation-speed-slider').value);

                    function recordFrame() {{
                        if (frame < nFrames) {{
                            updateScene(frame, true);
                            camera.position.copy(cameraState.position);
                            camera.quaternion.copy(cameraState.quaternion);
                            camera.zoom = cameraState.zoom;
                            camera.updateProjectionMatrix();

                            renderer.render(scene, camera);
                            gif.addFrame(renderer.domElement, {{ copy: true, delay: animationSpeed }});
                            frame++;
                            setTimeout(recordFrame, 10);
                        }} else {{
                            gif.render();
                        }}
                    }}
                    recordFrame();

                }} catch (error) {{
                    console.error("Failed to create GIF:", error);
                    btn.classList.remove('saving');
                    btn.querySelector('i').className = "fas fa-exclamation-triangle";
                    btn.querySelector('span').textContent = 'Error!';
                    btn.disabled = false;
                    if(workerUrl) URL.revokeObjectURL(workerUrl);
                }}
            }}
        </script>
    </body>
    </html>
    """
    if write_html:
        with open(write_html, 'w', encoding='utf-8') as f:
            f.write(html_template)

    return HTML(html_template)


def get_adjacency_matrix(atoms:ase.Atoms, covalent_radius_percent:float=108.)->np.ndarray:
    def _covalent_radii(element: str, percent: float):
        radius = covalent_radii.get(element, covalent_radii['default'])
        radius *= (percent / 100)
        return radius
    coordinates = atoms.positions
    symbols = atoms.get_chemical_symbols()
    L2_matrix = squareform(pdist(coordinates, 'euclidean'))
    radii_vector = np.array([_covalent_radii(symbol, covalent_radius_percent) for symbol in symbols])
    radii_sum_matrix = np.add.outer(radii_vector, radii_vector)
    adjacency_matrix = np.array(L2_matrix <= radii_sum_matrix).astype(int)
    np.fill_diagonal(adjacency_matrix, 0)
    return adjacency_matrix

def get_bonds_from_adjacency(adjacency_matrix: np.ndarray) -> List[List[int]]:
    bonds = []
    rows, cols = np.where(adjacency_matrix == 1)
    for i, j in zip(rows, cols):
        if i < j:
            bonds.append([i.item(), j.item()])
    return bonds

def atoms2rdkit_2d(atoms:ase.Atoms, covalent_radius_percent:float=108.)->np.ndarray:
    adjacency_matrix = get_adjacency_matrix(
        atoms=atoms,
        covalent_radius_percent=covalent_radius_percent
    )
    mol = Chem.RWMol()
    for atom in atoms:
        rdkit_atom = Chem.Atom(int(atom.number))
        mol.AddAtom(rdkit_atom)

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            if adjacency_matrix[i, j] == 1:
                mol.AddBond(i, j, Chem.rdchem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(mol)
        AllChem.Compute2DCoords(mol)
        coords = mol.GetConformer().GetPositions()
        return coords[:, :2]
    except Exception as e:
        print(f"RDKit Error: {e}")
        return atoms.get_positions()[:,:2]

def fragment_selector(
    atoms: ase.Atoms,
    covalent_radius_percent: float = 120.,
    atom_scale_3d: float = 0.5,
    atom_radius_2d: float = 8.0,
    bond_thickness_3d: float = 0.1,
    bond_thickness_2d: float = 2.0,
    highlight_color: str = '#FBBF24', # Gold/Amber for better contrast
    bg_color: str = '#2D3748', # Dark navy background
    width: int = 1200,
    height: int = 600,
    write_html: Optional[str] = None) -> HTML:
    """
    Visualizes ASE Atoms objects with linked 2D plot and 3D view in a Jupyter Notebook.
    - Adjust Covalent Radius % and recompute bonds from the sidebar.
    - Resize 2D/3D viewers.
    - Select and highlight atoms in 2D/3D views (rectangle, lasso).
    - Copy selected and residual atom index lists.
    """
    # 1. Prepare data
    positions_3d = atoms.get_positions()
    positions_2d = atoms2rdkit_2d(atoms, covalent_radius_percent=covalent_radius_percent)
    adjacency_matrix = get_adjacency_matrix(atoms, covalent_radius_percent)
    symbols = atoms.get_chemical_symbols()
    bonds = get_bonds_from_adjacency(adjacency_matrix)

    # Create a combined atom info dictionary for JS
    all_symbols = set(atomic_symbols2hex.keys()) | set(covalent_radii.keys())
    atom_info_js = {
        symbol: {
            "color": atomic_symbols2hex.get(symbol, atomic_symbols2hex['default']),
            "r": covalent_radii.get(symbol, covalent_radii['default'])
        } for symbol in all_symbols
    }

    # 2. Serialize data to JSON
    viewer_data = {
        "symbols": symbols,
        "positions_3d": positions_3d.tolist(),
        "positions_2d": positions_2d.tolist(),
        "bonds": bonds,
        "atom_info": atom_info_js, # Pass combined data
        "covalent_radii": covalent_radii, # Keep for JS re-computation
        "initial_covalent_radius_percent": covalent_radius_percent,
        "config": {
            "atom_scale_3d": atom_scale_3d,
            "atom_radius_2d": atom_radius_2d,
            "bond_thickness_3d": bond_thickness_3d,
            "bond_thickness_2d": bond_thickness_2d,
            "highlight_color": highlight_color,
            "bg_color": bg_color
        }
    }
    json_data = json.dumps(viewer_data)

    # 3. Generate HTML template (Modified: Navy Theme)
    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>2D-3D Linked Molecule Viewer</title>
        <style>
            /* --- MODIFIED: Navy Theme --- */
            body {{
                margin: 0;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                overflow: hidden;
                color: #F9FAFB; /* Light gray text */
            }}
            .main-container {{
                display: flex;
                flex-direction: row;
                width: {width}px;
                height: {height}px;
                background-color: #111827; /* Very dark blue */
                border-radius: 12px;
                border: 1px solid #374151; /* Muted gray-blue border */
                box-sizing: border-box;
            }}
            .sidebar {{
                width: 240px;
                background-color: #1F2937; /* Dark blue panel */
                padding: 15px;
                display: flex;
                flex-direction: column;
                gap: 20px;
                border-right: 1px solid #374151;
            }}
            .sidebar h3 {{
                font-size: 16px;
                margin: 0 0 10px 0;
                color: #9CA3AF; /* Medium gray header text */
                border-bottom: 1px solid #374151;
                padding-bottom: 10px;
            }}
            .control-group {{
                display: flex;
                flex-direction: column;
                gap: 8px;
                font-size: 14px;
            }}
            .control-group label {{
                cursor: pointer;
                display: flex;
                align-items: center;
                gap: 5px;
            }}
            .control-group input[type="radio"], .control-group input[type="range"] {{
                cursor: pointer;
            }}

            /* Toggle Button Styles */
            .toggle-button-group {{
                display: flex;
                gap: 5px;
                margin-bottom: 10px;
            }}
            .toggle-btn {{
                flex: 1;
                padding: 8px 12px;
                border: 1px solid #374151;
                background-color: #374151;
                color: #9CA3AF;
                border-radius: 6px;
                cursor: pointer;
                font-size: 13px;
                transition: all 0.2s ease;
                text-align: center;
                display: flex;
                align-items: center;
                justify-content: center;
                gap: 5px;
            }}
            .toggle-btn:hover {{
                background-color: #4B5563;
                color: #F9FAFB;
            }}
            .toggle-btn.active {{
                background-color: #3B82F6;
                color: white;
                border-color: #2563EB;
                box-shadow: 0 0 0 2px rgba(59, 130, 246, 0.3);
            }}
            .toggle-btn.active:hover {{
                background-color: #2563EB;
            }}

            input[type="range"] {{
                accent-color: #3B82F6;
            }}
            .slider-label {{
                display: flex;
                justify-content: space-between;
                font-size: 12px;
            }}
            #recompute-btn {{
                background-color: #3B82F6; /* Bright blue */
                color: white;
                border: none;
                padding: 8px 12px;
                border-radius: 6px;
                cursor: pointer;
                font-size: 14px;
                transition: background-color 0.2s;
            }}
            #recompute-btn:hover {{
                background-color: #2563EB; /* Darker bright blue */
            }}
            .content-wrapper {{ flex: 1; display: flex; flex-direction: column; min-width: 0; }}
            .viewer-area {{ display: flex; flex: 1; gap: 10px; min-height: 0; padding: 10px; }}
            .panel {{
                flex: 1;
                border: 1px solid #374151;
                border-radius: 8px;
                overflow: hidden;
                position: relative;
            }}
            #plot-2d-container, #canvas-3d-container {{ min-width: 200px; background-color: {bg_color}; }}
            #canvas-2d {{ display: block; }}
            #selection-canvas-3d {{ position: absolute; top: 0; left: 0; pointer-events: none; }}
            .console-area-wrapper {{ height: 80px; padding: 0 10px 10px 10px; }}
            .console-area {{
                height: 100%; box-sizing: border-box; background-color: #1F2937;
                border-radius: 8px; border: 1px solid #374151; font-family: 'Menlo', 'Monaco', monospace;
                font-size: 14px; overflow-y: auto; flex-direction: column; align-items: stretch;
                justify-content: space-around; padding: 5px 10px; display: flex;
            }}
            .console-line {{ display: flex; justify-content: space-between; align-items: center; gap: 10px; }}
            #console-output, #console-output-residual {{
                flex-grow: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
            }}
            #copy-btn, #copy-btn-residual {{
                background-color: #374151; color: #F9FAFB; border: none; padding: 6px 10px;
                border-radius: 5px; cursor: pointer; transition: background-color 0.2s; flex-shrink: 0;
            }}
            #copy-btn:hover, #copy-btn-residual:hover {{ background-color: #3B82F6; }}
        </style>
    </head>
    <body>
        <div class="main-container">
            <div class="sidebar">
                 <div class="control-group">
                    <h3>3D Mode</h3>
                    <div class="toggle-button-group">
                        <button class="toggle-btn active" data-mode="rotate">🔄 Rotate</button>
                        <button class="toggle-btn" data-mode="select">👆 Select</button>
                    </div>
                </div>
                <div class="control-group">
                    <h3>Select Mode</h3>
                    <div class="toggle-button-group">
                        <button class="toggle-btn" data-tool="lasso">➰ Lasso</button>
                        <button class="toggle-btn active" data-tool="rectangle">🔲 Rectangle</button>
                    </div>
                </div>
                <div class="control-group">
                    <h3>Bond Calculation</h3>
                    <label for="covalent-radius-slider">Covalent Radius (%): <span id="slider-value">{covalent_radius_percent}</span>%</label>
                    <input type="range" id="covalent-radius-slider" min="50" max="150" value="{covalent_radius_percent}">
                    <button id="recompute-btn">Recompute Bonds</button>
                </div>
            </div>
            <div class="content-wrapper">
                <div class="viewer-area">
                    <div id="plot-2d-container" class="panel"></div>
                    <div id="canvas-3d-container" class="panel"></div>
                </div>
                <div class="console-area-wrapper">
                    <div class="console-area">
                        <div class="console-line">
                            <span id="console-output" title="Selected atoms: []">Selected atoms: []</span>
                            <button id="copy-btn" title="Copy selected index list">Copy</button>
                        </div>
                        <div class="console-line">
                            <span id="console-output-residual" title="Residual atoms: []">Residual atoms: []</span>
                            <button id="copy-btn-residual" title="Copy residual index list">Copy</button>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>

        <script>
            // --- Data and Config ---
            const viewerData = {json_data};
            const config = viewerData.config;
            // Atom info is now passed from Python
            const atomInfo = viewerData.atom_info;

            // --- Global State ---
            let selectedIndices = new Set();
            let primaryMode = 'rotate'; // 'rotate' or 'select'
            let selectTool = 'rectangle'; // 'rectangle' or 'lasso'

            // --- 2D Plot Globals ---
            let canvas2D, ctx2D, positions2D_canvas = [];
            let isDrawing2D = false, selectionRect = {{}}, lassoPath = [];

            // --- 3D View Globals ---
            let scene3D, camera3D, renderer3D, controls3D, atomGroup3D, raycaster3D;
            let canvas3DOverlay, ctx3DOverlay;
            let isDrawing3D = false, selectionRect3D = {{}}, lassoPath3D = [];

            // --- UI Elements ---
            const consoleOutput = document.getElementById('console-output');
            const consoleOutputResidual = document.getElementById('console-output-residual');
            const container2D = document.getElementById('plot-2d-container');
            const container3D = document.getElementById('canvas-3d-container');
            const copyBtn = document.getElementById('copy-btn');
            const copyBtnResidual = document.getElementById('copy-btn-residual');
            const covalentRadiusSlider = document.getElementById('covalent-radius-slider');
            const sliderValueSpan = document.getElementById('slider-value');
            const recomputeBtn = document.getElementById('recompute-btn');

            // --- Initialization ---
            function init() {{
                init2D();
                init3D();
                setupEventListeners();
                controls3D.enabled = true;
                updateSelection();
            }}

            // --- Bond Recomputation ---
            function recomputeBondsJS(percent) {{
                const positions = viewerData.positions_3d;
                const symbols = viewerData.symbols;
                const covalentRadii = viewerData.covalent_radii; // Uses the full radii dict
                const defaultRadius = covalentRadii['default'] || 0.7;
                const getRadius = (symbol) => (covalentRadii[symbol] || defaultRadius) * (percent / 100);
                const radiiVector = symbols.map(s => getRadius(s));
                const newBonds = [];
                for (let i = 0; i < positions.length; i++) {{
                    for (let j = i + 1; j < positions.length; j++) {{
                        const p1 = positions[i], p2 = positions[j];
                        const dist = Math.sqrt(Math.pow(p1[0]-p2[0],2) + Math.pow(p1[1]-p2[1],2) + Math.pow(p1[2]-p2[2],2));
                        const radiiSum = radiiVector[i] + radiiVector[j];
                        if (dist <= radiiSum) newBonds.push([i, j]);
                    }}
                }}
                viewerData.bonds = newBonds;
                draw2D();
                draw3D();
            }}

            // --- 2D Plot ---
            function init2D() {{
                canvas2D = document.createElement('canvas');
                canvas2D.id = 'canvas-2d';
                ctx2D = canvas2D.getContext('2d');
                container2D.appendChild(canvas2D);
                const resizeObserver = new ResizeObserver(() => setTimeout(resizeCanvas2D, 16));
                resizeObserver.observe(container2D);
                resizeCanvas2D();
            }}

            function resizeCanvas2D() {{
                const rect = container2D.getBoundingClientRect();
                if (rect.width === 0 || rect.height === 0) return;
                const dpr = window.devicePixelRatio || 1;
                canvas2D.width = rect.width * dpr;
                canvas2D.height = rect.height * dpr;
                canvas2D.style.width = rect.width + 'px';
                canvas2D.style.height = rect.height + 'px';
                ctx2D.scale(dpr, dpr);
                map2DCoordinates(rect.width, rect.height);
                draw2D();
            }}

            function map2DCoordinates(canvasWidth, canvasHeight) {{
                const {{positions_2d}} = viewerData;
                if (positions_2d.length === 0) return;
                const min = [Infinity, Infinity], max = [-Infinity, -Infinity];
                positions_2d.forEach(p => {{
                    min[0] = Math.min(min[0], p[0]); min[1] = Math.min(min[1], p[1]);
                    max[0] = Math.max(max[0], p[0]); max[1] = Math.max(max[1], p[1]);
                }});
                const dataWidth = max[0] - min[0] || 1, dataHeight = max[1] - min[1] || 1;
                const padding = 20;
                const scaleX = (canvasWidth - padding * 2) / dataWidth;
                const scaleY = (canvasHeight - padding * 2) / dataHeight;
                const scale2D = Math.min(scaleX, scaleY);
                const offsetX = (canvasWidth - dataWidth * scale2D) / 2 - min[0] * scale2D;
                const offsetY = (canvasHeight - dataHeight * scale2D) / 2 - min[1] * scale2D;
                positions2D_canvas = positions_2d.map(p => ({{ x: p[0] * scale2D + offsetX, y: p[1] * scale2D + offsetY }}));
            }}

            function draw2D() {{
                if (!ctx2D) return;
                const rect = container2D.getBoundingClientRect();
                ctx2D.clearRect(0, 0, rect.width, rect.height);
                ctx2D.strokeStyle = '#6B7280'; ctx2D.lineWidth = config.bond_thickness_2d;
                viewerData.bonds.forEach(bond => {{
                    if(!positions2D_canvas[bond[0]] || !positions2D_canvas[bond[1]]) return;
                    const p1 = positions2D_canvas[bond[0]], p2 = positions2D_canvas[bond[1]];
                    ctx2D.beginPath(); ctx2D.moveTo(p1.x, p1.y); ctx2D.lineTo(p2.x, p2.y); ctx2D.stroke();
                }});
                positions2D_canvas.forEach((pos, i) => {{
                    const symbol = viewerData.symbols[i];
                    const info = atomInfo[symbol] || atomInfo['default'];
                    ctx2D.beginPath(); ctx2D.arc(pos.x, pos.y, config.atom_radius_2d, 0, 2 * Math.PI);
                    ctx2D.fillStyle = info.color; ctx2D.fill();
                    if (selectedIndices.has(i)) {{
                        ctx2D.strokeStyle = config.highlight_color; ctx2D.lineWidth = 3; ctx2D.stroke();
                    }}
                }});
                if (isDrawing2D) {{
                    ctx2D.fillStyle = 'rgba(251, 191, 36, 0.2)';
                    ctx2D.strokeStyle = 'rgba(251, 191, 36, 0.8)';
                    ctx2D.lineWidth = 1;
                    if (selectTool === 'rectangle') {{
                        const w = selectionRect.x2 - selectionRect.x1, h = selectionRect.y2 - selectionRect.y1;
                        ctx2D.fillRect(selectionRect.x1, selectionRect.y1, w, h); ctx2D.strokeRect(selectionRect.x1, selectionRect.y1, w, h);
                    }} else if (selectTool === 'lasso' && lassoPath.length > 1) {{
                        ctx2D.beginPath(); ctx2D.moveTo(lassoPath[0].x, lassoPath[0].y);
                        for (let i = 1; i < lassoPath.length; i++) ctx2D.lineTo(lassoPath[i].x, lassoPath[i].y);
                        ctx2D.stroke();
                    }}
                }}
            }}

            // --- 3D View ---
            function init3D() {{
                scene3D = new THREE.Scene();
                scene3D.background = new THREE.Color(config.bg_color);
                camera3D = new THREE.PerspectiveCamera(75, container3D.clientWidth / container3D.clientHeight, 0.1, 1000);
                const com = new THREE.Vector3();
                if (viewerData.positions_3d.length > 0) {{
                    viewerData.positions_3d.forEach(p => com.add(new THREE.Vector3(...p)));
                    com.divideScalar(viewerData.positions_3d.length);
                }}
                camera3D.position.z = com.z + 20;
                renderer3D = new THREE.WebGLRenderer({{ antialias: true }});
                renderer3D.setPixelRatio(window.devicePixelRatio);
                container3D.appendChild(renderer3D.domElement);
                canvas3DOverlay = document.createElement('canvas');
                canvas3DOverlay.id = 'selection-canvas-3d';
                ctx3DOverlay = canvas3DOverlay.getContext('2d');
                container3D.appendChild(canvas3DOverlay);
                scene3D.add(new THREE.AmbientLight(0xcccccc, 0.8));
                const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
                dirLight.position.set(5, 10, 7.5); scene3D.add(dirLight);
                controls3D = new THREE.OrbitControls(camera3D, renderer3D.domElement);
                controls3D.target.copy(com);
                raycaster3D = new THREE.Raycaster();
                const resizeObserver = new ResizeObserver(() => {{
                    const rect = container3D.getBoundingClientRect();
                    if (rect.height === 0 || rect.width === 0) return;
                    camera3D.aspect = rect.width / rect.height;
                    camera3D.updateProjectionMatrix();
                    renderer3D.setSize(rect.width, rect.height);
                    const dpr = window.devicePixelRatio || 1;
                    canvas3DOverlay.width = rect.width * dpr;
                    canvas3DOverlay.height = rect.height * dpr;
                    canvas3DOverlay.style.width = rect.width + 'px';
                    canvas3DOverlay.style.height = rect.height + 'px';
                    ctx3DOverlay.scale(dpr, dpr);
                }});
                resizeObserver.observe(container3D);
                const initialRect = container3D.getBoundingClientRect();
                renderer3D.setSize(initialRect.width, initialRect.height);
                atomGroup3D = new THREE.Group();
                scene3D.add(atomGroup3D);
                draw3D();
                function animate() {{ requestAnimationFrame(animate); controls3D.update(); renderer3D.render(scene3D, camera3D); }}
                animate();
            }}

            function draw3D() {{
                while (atomGroup3D.children.length) {{
                    const child = atomGroup3D.children[0];
                    atomGroup3D.remove(child);
                    if (child.geometry) child.geometry.dispose();
                    if (child.material && typeof child.material.dispose === 'function') child.material.dispose();
                }}
                const {{ positions_3d, symbols, bonds }} = viewerData;
                const materials = {{}};
                const highlightMaterial = new THREE.MeshStandardMaterial({{ color: config.highlight_color, emissive: config.highlight_color, emissiveIntensity: 0.5 }});
                symbols.forEach(s => {{ if (!materials[s]) materials[s] = new THREE.MeshLambertMaterial({{ color: (atomInfo[s] || atomInfo['default']).color }}); }});
                positions_3d.forEach((pos, i) => {{
                    const info = atomInfo[symbols[i]] || atomInfo['default'];
                    const geometry = new THREE.SphereGeometry(info.r * config.atom_scale_3d, 32, 32);
                    const material = selectedIndices.has(i) ? highlightMaterial : materials[symbols[i]];
                    const sphere = new THREE.Mesh(geometry, material);
                    sphere.position.set(pos[0], pos[1], pos[2]);
                    sphere.userData.index = i;
                    atomGroup3D.add(sphere);
                }});
                bonds.forEach(bond => {{
                    const p1 = new THREE.Vector3(...positions_3d[bond[0]]), p2 = new THREE.Vector3(...positions_3d[bond[1]]);
                    const mid = p1.clone().add(p2).multiplyScalar(0.5);
                    const mat1 = selectedIndices.has(bond[0]) ? highlightMaterial : materials[symbols[bond[0]]];
                    const mat2 = selectedIndices.has(bond[1]) ? highlightMaterial : materials[symbols[bond[1]]];
                    atomGroup3D.add(new THREE.Mesh(new THREE.TubeGeometry(new THREE.LineCurve3(p1, mid), 1, config.bond_thickness_3d, 8, false), mat1));
                    atomGroup3D.add(new THREE.Mesh(new THREE.TubeGeometry(new THREE.LineCurve3(mid, p2), 1, config.bond_thickness_3d, 8, false), mat2));
                }});
            }}

            function draw3DSelectionTool() {{
                if (!ctx3DOverlay) return;
                const rect = canvas3DOverlay.getBoundingClientRect();
                ctx3DOverlay.clearRect(0, 0, rect.width, rect.height);
                if (!isDrawing3D) return;
                ctx3DOverlay.fillStyle = 'rgba(251, 191, 36, 0.2)';
                ctx3DOverlay.strokeStyle = 'rgba(251, 191, 36, 0.8)';
                ctx3DOverlay.lineWidth = 1;
                if (selectTool === 'rectangle') {{
                    const w = selectionRect3D.x2 - selectionRect3D.x1, h = selectionRect3D.y2 - selectionRect3D.y1;
                    ctx3DOverlay.fillRect(selectionRect3D.x1, selectionRect3D.y1, w, h);
                    ctx3DOverlay.strokeRect(selectionRect3D.x1, selectionRect3D.y1, w, h);
                }} else if (selectTool === 'lasso' && lassoPath3D.length > 1) {{
                    ctx3DOverlay.beginPath();
                    ctx3DOverlay.moveTo(lassoPath3D[0].x, lassoPath3D[0].y);
                    for (let i = 1; i < lassoPath3D.length; i++) ctx3DOverlay.lineTo(lassoPath3D[i].x, lassoPath3D[i].y);
                    ctx3DOverlay.stroke();
                }}
            }}

            // --- Event Handling & Selection ---
            function setupEventListeners() {{
                // Toggle button handlers
                document.querySelectorAll('[data-mode]').forEach(btn => {{
                    btn.addEventListener('click', (e) => {{
                        primaryMode = e.target.dataset.mode;
                        controls3D.enabled = primaryMode === 'rotate';
                        canvas3DOverlay.style.pointerEvents = (primaryMode === 'select') ? 'auto' : 'none';

                        // Update toggle button states
                        document.querySelectorAll('[data-mode]').forEach(b => b.classList.remove('active'));
                        e.target.classList.add('active');
                    }});
                }});

                document.querySelectorAll('[data-tool]').forEach(btn => {{
                    btn.addEventListener('click', (e) => {{
                        selectTool = e.target.dataset.tool;

                        // Update toggle button states
                        document.querySelectorAll('[data-tool]').forEach(b => b.classList.remove('active'));
                        e.target.classList.add('active');
                    }});
                }});

                covalentRadiusSlider.addEventListener('input', (e) => sliderValueSpan.textContent = e.target.value);
                recomputeBtn.addEventListener('click', () => recomputeBondsJS(parseFloat(covalentRadiusSlider.value)));
                copyBtn.addEventListener('click', () => {{
                    const textToCopy = consoleOutput.title.replace('Selected atoms: ', '');
                    navigator.clipboard.writeText(textToCopy).then(() => {{ copyBtn.textContent = 'Copied!'; setTimeout(() => copyBtn.textContent = 'Copy', 1500); }});
                }});
                copyBtnResidual.addEventListener('click', () => {{
                    const textToCopy = consoleOutputResidual.title.replace('Residual atoms: ', '');
                    navigator.clipboard.writeText(textToCopy).then(() => {{ copyBtnResidual.textContent = 'Copied!'; setTimeout(() => copyBtnResidual.textContent = 'Copy', 1500); }});
                }});
                canvas2D.addEventListener('mousedown', handle2DMouseDown);
                canvas2D.addEventListener('mousemove', handle2DMouseMove);
                canvas2D.addEventListener('mouseup', handle2DMouseUp);
                canvas2D.addEventListener('mouseleave', () => {{ if (isDrawing2D) {{ isDrawing2D = false; lassoPath = []; draw2D(); }} }});
                canvas3DOverlay.addEventListener('mousedown', handle3DMouseDown);
                canvas3DOverlay.addEventListener('mousemove', handle3DMouseMove);
                canvas3DOverlay.addEventListener('mouseup', handle3DMouseUp);
                canvas3DOverlay.addEventListener('mouseleave', () => {{ if (isDrawing3D) {{ isDrawing3D = false; lassoPath3D = []; draw3DSelectionTool(); }} }});
            }}
            function handle2DMouseDown(e) {{ const rect = canvas2D.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; isDrawing2D = true; if (selectTool === 'rectangle') selectionRect = {{ x1: x, y1: y, x2: x, y2: y }}; else lassoPath = [{{x, y}}]; }}
            function handle2DMouseMove(e) {{ if (!isDrawing2D) return; const rect = canvas2D.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; if (selectTool === 'rectangle') {{ selectionRect.x2 = x; selectionRect.y2 = y; }} else lassoPath.push({{x, y}}); draw2D(); }}
            function handle2DMouseUp(e) {{ if (!isDrawing2D) return; isDrawing2D = false; selectAtomsIn2D(e); lassoPath = []; updateSelection(); }}
            function handle3DMouseDown(e) {{ if (primaryMode === 'rotate') return; const rect = canvas3DOverlay.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; isDrawing3D = true; if (selectTool === 'rectangle') selectionRect3D = {{ x1: x, y1: y, x2: x, y2: y }}; else lassoPath3D = [{{x, y}}]; }}
            function handle3DMouseMove(e) {{ if (!isDrawing3D) return; const rect = canvas3DOverlay.getBoundingClientRect(); const x = e.clientX - rect.left, y = e.clientY - rect.top; if (selectTool === 'rectangle') {{ selectionRect3D.x2 = x; selectionRect3D.y2 = y; }} else lassoPath3D.push({{x, y}}); draw3DSelectionTool(); }}
            function handle3DMouseUp(e) {{ if (!isDrawing3D) return; isDrawing3D = false; selectAtomsIn3D(e); lassoPath3D = []; draw3DSelectionTool(); }}
            function selectAtomsIn2D(event) {{ const currentSelection = new Set(selectedIndices); if (!event.shiftKey) selectedIndices.clear(); const x_min = Math.min(selectionRect.x1, selectionRect.x2), x_max = Math.max(selectionRect.x1, selectionRect.x2); const y_min = Math.min(selectionRect.y1, selectionRect.y2), y_max = Math.max(selectionRect.y1, selectionRect.y2); const isClick = Math.hypot(x_max - x_min, y_max - y_min) < 5; if (isClick) {{ let closestIndex = -1, minDistance = Infinity; positions2D_canvas.forEach((pos, i) => {{ const dist = Math.hypot(pos.x - x_min, pos.y - y_min); if (dist < config.atom_radius_2d * 1.5 && dist < minDistance) {{ minDistance = dist; closestIndex = i; }} }}); if (closestIndex !== -1) {{ const index = closestIndex; if (event.shiftKey) {{ if (currentSelection.has(index)) selectedIndices.delete(index); else selectedIndices.add(index); }} else {{ if (currentSelection.has(index) && currentSelection.size === 1) {{ /* deselect */ }} else {{ selectedIndices.clear(); selectedIndices.add(index); }} }} }} }} else if (selectTool === 'rectangle') {{ positions2D_canvas.forEach((pos, i) => {{ if (pos.x >= x_min && pos.x <= x_max && pos.y >= y_min && pos.y <= y_max) {{ if (event.shiftKey && currentSelection.has(i)) selectedIndices.delete(i); else selectedIndices.add(i); }} }}); }} else if (selectTool === 'lasso') {{ if (lassoPath.length < 3) return; positions2D_canvas.forEach((pos, i) => {{ if (isPointInPolygon(pos, lassoPath)) {{ if (event.shiftKey && currentSelection.has(i)) selectedIndices.delete(i); else selectedIndices.add(i); }} }}); }} }}
            function selectAtomsIn3D(event) {{ const currentSelection = new Set(selectedIndices); if (!event.shiftKey) selectedIndices.clear(); const rect = canvas3DOverlay.getBoundingClientRect(); const x_min = Math.min(selectionRect3D.x1, selectionRect3D.x2), x_max = Math.max(selectionRect3D.x1, selectionRect3D.x2); const y_min = Math.min(selectionRect3D.y1, selectionRect3D.y2), y_max = Math.max(selectionRect3D.y1, selectionRect3D.y2); const isClick = Math.hypot(x_max - x_min, y_max - y_min) < 5; if (isClick) {{ const mouse = new THREE.Vector2(); mouse.x = (selectionRect3D.x1 / rect.width) * 2 - 1; mouse.y = -(selectionRect3D.y1 / rect.height) * 2 + 1; raycaster3D.setFromCamera(mouse, camera3D); const intersects = raycaster3D.intersectObjects(atomGroup3D.children.filter(c => c.isMesh && c.userData.index !== undefined)); if (intersects.length > 0) {{ const index = intersects[0].object.userData.index; if (event.shiftKey) {{ if (currentSelection.has(index)) selectedIndices.delete(index); else selectedIndices.add(index); }} else {{ if (currentSelection.has(index) && currentSelection.size === 1) {{ /* deselect */ }} else {{ selectedIndices.clear(); selectedIndices.add(index); }} }} }} }} else {{ const frustum = new THREE.Frustum(); frustum.setFromProjectionMatrix(new THREE.Matrix4().multiplyMatrices(camera3D.projectionMatrix, camera3D.matrixWorldInverse)); viewerData.positions_3d.forEach((pos, i) => {{ const worldPoint = new THREE.Vector3(...pos); if (!frustum.containsPoint(worldPoint)) return; const screenPoint = worldPoint.clone().project(camera3D); const x = (screenPoint.x * 0.5 + 0.5) * rect.width; const y = (-screenPoint.y * 0.5 + 0.5) * rect.height; let inShape = false; if (selectTool === 'rectangle') {{ inShape = (x >= x_min && x <= x_max && y >= y_min && y <= y_max); }} else if (selectTool === 'lasso' && lassoPath3D.length > 2) {{ inShape = isPointInPolygon({{x, y}}, lassoPath3D); }} if (inShape) {{ if (event.shiftKey && currentSelection.has(i)) selectedIndices.delete(i); else selectedIndices.add(i); }} }}); }} updateSelection(); }}
            function isPointInPolygon(point, polygon) {{ let x = point.x, y = point.y, inside = false; for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {{ let xi = polygon[i].x, yi = polygon[i].y, xj = polygon[j].x, yj = polygon[j].y; let intersect = ((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi); if (intersect) inside = !inside; }} return inside; }}
            function updateSelection() {{ const selectedArray = Array.from(selectedIndices).sort((a,b)=>a-b); const totalAtoms = viewerData.symbols.length; const allIndices = new Set(Array.from(Array(totalAtoms).keys())); const residualArray = [...allIndices].filter(i => !selectedIndices.has(i)); const selectedText = `[${{selectedArray.join(', ')}}]`; consoleOutput.textContent = `Selected atoms: ${{selectedText}}`; consoleOutput.title = `Selected atoms: ${{selectedText}}`; const residualText = `[${{residualArray.join(', ')}}]`; consoleOutputResidual.textContent = `Residual atoms: ${{residualText}}`; consoleOutputResidual.title = `Residual atoms: ${{residualText}}`; draw2D(); draw3D(); }}

            // --- Start the application ---
            init();
        </script>
    </body>
    </html>
    """
    if write_html:
        with open(write_html, 'w', encoding='utf-8') as f:
            f.write(html_template)

    return HTML(html_template)


def normal_mode_viewer(
    atoms: Atoms,
    vib,  # ASE Vibrations object
    amplitude: float = 0.75,
    n_frames: int = 40,
    initial_mode_index: int = 1,
    rotation_mode: str = "orbit",
    bg_color: str = "#2D3748",
    show_bonds: bool = True,
    show_cell: bool = True,
    show_vectors: bool = False,
    atom_scale: float = 1.0,
    bond_cutoff_scale: float = 1.2,
    bond_thickness: float = 0.08,
    animation_fps: int = 30,
    width: int = 800,
    height: int = 600,
    write_html: Optional[str] = None,
) -> HTML:
    """
    Visualizes Normal Modes in 3D within a Jupyter Notebook.

    Args:
        atoms: ASE Atoms object of the equilibrium structure.
        vib: ASE Vibrations object.
        amplitude: Initial amplitude for the vibrational animation.
        n_frames: Number of frames per mode in the animation.
        rotation_mode: Camera rotation mode. 'orbit' or 'trackball'.
        initial_mode_index: Index of the vibrational mode to display initially (default: 1).
        bg_color: Background color (HEX code).
        show_bonds: Initial state of the 'Bonds' checkbox.
        show_cell: Initial state of the 'Unit Cell' checkbox.
        show_vectors: Initial state of the 'Mode Vectors' checkbox.
        atom_scale: Initial scale value for atom radii.
        bond_cutoff_scale: Initial scale factor for bond distance calculation.
        bond_thickness: Initial thickness value for bonds.
        animation_fps: Frames per second for the animation.
        width: Width of the viewer (px).
        height: Height of the viewer (px).
        write_html: HTML file to save the output to.

    Returns:
        IPython.display.HTML object.
    """

    assert initial_mode_index > 0, "initial_mode_index starts from 1"
    initial_mode_index -= 1 # internally starts from 0

    # Get vibrational data
    mode_vectors = vib.get_vibrations().get_modes(all_atoms=True)
    freqs = vib.get_frequencies()

    # Initial atom data
    initial_atoms_data = {
        "symbols": atoms.get_chemical_symbols(),
        "positions": atoms.get_positions().tolist(),
        "cell": atoms.get_cell().tolist(),
    }
    has_initial_cell = atoms.get_cell().any()

    # Serialize data for JavaScript
    json_initial_atoms = json.dumps(initial_atoms_data)
    json_mode_vectors = json.dumps(mode_vectors.tolist())

    # Process frequency data for display
    formatted_freqs = []
    is_imaginary_list = [bool(np.iscomplexobj(f) and f.imag != 0) for f in freqs]
    for i, f in enumerate(freqs):
        if is_imaginary_list[i]:
            formatted_freqs.append(f"{f.imag:.2f}i")
        else:
            formatted_freqs.append(f"{f.real:.2f}")

    json_freqs = json.dumps(formatted_freqs)
    json_is_imaginary = json.dumps(is_imaginary_list)

    # External atom data
    json_atomic_colors = json.dumps(atomic_symbols2hex)
    json_covalent_radii = json.dumps(covalent_radii)

    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Normal Mode Viewer</title>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
        <style>
            body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; overflow: hidden;}}
            .viewer-container {{ position: relative; width: {width}px; height: {height}px; display: flex; background-color: #2D3748; border-radius: 8px; overflow: hidden; border: 1px solid #4A5568; }}

            .sidebar {{
                width: 260px;
                padding: 20px;
                background-color: #1A202C;
                color: #E2E8F0;
                display: flex;
                flex-direction: column;
                gap: 16px;
                overflow-y: auto;
                transition: width 0.3s ease, padding 0.3s ease;
                flex-shrink: 0;
                -ms-overflow-style: none;
                scrollbar-width: none;
            }}
            .sidebar::-webkit-scrollbar {{ display: none; }}

            .sidebar.collapsed {{
                width: 0;
                padding: 0;
                overflow: hidden;
            }}

            .sidebar h2 {{ font-size: 18px; margin: 0; padding-bottom: 10px; border-bottom: 1px solid #4A5568; color: #A0AEC0; }}
            .sidebar h3 {{ font-size: 1rem; margin: 0; }}
            .sidebar .control-group, .sidebar details {{ display: flex; flex-direction: column; gap: 12px; }}
            select {{ background-color: #2D3748; color: #E2E8F0; border: 1px solid #4A5568; border-radius: 4px; padding: 4px; width: 100%; }}
            input[type="color"] {{ width: 100%; height: 25px; border: 1px solid #4A5568; border-radius: 4px; padding: 2px; background-color: #2D3748;}}
            .animation-controls {{ display: flex; align-items: center; gap: 10px; }}
            .copy-button {{ background-color: #4A5568; border: none; border-radius: 5px; color: white; padding: 6px 10px; cursor: pointer; font-size: 12px; text-align: center; transition: background-color 0.2s; }}
            .copy-button:hover {{ background-color: #636e82; }}
            #play-pause-btn {{ background-color: #4A5568; border: none; border-radius: 50%; width: 32px; height: 32px; cursor: pointer; display: flex; justify-content: center; align-items: center; flex-shrink: 0; transition: background-color 0.2s; }}
            #play-pause-btn:hover {{ background-color: #636e82; }}
            #play-pause-btn svg {{ fill: white; }}
            input[type="range"] {{ width: 100%; }}
            .slider-group {{ display: flex; flex-direction: column; gap: 4px;}}
            .slider-label, .slider-value {{ font-size: 12px; color: #A0AEC0; }}

            #toggle-sidebar-btn {{ position: absolute; top: 15px; left: 15px; z-index: 100; background: rgba(45, 55, 72, 0.9); color: #E2E8F0; border: 1px solid #4A5568; border-radius: 6px; width: 36px; height: 36px; cursor: pointer; display: flex; flex-direction: column; justify-content: center; align-items: center; gap: 3px; transition: all 0.3s ease; backdrop-filter: blur(5px); user-select: none; }}
            #toggle-sidebar-btn:hover {{ background: rgba(45, 55, 72, 1); border-color: #718096; }}
            .hamburger-line {{ width: 18px; height: 2px; background-color: #E2E8F0; border-radius: 1px; transition: all 0.3s ease; }}
            #toggle-sidebar-btn.open .hamburger-line:nth-child(1) {{ transform: rotate(45deg) translate(4px, 4px); }}
            #toggle-sidebar-btn.open .hamburger-line:nth-child(2) {{ opacity: 0; }}
            #toggle-sidebar-btn.open .hamburger-line:nth-child(3) {{ transform: rotate(-45deg) translate(4px, -4px); }}

            #canvas-container {{ flex: 1; min-width: 0; height: 100%; cursor: grab; background-color: {bg_color}; }}
            #top-right-controls {{ position: absolute; top: 15px; right: 15px; z-index: 10; }}
            details {{ border-top: 1px solid #4A5568; padding-top: 12px; }}
            summary {{ cursor: pointer; font-size: 1rem; font-weight: bold; padding: 5px 0; list-style-position: inside; }}
            summary::marker {{ color: #A0AEC0; }}
            .details-content {{ padding-top: 10px; display: flex; flex-direction: column; gap: 12px; }}
            .toggle-switch {{ position: relative; display: flex; align-items: center; justify-content: space-between; width: 100%;}}
            .toggle-switch-label {{ font-size: 14px; cursor: pointer; user-select: none; }}
            .toggle-switch input {{ opacity: 0; width: 0; height: 0; }}
            .switch-slider {{ position: relative; cursor: pointer; width: 34px; height: 20px; background-color: #4A5568; border-radius: 34px; transition: .4s; }}
            .switch-slider:before {{ position: absolute; content: ""; height: 14px; width: 14px; left: 3px; bottom: 3px; background-color: white; border-radius: 50%; transition: .4s; }}
            input:checked + .switch-slider {{ background-color: #4299E1; }}
            input:checked + .switch-slider:before {{ transform: translateX(14px); }}

            #save-gif-btn {{ background-color: rgba(74, 85, 104, 0.7); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 8px; color: #E2E8F0; padding: 10px 16px; cursor: pointer; font-size: 14px; font-weight: 600; display: flex; align-items: center; gap: 8px; transition: all 0.3s ease; backdrop-filter: blur(5px); }}
            #save-gif-btn:hover {{ background-color: rgba(99, 110, 130, 0.85); transform: translateY(-2px); }}
            #save-gif-btn:active {{ transform: translateY(0); }}
            #save-gif-btn.saving {{ background: #6B7280; cursor: not-allowed; transform: none; }}
        </style>
    </head>
    <body>
        <div class="viewer-container">
            <div id="toggle-sidebar-btn" class="open">
                <div class="hamburger-line"></div><div class="hamburger-line"></div><div class="hamburger-line"></div>
            </div>

            <div class="sidebar" id="sidebar">
                <h2>Normal Mode Viewer</h2>
                <div class="control-group">
                    <h3>Vibrational Mode</h3>
                    <select id="mode-select"></select>
                </div>
                <div id="animation-section" class="control-group">
                    <h3>Animation</h3>
                    <div class="animation-controls">
                        <button id="play-pause-btn">
                            <svg id="play-icon" width="14" height="14" viewBox="0 0 24 24"><path d="M8 5v14l11-7z"/></svg>
                            <svg id="pause-icon" width="14" height="14" viewBox="0 0 24 24" style="display:none;"><path d="M6 19h4V5H6v14zm8-14v14h4V5h-4z"/></svg>
                        </button>
                        <input type="range" id="frame-slider" min="0" max="{n_frames - 1}" value="0">
                    </div>
                    <span id="frame-label" class="slider-value">Frame: 0 / {n_frames - 1}</span>
                    <button id="copy-xyz-btn" class="copy-button">Copy Trajectory (XYZ)</button>
                </div>
                <details open>
                    <summary>Display Options</summary>
                    <div class="details-content">
                        <label class="toggle-switch">
                            <span class="toggle-switch-label">Bonds</span>
                            <input type="checkbox" id="bonds-toggle" {'checked' if show_bonds else ''}>
                            <span class="switch-slider"></span>
                        </label>
                        <label class="toggle-switch" style="display: {'flex' if has_initial_cell else 'none'}">
                            <span class="toggle-switch-label">Unit Cell</span>
                            <input type="checkbox" id="cell-toggle" {'checked' if show_cell else ''}>
                            <span class="switch-slider"></span>
                        </label>
                        <label class="toggle-switch">
                            <span class="toggle-switch-label">Mode Vectors</span>
                            <input type="checkbox" id="vectors-toggle" {'checked' if show_vectors else ''}>
                            <span class="switch-slider"></span>
                        </label>
                    </div>
                </details>
                <details>
                    <summary>Advanced Settings</summary>
                    <div class="details-content">
                        <div class="slider-group">
                            <label for="style-select" class="slider-label">Style</label>
                            <select id="style-select">
                                <option value="lambert">Standard</option><option value="glossy">Glossy</option>
                                <option value="metallic">Metallic</option><option value="toon">Toon</option>
                                <option value="2d">2D</option><option value="neon">Neon</option>
                            </select>
                        </div>
                        <div class="slider-group">
                            <label for="rotation-mode-select" class="slider-label">Rotation Mode</label>
                            <select id="rotation-mode-select">
                                <option value="orbit" {'selected' if rotation_mode == 'orbit' else ''}>Orbit</option>
                                <option value="trackball" {'selected' if rotation_mode == 'trackball' else ''}>Trackball</option>
                            </select>
                        </div>
                        <div class="slider-group">
                            <label for="amplitude-slider" class="slider-label">Displacement Amplitude</label>
                            <input type="range" id="amplitude-slider" min="0.1" max="2.5" step="0.05" value="{amplitude}">
                            <span id="amplitude-label" class="slider-value">Scale: {amplitude:.2f}</span>
                        </div>
                        <div class="slider-group"><label for="atom-scale-slider" class="slider-label">Atom Size</label><input type="range" id="atom-scale-slider" min="0.1" max="2.0" step="0.05" value="{atom_scale}"><span id="atom-scale-label" class="slider-value">Scale: {atom_scale:.2f}</span></div>
                        <div class="slider-group"><label for="bond-thickness-slider" class="slider-label">Bond Thickness</label><input type="range" id="bond-thickness-slider" min="0.02" max="0.2" step="0.01" value="{bond_thickness}"><span id="bond-thickness-label" class="slider-value">Value: {bond_thickness:.2f}</span></div>
                        <div class="slider-group"><label for="bond-cutoff-slider" class="slider-label">Bond Cutoff Scale</label><input type="range" id="bond-cutoff-slider" min="0.8" max="1.5" step="0.01" value="{bond_cutoff_scale}"><span id="bond-cutoff-label" class="slider-value">Scale: {bond_cutoff_scale:.2f}</span></div>
                        <div class="slider-group"><label for="animation-speed-slider" class="slider-label">Animation Speed (FPS)</label><input type="range" id="animation-speed-slider" min="1" max="180" step="1" value="{animation_fps}"><span id="animation-speed-label" class="slider-value">FPS: {animation_fps}</span></div>
                        <div class="slider-group"><label for="bg-color-picker" class="slider-label">Background Color</label><input type="color" id="bg-color-picker" value="{bg_color}"></div>
                    </div>
                </details>
            </div>
            <div id="canvas-container"></div>
            <div id="top-right-controls">
                <button id="save-gif-btn">
                    <i class="fa fa-download"></i>
                    <span>Save as GIF</span>
                </button>
            </div>
        </div>

        <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/gif.js/0.2.0/gif.js"></script>

        <script>
            // --- DATA FROM PYTHON ---
            const initialAtomsData = {json_initial_atoms};
            const modeVectors = {json_mode_vectors};
            const frequencies = {json_freqs};
            const isImaginary = {json_is_imaginary};
            const nSteps = {n_frames};
            const atomicColorsHex = {json_atomic_colors};
            const covalentRadii = {json_covalent_radii};

            // --- GLOBAL STATE ---
            let scene, camera, renderer, controls;
            let atomGroup, bondGroup, cellGroup, vectorGroup;
            let animationInterval = null, isPlaying = false;
            let currentFrameIndex = 0;
            let currentModeIndex = {initial_mode_index};
            let atomInfo = {{}};
            let currentAmplitude = {amplitude};
            let currentAtomScale = {atom_scale};
            let currentBondCutoffFactor = {bond_cutoff_scale};
            let currentBondThickness = {bond_thickness};
            let currentAnimationSpeed = 1000 / {animation_fps};
            let currentStyle = 'lambert';

            // --- UI ELEMENTS ---
            const modeSelect = document.getElementById('mode-select');
            const frameSlider = document.getElementById('frame-slider');
            const copyXyzBtn = document.getElementById('copy-xyz-btn');
            const saveGifBtn = document.getElementById('save-gif-btn');

            // --- INITIALIZATION ---
            function init() {{
                buildAtomInfo();
                const container = document.getElementById('canvas-container');
                scene = new THREE.Scene();
                scene.background = new THREE.Color("{bg_color}");
                camera = new THREE.PerspectiveCamera(75, container.clientWidth / container.clientHeight, 0.1, 1000);
                camera.position.z = 20;
                renderer = new THREE.WebGLRenderer({{ antialias: true, preserveDrawingBuffer: true }});
                renderer.setSize(container.clientWidth, container.clientHeight);
                renderer.setPixelRatio(window.devicePixelRatio);
                container.appendChild(renderer.domElement);
                const ambientLight = new THREE.AmbientLight(0xcccccc, 0.8);
                scene.add(ambientLight);
                const directionalLight = new THREE.DirectionalLight(0xffffff, 0.6);
                directionalLight.position.set(5, 10, 7.5);
                scene.add(directionalLight);
                atomGroup = new THREE.Group();
                bondGroup = new THREE.Group();
                cellGroup = new THREE.Group();
                vectorGroup = new THREE.Group();
                scene.add(atomGroup, bondGroup, cellGroup, vectorGroup);
                setupControls(document.getElementById('rotation-mode-select').value);
                setupUI();
                switchMode(currentModeIndex);
                window.addEventListener('resize', onWindowResize);
                animate();
            }}

            function buildAtomInfo() {{
                const symbols = new Set(initialAtomsData.symbols);
                symbols.add('default');
                for (const symbol of symbols) {{
                    const color = atomicColorsHex[symbol] || atomicColorsHex['default'];
                    const radius = covalentRadii[symbol] || covalentRadii['default'];
                    atomInfo[symbol] = {{ color: new THREE.Color(color).getHex(), radius: radius }};
                }}
            }}

            function setupControls(mode) {{
                if (controls) controls.dispose();
                const container = document.getElementById('canvas-container');
                if (mode === 'trackball') {{
                    controls = new THREE.TrackballControls(camera, container);
                    controls.rotateSpeed = 5.0;
                }} else {{ // 'orbit'
                    controls = new THREE.OrbitControls(camera, container);
                    controls.enableDamping = true;
                }}
            }}

            // --- UI SETUP & EVENT LISTENERS ---
            function setupUI() {{
                // Sidebar Toggle Button Logic
                const toggleBtn = document.getElementById('toggle-sidebar-btn');
                const sidebar = document.getElementById('sidebar');
                toggleBtn.addEventListener('click', () => {{
                    sidebar.classList.toggle('collapsed');
                    toggleBtn.classList.toggle('open');
                    // Resize canvas after transition finishes
                    setTimeout(() => {{ onWindowResize(); }}, 310);
                }});

                frequencies.forEach((freq, i) => {{
                    const option = document.createElement('option');
                    option.value = i;
                    const indicator = isImaginary[i] ? '🔴' : '🟢';
                    const modeStr = `Mode ${{String(i + 1).padEnd(3, ' ')}} :`;
                    option.textContent = `${{modeStr}} ${{indicator}} ${{freq}} cm⁻¹`;
                    modeSelect.appendChild(option);
                }});

                // 드롭다운 메뉴의 초기 선택값 설정
                modeSelect.value = currentModeIndex;

                modeSelect.addEventListener('change', (e) => switchMode(parseInt(e.target.value)));
                document.getElementById('rotation-mode-select').addEventListener('change', (e) => setupControls(e.target.value));
                document.getElementById('style-select').addEventListener('change', (e) => {{ currentStyle = e.target.value; updateScene(currentFrameIndex, false); }});
                document.getElementById('bg-color-picker').addEventListener('input', (e) => scene.background.set(e.target.value));
                document.getElementById('bonds-toggle').addEventListener('change', (e) => {{ bondGroup.visible = e.target.checked; }});
                document.getElementById('cell-toggle').addEventListener('change', (e) => {{ cellGroup.visible = e.target.checked; }});
                document.getElementById('vectors-toggle').addEventListener('change', (e) => {{ vectorGroup.visible = e.target.checked; }});
                document.getElementById('amplitude-slider').addEventListener('input', (e) => {{
                    currentAmplitude = parseFloat(e.target.value);
                    document.getElementById('amplitude-label').textContent = `Scale: ${{currentAmplitude.toFixed(2)}}`;
                    updateScene(currentFrameIndex, false);
                    drawModeVectors();
                }});
                document.getElementById('atom-scale-slider').addEventListener('input', (e) => {{ currentAtomScale = parseFloat(e.target.value); document.getElementById('atom-scale-label').textContent = `Scale: ${{currentAtomScale.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                document.getElementById('bond-thickness-slider').addEventListener('input', (e) => {{ currentBondThickness = parseFloat(e.target.value); document.getElementById('bond-thickness-label').textContent = `Value: ${{currentBondThickness.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                document.getElementById('bond-cutoff-slider').addEventListener('input', (e) => {{ currentBondCutoffFactor = parseFloat(e.target.value); document.getElementById('bond-cutoff-label').textContent = `Scale: ${{currentBondCutoffFactor.toFixed(2)}}`; updateScene(currentFrameIndex, false); }});
                frameSlider.addEventListener('input', (e) => updateScene(parseInt(e.target.value)));
                document.getElementById('play-pause-btn').addEventListener('click', togglePlay);
                copyXyzBtn.addEventListener('click', copyTrajectoryToClipboard);
                saveGifBtn.addEventListener('click', saveAsGif);
                document.getElementById('animation-speed-slider').addEventListener('input', (e) => {{
                    const fps = parseInt(e.target.value);
                    document.getElementById('animation-speed-label').textContent = `FPS: ${{fps}}`;
                    currentAnimationSpeed = 1000 / fps;
                    if (isPlaying) {{
                        clearInterval(animationInterval);
                        animationInterval = setInterval(playNextFrame, currentAnimationSpeed);
                    }}
                }});
            }}

            function onWindowResize() {{
                const container = document.getElementById('canvas-container');
                if (container.clientWidth > 0 && container.clientHeight > 0) {{
                    camera.aspect = container.clientWidth / container.clientHeight;
                    camera.updateProjectionMatrix();
                    renderer.setSize(container.clientWidth, container.clientHeight);
                }}
            }}

            function animate() {{
                requestAnimationFrame(animate);
                controls.update();
                renderer.render(scene, camera);
            }}

            // --- The rest of the JavaScript functions (switchMode, updateScene, createAtom, etc.) remain unchanged ---
            // (The full script from the original prompt is included here for completeness)
            function switchMode(modeIdx) {{
                if (isPlaying) togglePlay();
                currentModeIndex = modeIdx;
                modeSelect.value = currentModeIndex;
                updateScene(0);
                drawModeVectors();
            }}

            function updateScene(frameIdx, updateSlider=true) {{
                currentFrameIndex = frameIdx;
                clearGroup(atomGroup); clearGroup(bondGroup); clearGroup(cellGroup);
                const time_step = (2 * Math.PI / nSteps) * frameIdx;
                const displacement_factor = currentAmplitude * Math.sin(time_step);
                const {{ symbols, positions: initialPositions, cell }} = initialAtomsData;
                const currentModeVector = modeVectors[currentModeIndex];
                const positions = initialPositions.map((initialPos, i) => {{
                    const displacement = new THREE.Vector3(...currentModeVector[i]).multiplyScalar(displacement_factor);
                    return new THREE.Vector3(...initialPos).add(displacement);
                }});
                positions.forEach((pos, i) => atomGroup.add(createAtom(pos, symbols[i])));
                for (let i = 0; i < positions.length; i++) {{
                    for (let j = i + 1; j < positions.length; j++) {{
                        const r_i = (atomInfo[symbols[i]] || atomInfo['default']).radius;
                        const r_j = (atomInfo[symbols[j]] || atomInfo['default']).radius;
                        const cutoff = (r_i + r_j) * currentBondCutoffFactor;
                        if (positions[i].distanceTo(positions[j]) < cutoff) {{
                            createBond(positions[i], positions[j], symbols[i], symbols[j]);
                        }}
                    }}
                }}
                if (cell && cell.flat().some(v => v !== 0)) drawCell(cell);
                bondGroup.visible = document.getElementById('bonds-toggle').checked;
                cellGroup.visible = document.getElementById('cell-toggle').checked;
                if (updateSlider) frameSlider.value = currentFrameIndex;
                document.getElementById('frame-label').textContent = `Frame: ${{currentFrameIndex}} / ${{nSteps - 1}}`;
            }}
            function getMaterial(color, style) {{
                switch(style) {{
                    case 'glossy': return new THREE.MeshPhongMaterial({{ color: color, shininess: 100, specular: 0x222222 }});
                    case 'metallic': return new THREE.MeshStandardMaterial({{ color: color, metalness: 0.3, roughness: 0.4 }});
                    case 'toon': return new THREE.MeshToonMaterial({{ color: color }});
                    case 'neon': return new THREE.MeshBasicMaterial({{ color: new THREE.Color(color).multiplyScalar(1.5) }});
                    case 'lambert':
                    default: return new THREE.MeshLambertMaterial({{ color: color }});
                }}
            }}
            function createAtom(pos, symbol) {{
                const info = atomInfo[symbol] || atomInfo['default'];
                const scaledRadius = info.radius * currentAtomScale;
                if (currentStyle === '2d') {{
                    const canvas = document.createElement('canvas');
                    const context = canvas.getContext('2d');
                    canvas.width = 128; canvas.height = 128;
                    context.beginPath(); context.arc(64, 64, 60, 0, 2 * Math.PI, false); context.fillStyle = 'black'; context.fill();
                    context.beginPath(); context.arc(64, 64, 56, 0, 2 * Math.PI, false); context.fillStyle = 'white'; context.fill();
                    context.font = 'bold 50px Arial'; context.fillStyle = 'black'; context.textAlign = 'center'; context.textBaseline = 'middle';
                    context.fillText(symbol, 64, 64);
                    const texture = new THREE.CanvasTexture(canvas);
                    const spriteMaterial = new THREE.SpriteMaterial({{ map: texture, transparent: true, depthTest: false }});
                    const sprite = new THREE.Sprite(spriteMaterial);
                    sprite.position.copy(pos);
                    sprite.scale.set(scaledRadius * 2.5, scaledRadius * 2.5, 1);
                    return sprite;
                }}
                const geometry = new THREE.SphereGeometry(scaledRadius, 32, 32);
                const material = getMaterial(info.color, currentStyle);
                const sphere = new THREE.Mesh(geometry, material);
                sphere.position.copy(pos);
                if (currentStyle === 'toon') {{
                    const outlineGeometry = geometry.clone();
                    outlineGeometry.scale(1.1, 1.1, 1.1);
                    const outlineMaterial = new THREE.MeshBasicMaterial({{ color: 0x000000, side: THREE.BackSide }});
                    const outline = new THREE.Mesh(outlineGeometry, outlineMaterial);
                    sphere.add(outline);
                }}
                return sphere;
            }}
            function createBond(p1, p2, sym1, sym2) {{
                if (p1.distanceTo(p2) <= 0.1) return;
                switch (currentStyle) {{
                    case '2d': {{{{
                        const dir = p2.clone().sub(p1).normalize();
                        const geometry = new THREE.CylinderGeometry(currentBondThickness, currentBondThickness, p1.distanceTo(p2), 8, 1);
                        const material = new THREE.MeshBasicMaterial({{ color: 0x666666 }});
                        const bond = new THREE.Mesh(geometry, material);
                        bond.position.copy(p1).add(p2).multiplyScalar(0.5);
                        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
                        bondGroup.add(bond);
                        break;
                    }}}}
                    case 'toon': {{{{
                        const path = new THREE.Line3(p1, p2);
                        const distance = path.distance();
                        const geometry = new THREE.CylinderGeometry(currentBondThickness, currentBondThickness, distance, 8);
                        const material = new THREE.MeshToonMaterial({{ color: 0x000000 }});
                        const bond = new THREE.Mesh(geometry, material);
                        bond.position.lerpVectors(p1, p2, 0.5);
                        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), path.delta(new THREE.Vector3()).normalize());
                        bondGroup.add(bond);
                        break;
                    }}}}
                    default: {{{{
                        const info1 = atomInfo[sym1] || atomInfo.default;
                        const info2 = atomInfo[sym2] || atomInfo.default;
                        const midPoint = p1.clone().add(p2).multiplyScalar(0.5);
                        const bond1 = createHalfBond(p1, midPoint, info1.color);
                        const bond2 = createHalfBond(midPoint, p2, info2.color);
                        if (bond1) bondGroup.add(bond1);
                        if (bond2) bondGroup.add(bond2);
                        break;
                    }}}}
                }}
            }}
            function createHalfBond(start, end, color) {{
                if (start.distanceTo(end) <= 0) return null;
                const path = new THREE.LineCurve3(start, end);
                const geometry = new THREE.TubeGeometry(path, 1, currentBondThickness, 8, false);
                const material = (currentStyle === 'neon') ? new THREE.MeshBasicMaterial({{ color: 0xffffff }}) : getMaterial(color, currentStyle);
                return new THREE.Mesh(geometry, material);
            }}
            function drawCell(cell) {{
                const [v1, v2, v3] = cell.map(c => new THREE.Vector3(...c));
                const points = [ new THREE.Vector3(0,0,0), v1, v2, v3, v1.clone().add(v2), v1.clone().add(v3), v2.clone().add(v3), v1.clone().add(v2).add(v3) ];
                const edges = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]];
                const lineMaterial = new THREE.LineBasicMaterial({{ color: 0xAAAAAA }});
                edges.forEach(edge => {{
                    const geometry = new THREE.BufferGeometry().setFromPoints([points[edge[0]], points[edge[1]]]);
                    cellGroup.add(new THREE.Line(geometry, lineMaterial));
                }});
            }}
            function drawModeVectors() {{
                clearGroup(vectorGroup);
                const initialPositions = initialAtomsData.positions.map(p => new THREE.Vector3(...p));
                const currentModeVec = modeVectors[currentModeIndex];
                initialPositions.forEach((pos, i) => {{
                    const dir = new THREE.Vector3(...currentModeVec[i]);
                    const baseLength = dir.length();
                    if (baseLength < 0.01) return;
                    dir.normalize();
                    const scaledLength = baseLength * currentAmplitude * 3.0;
                    const headLength = Math.max(scaledLength * 0.25, 0.15);
                    const headWidth = Math.max(scaledLength * 0.15, 0.1);
                    const forwardArrow = new THREE.ArrowHelper(dir.clone(), pos, scaledLength, 0x0096FF, headLength, headWidth);
                    vectorGroup.add(forwardArrow);
                    const backwardArrow = new THREE.ArrowHelper(dir.clone().negate(), pos, scaledLength, 0xff0000, headLength, headWidth);
                    vectorGroup.add(backwardArrow);
                }});
                vectorGroup.visible = document.getElementById('vectors-toggle').checked;
            }}
            function clearGroup(group) {{
                while (group.children.length > 0) {{
                    const child = group.children[0]; group.remove(child);
                    if(child.geometry) child.geometry.dispose();
                    if(child.material) {{
                        if (child.material.map) child.material.map.dispose();
                        if (Array.isArray(child.material)) {{
                            child.material.forEach(m => m.dispose());
                        }} else {{
                            child.material.dispose();
                        }}
                    }}
                }}
            }}
            function playNextFrame() {{
                let nextFrame = (currentFrameIndex + 1) % nSteps;
                updateScene(nextFrame);
            }}
            function togglePlay() {{
                isPlaying = !isPlaying;
                const playIcon = document.getElementById('play-icon');
                const pauseIcon = document.getElementById('pause-icon');
                if (isPlaying) {{
                    playIcon.style.display = 'none'; pauseIcon.style.display = 'block';
                    animationInterval = setInterval(playNextFrame, currentAnimationSpeed);
                }} else {{
                    playIcon.style.display = 'block'; pauseIcon.style.display = 'none';
                    clearInterval(animationInterval);
                }}
            }}
            async function copyTrajectoryToClipboard() {{
                const btn = document.getElementById('copy-xyz-btn');
                const originalText = btn.textContent;
                let xyzString = '';
                try {{
                    const nAtoms = initialAtomsData.symbols.length;
                    for (let i = 0; i < nSteps; i++) {{
                        const t = (2 * Math.PI / nSteps) * i;
                        const df = currentAmplitude * Math.sin(t);
                        xyzString += `${{nAtoms}}\\nFrame ${{i}}\\n`;
                        initialAtomsData.symbols.forEach((s, j) => {{
                            const p = initialAtomsData.positions[j];
                            const v = modeVectors[currentModeIndex][j];
                            xyzString += `${{s.padEnd(4)}} ${{(p[0] + v[0] * df).toFixed(8).padStart(12)}} ${{(p[1] + v[1] * df).toFixed(8).padStart(12)}} ${{(p[2] + v[2] * df).toFixed(8).padStart(12)}}\\n`;
                        }});
                    }}
                    await navigator.clipboard.writeText(xyzString);
                    btn.textContent = 'Copied!';
                }} catch (err) {{
                    console.error('Failed to copy trajectory: ', err);
                    btn.textContent = 'Copy Failed!';
                }} finally {{
                    setTimeout(() => {{ btn.textContent = originalText; }}, 2000);
                }}
            }}
            async function saveAsGif() {{
                const btn = document.getElementById('save-gif-btn');
                btn.classList.add('saving');
                btn.querySelector('span').textContent = 'Recording...';
                btn.disabled = true;
                let workerUrl;
                try {{
                    const workerScriptResponse = await fetch('https://cdnjs.cloudflare.com/ajax/libs/gif.js/0.2.0/gif.worker.js');
                    if (!workerScriptResponse.ok) throw new Error('Failed to fetch worker script');
                    const workerScriptText = await workerScriptResponse.text();
                    const workerBlob = new Blob([workerScriptText], {{ type: 'application/javascript' }});
                    workerUrl = URL.createObjectURL(workerBlob);
                    const wasPlaying = isPlaying;
                    if (wasPlaying) togglePlay();
                    const gif = new GIF({{
                        workers: 2, quality: 10,
                        width: renderer.domElement.width, height: renderer.domElement.height,
                        workerScript: workerUrl, background: scene.background.getStyle(),
                    }});
                    gif.on('finished', function(blob) {{
                        const url = URL.createObjectURL(blob);
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = `mode_${{currentModeIndex + 1}}.gif`;
                        a.click();
                        URL.revokeObjectURL(url);
                        URL.revokeObjectURL(workerUrl);
                        btn.classList.remove('saving');
                        btn.querySelector('span').textContent = 'Save as GIF';
                        btn.disabled = false;
                        if (wasPlaying) togglePlay();
                    }});
                    const cameraState = camera.clone();
                    let frame = 0;
                    function recordFrame() {{
                        if (frame < nSteps) {{
                            updateScene(frame, true);
                            camera.position.copy(cameraState.position);
                            camera.quaternion.copy(cameraState.quaternion);
                            camera.zoom = cameraState.zoom;
                            camera.updateProjectionMatrix();
                            renderer.render(scene, camera);
                            gif.addFrame(renderer.domElement, {{ copy: true, delay: currentAnimationSpeed }});
                            frame++;
                            setTimeout(recordFrame, 10);
                        }} else {{
                            gif.render();
                        }}
                    }}
                    recordFrame();
                }} catch (error) {{
                    console.error("Failed to create GIF:", error);
                    btn.classList.remove('saving');
                    btn.querySelector('span').textContent = 'Error!';
                    btn.disabled = false;
                    if(workerUrl) URL.revokeObjectURL(workerUrl);
                }}
            }}

            // --- START ---
            init();
        </script>
    </body>
    </html>
    """

    if write_html:
        assert isinstance(write_html, str), "write_html must be a string"
        with open(write_html, "w", encoding="utf-8") as f:
            f.write(html_template)

    return HTML(html_template)
