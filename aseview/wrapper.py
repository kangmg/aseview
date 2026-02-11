"""
Wrapper module for molecular visualization
"""
import json
import os
from typing import Dict, Any, List, Union
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

        # Add forces if present (from calculator or arrays)
        try:
            if atoms.calc is not None:
                forces = atoms.get_forces()
                data["forces"] = forces.tolist()
        except Exception:
            pass

        # Fallback to arrays if not from calculator
        if "forces" not in data and hasattr(atoms, 'arrays') and 'forces' in atoms.arrays:
            data["forces"] = atoms.arrays['forces'].tolist()

        # Add energy if available (from calculator or info dict)
        try:
            if atoms.calc is not None:
                data["energy"] = float(atoms.get_potential_energy())
        except Exception:
            pass

        # Check info dict for energy (common in trajectory files)
        if "energy" not in data and hasattr(atoms, 'info'):
            for key in ['energy', 'Energy', 'E', 'total_energy']:
                if key in atoms.info:
                    data["energy"] = float(atoms.info[key])
                    break

        # Add name if available in info dict
        if hasattr(atoms, 'info') and 'name' in atoms.info:
            data["name"] = str(atoms.info['name'])

        # Add charges if available (from arrays or info dict)
        charges = None
        if hasattr(atoms, 'arrays') and 'charges' in atoms.arrays:
            charges = atoms.arrays['charges']
        elif hasattr(atoms, 'info') and 'charges' in atoms.info:
            charges = atoms.info['charges']

        if charges is not None:
            try:
                import numpy as np
                charges_array = np.asarray(charges, dtype=float).flatten()
                if len(charges_array) == len(symbols):
                    data["charges"] = charges_array.tolist()
            except (TypeError, ValueError):
                pass

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
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str, List]):
        """
        Initialize the viewer with molecular data.
        
        Args:
            data: Can be an ASE Atoms object, a dictionary with molecular data,
                  a path to a JSON file containing molecular data, or a list of any of these.
        """
        self.data = self._process_data(data)
    
    def _process_data(self, data):
        """Process input data into a standardized format."""
        if isinstance(data, str):
            # Assume it's a path to a JSON file
            with open(data, 'r') as f:
                return json.load(f)
        elif isinstance(data, Atoms):
            return MolecularData.from_atoms(data)
        elif isinstance(data, list):
            # If it's a list, process each element
            processed_list = []
            for item in data:
                if isinstance(item, Atoms):
                    processed_list.append(MolecularData.from_atoms(item))
                elif isinstance(item, dict):
                    processed_list.append(item)
                else:
                    raise ValueError(f"Unsupported data type in list: {type(item)}")
            return processed_list
        elif isinstance(data, dict):
            return data
        else:
            raise ValueError(f"Unsupported data type: {type(data)}")
    
    def _get_js_data(self) -> str:
        """Convert molecular data to JavaScript-compatible JSON string."""
        # Ensure data is a list for consistency
        if not isinstance(self.data, list):
            data_to_convert = [self.data]
        else:
            data_to_convert = self.data
            
        return json.dumps(data_to_convert)
    
    def show(self, width='100%', height=600):
        """Display the viewer in a Jupyter notebook."""
        try:
            from IPython.display import display, HTML

            html = self._generate_html()

            # Escape HTML for srcdoc attribute (handles quotes and special chars)
            # This avoids data: URI size limits that cause failures with large molecules
            escaped_html = (html
                .replace('&', '&amp;')
                .replace('"', '&quot;')
                .replace('<', '&lt;')
                .replace('>', '&gt;')
            )

            # Use srcdoc instead of data: URI to avoid size limits
            # srcdoc supports much larger content than data: URIs
            width_style = width if isinstance(width, str) else f'{width}px'
            iframe_html = f'''
<div style="width: {width_style}; height: {height}px; position: relative;">
    <iframe
        srcdoc="{escaped_html}"
        width="100%"
        height="100%"
        style="border: 1px solid #ddd; border-radius: 4px; display: block;"
        sandbox="allow-scripts allow-same-origin"
        scrolling="no">
    </iframe>
</div>
'''
            display(HTML(iframe_html))
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
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str, List], **kwargs):
        super().__init__(data)
        self.settings = {
            "bondThreshold": 1.2,  # Scale factor for covalent radii sum
            "bondThickness": 0.1,
            "atomSize": 0.4,
            "animationSpeed": 30,
            "forceScale": 1.0,
            "backgroundColor": "#1f2937",
            "style": "Cartoon",
            "showCell": True,
            "showBond": True,
            "showShadow": False,
            "showShading": True,
            "showEnergyPlot": False,
            "showForces": False,
            "colorBy": "Element",  # "Element" or "Charge"
            "normalizeCharges": False,  # Normalize charges to -1 to 1 range
            "chargeColormap": "coolwarm",  # Colormap for charge visualization
            "viewMode": "Perspective",
            "rotationMode": "TrackBall",
            "selectionMode": "Lasso",
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the molecular viewer."""
        js_data = self._get_js_data()

        # Load the HTML template (prefer packaged template, fallback to static copy)
        try:
            template_candidates = [
                os.path.join(os.path.dirname(__file__), "..", "static", "molecular_viewer.html"),
                os.path.join(os.path.dirname(__file__), "templates", "molecular_viewer.html"),
            ]

            html = None
            for candidate in template_candidates:
                if os.path.exists(candidate):
                    with open(candidate, 'r', encoding='utf-8') as f:
                        html = f.read()
                    break

            if html is None:
                return self._generate_simple_html()
        except FileNotFoundError:
            return self._generate_simple_html()

        # Inline all JavaScript dependencies for Jupyter compatibility
        vendor_dir = os.path.join(os.path.dirname(__file__), "static", "js", "vendor")
        
        # Inline Three.js
        three_path = os.path.join(vendor_dir, "three.min.js")
        if os.path.exists(three_path):
            with open(three_path, 'r', encoding='utf-8') as f:
                three_js = f.read()
            html = html.replace(
                '<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>',
                f'<script>{three_js}</script>'
            )
        
        # Inline OrbitControls
        orbit_path = os.path.join(vendor_dir, "OrbitControls.js")
        if os.path.exists(orbit_path):
            with open(orbit_path, 'r', encoding='utf-8') as f:
                orbit_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"></script>',
                f'<script>{orbit_js}</script>'
            )
        
        # Inline TrackballControls
        trackball_path = os.path.join(vendor_dir, "TrackballControls.js")
        if os.path.exists(trackball_path):
            with open(trackball_path, 'r', encoding='utf-8') as f:
                trackball_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>',
                f'<script>{trackball_js}</script>'
            )
        
        # Inline styles.js
        styles_path = os.path.join(os.path.dirname(__file__), "static", "js", "styles.js")
        if os.path.exists(styles_path):
            with open(styles_path, 'r', encoding='utf-8') as style_file:
                styles_js = style_file.read()

            external_tag = '<script src="/static/js/styles.js"></script>'
            if external_tag in html:
                inline_styles_tag = f"<script>\n{styles_js}\n</script>"
                html = html.replace(external_tag, inline_styles_tag)
        
        # Inline gifshot.js for GIF export
        gifshot_path = os.path.join(vendor_dir, "gifshot.js")
        if os.path.exists(gifshot_path):
            with open(gifshot_path, 'r', encoding='utf-8') as f:
                gifshot_js = f.read()
            # Add gifshot.js before closing </head> tag
            html = html.replace('</head>', f'<script>{gifshot_js}</script>\n</head>')

        # Inject molecular data
        if '{{molecular_data}}' in html:
            html = html.replace('{{molecular_data}}', js_data)
            return html

        # If template does not expose placeholder, call setMolecularData at load
        # Use both DOMContentLoaded and immediate execution for iframe compatibility
        settings_json = json.dumps(self.settings)
        script = f"""
<script>
    (function() {{
        var MAX_RETRIES = 100;  // 5 seconds max (100 * 50ms)
        var retryCount = 0;

        function initData() {{
            console.log('Initializing molecular data...');

            // Apply settings from Python (ensure settings exists)
            if (typeof settings === 'undefined') {{
                window.settings = {{}};
            }}
            const pythonSettings = {settings_json};
            Object.assign(settings, pythonSettings);

            // Wait for renderer to be initialized with timeout
            function trySetData() {{
                retryCount++;
                if (typeof renderer !== 'undefined' && renderer !== null) {{
                    if (typeof setMolecularData === 'function') {{
                        setMolecularData({js_data});
                    }} else {{
                        console.error('setMolecularData function not found!');
                    }}
                }} else if (retryCount < MAX_RETRIES) {{
                    setTimeout(trySetData, 50);
                }} else {{
                    console.error('Renderer initialization timeout after ' + (MAX_RETRIES * 50) + 'ms');
                }}
            }}

            trySetData();
        }}

        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', initData);
        }} else {{
            initData();
        }}
    }})();
</script>
</body>
"""
        return html.replace('</body>', script)
    
    def _generate_simple_html(self) -> str:
        """Generate a simple HTML fallback."""
        return f"""
<!DOCTYPE html>
<html>
<head>
    <title>Molecular Viewer</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ font-family: Arial, sans-serif; background: #111827; color: #f9fafb; display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        .container {{ text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Molecular Viewer</h1>
        <p>Molecular data loaded successfully</p>
        <div>
            <p>Number of molecules: {len(self.data) if isinstance(self.data, list) else 1}</p>
            <p>Atoms in first molecule: {len(self.data[0]['symbols']) if (isinstance(self.data, list) and self.data and 'symbols' in self.data[0]) else (len(self.data['symbols']) if (isinstance(self.data, dict) and 'symbols' in self.data) else 0)}</p>
        </div>
    </div>
    <script>
        function setMolecularData(data) {{
            console.log('Molecular data set:', data);
        }}
        setMolecularData({self._get_js_data()});
    </script>
</body>
</html>
        """


class NormalViewer(BaseViewer):
    """Normal mode viewer for molecular vibrations."""

    def __init__(
        self,
        atoms: Union[Atoms, Dict[str, Any]],
        vibrations=None,
        mode_vectors: List = None,
        frequencies: List = None,
        n_frames: int = 30,
        **kwargs
    ):
        """
        Initialize NormalViewer for visualizing vibrational modes.

        Args:
            atoms: ASE Atoms object of equilibrium structure
            vibrations: ASE Vibrations or VibrationsData object (optional)
            mode_vectors: List of mode displacement vectors (optional, extracted from vibrations if not provided)
            frequencies: List of frequencies in cm^-1 (optional, extracted from vibrations if not provided)
            n_frames: Number of animation frames per cycle
            **kwargs: Additional viewer settings
        """
        # Process atoms data
        if isinstance(atoms, Atoms):
            self.atoms_data = MolecularData.from_atoms(atoms)
        elif isinstance(atoms, dict):
            self.atoms_data = atoms
        else:
            raise ValueError(f"Unsupported atoms type: {type(atoms)}")

        # Store as data for compatibility
        self.data = self.atoms_data

        # Extract vibration data
        self.mode_vectors = None
        self.frequencies = None
        self.is_imaginary = None

        if vibrations is not None:
            self._extract_vibration_data(vibrations)
        elif mode_vectors is not None and frequencies is not None:
            self.mode_vectors = mode_vectors if isinstance(mode_vectors, list) else mode_vectors.tolist()
            self.frequencies, self.is_imaginary = self._format_frequencies(frequencies)

        self.n_frames = n_frames

        self.settings = {
            "bondThreshold": 1.2,
            "bondThickness": 0.1,
            "atomSize": 0.4,
            "animationSpeed": 30,
            "backgroundColor": "#1f2937",
            "style": "cartoon",
            "showCell": True,
            "showBond": True,
            "showShadow": False,
            "displacementAmplitude": 0.75,
            "showModeVector": False,
            "initialModeIndex": 0,
            "nFrames": n_frames,
            **kwargs
        }

    def _extract_vibration_data(self, vib):
        """Extract mode vectors and frequencies from ASE Vibrations object."""
        import numpy as np

        # Try to get VibrationsData
        try:
            from ase.vibrations import Vibrations
            from ase.vibrations.data import VibrationsData
        except ImportError:
            raise ImportError("ASE vibrations module is required")

        if hasattr(vib, 'get_vibrations'):
            vib_data = vib.get_vibrations()
        else:
            vib_data = vib

        # Get mode vectors (all atoms)
        modes = vib_data.get_modes(all_atoms=True)
        self.mode_vectors = modes.tolist()

        # Get frequencies and format them
        freqs = vib_data.get_frequencies()
        self.frequencies, self.is_imaginary = self._format_frequencies(freqs)

    def _format_frequencies(self, freqs):
        """Format frequencies for display, handling imaginary values."""
        import numpy as np

        formatted = []
        is_imaginary = []

        for f in freqs:
            if np.iscomplexobj(f) and f.imag != 0:
                formatted.append(f"{f.imag:.2f}i")
                is_imaginary.append(True)
            else:
                val = f.real if np.iscomplexobj(f) else float(f)
                formatted.append(f"{val:.2f}")
                is_imaginary.append(False)

        return formatted, is_imaginary

    def _get_js_data(self) -> str:
        """Get molecular data as JSON for JavaScript."""
        return json.dumps(self.atoms_data)

    def _get_vibration_data(self) -> dict:
        """Get vibration data for JavaScript."""
        return {
            "modeVectors": self.mode_vectors,
            "frequencies": self.frequencies,
            "isImaginary": self.is_imaginary,
            "nFrames": self.n_frames
        }

    def _generate_html(self) -> str:
        """Generate the HTML content for the normal mode viewer."""
        atoms_json = self._get_js_data()

        # Prepare vibration data
        if self.mode_vectors is not None:
            vib_data = self._get_vibration_data()
            vib_json = json.dumps(vib_data)
            has_vibrations = "true"
        else:
            vib_json = "null"
            has_vibrations = "false"

        # Load normal_viewer template
        try:
            template_path = os.path.join(os.path.dirname(__file__), "templates", "normal_viewer.html")
            with open(template_path, 'r', encoding='utf-8') as f:
                html = f.read()
        except FileNotFoundError:
            return self._generate_simple_html("Normal Mode Viewer")

        # Inline JavaScript dependencies
        vendor_dir = os.path.join(os.path.dirname(__file__), "static", "js", "vendor")

        # Inline Three.js
        three_path = os.path.join(vendor_dir, "three.min.js")
        if os.path.exists(three_path):
            with open(three_path, 'r', encoding='utf-8') as f:
                three_js = f.read()
            html = html.replace(
                '<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>',
                f'<script>{three_js}</script>'
            )

        # Inline OrbitControls
        orbit_path = os.path.join(vendor_dir, "OrbitControls.js")
        if os.path.exists(orbit_path):
            with open(orbit_path, 'r', encoding='utf-8') as f:
                orbit_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"></script>',
                f'<script>{orbit_js}</script>'
            )

        # Inline TrackballControls
        trackball_path = os.path.join(vendor_dir, "TrackballControls.js")
        if os.path.exists(trackball_path):
            with open(trackball_path, 'r', encoding='utf-8') as f:
                trackball_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>',
                f'<script>{trackball_js}</script>'
            )

        # Inline styles.js
        styles_path = os.path.join(os.path.dirname(__file__), "static", "js", "styles.js")
        if os.path.exists(styles_path):
            with open(styles_path, 'r', encoding='utf-8') as f:
                styles_js = f.read()
            external_tag = '<script src="/static/js/styles.js"></script>'
            if external_tag in html:
                html = html.replace(external_tag, f'<script>\n{styles_js}\n</script>')

        # Inline gifshot.js for GIF export
        gifshot_path = os.path.join(vendor_dir, "gifshot.js")
        if os.path.exists(gifshot_path):
            with open(gifshot_path, 'r', encoding='utf-8') as f:
                gifshot_js = f.read()
            html = html.replace('</head>', f'<script>{gifshot_js}</script>\n</head>')

        # Inject data and settings
        settings_json = json.dumps(self.settings)
        script = f"""
<script>
    (function() {{
        var MAX_RETRIES = 100;
        var retryCount = 0;

        // Inject data from Python
        window.equilibriumAtoms = {atoms_json};
        window.vibrationData = {vib_json};
        window.hasVibrations = {has_vibrations};

        function initData() {{
            if (typeof settings === 'undefined') {{
                window.settings = {{}};
            }}
            const pythonSettings = {settings_json};
            Object.assign(settings, pythonSettings);

            function trySetData() {{
                retryCount++;
                if (typeof renderer !== 'undefined' && renderer !== null) {{
                    if (typeof initNormalModeViewer === 'function') {{
                        initNormalModeViewer(window.equilibriumAtoms, window.vibrationData);
                    }} else if (typeof setMolecularData === 'function') {{
                        // Fallback for templates without vibration support
                        setMolecularData([window.equilibriumAtoms]);
                    }}
                }} else if (retryCount < MAX_RETRIES) {{
                    setTimeout(trySetData, 50);
                }}
            }}

            trySetData();
        }}

        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', initData);
        }} else {{
            initData();
        }}
    }})();
</script>
</body>
"""
        return html.replace('</body>', script)

    def _generate_simple_html(self, title: str) -> str:
        """Generate a simple HTML fallback."""
        return f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <meta charset="UTF-8">
    <style>
        body {{ font-family: Arial, sans-serif; background: #111827; color: #f9fafb;
                display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        .container {{ text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <p>Data loaded successfully</p>
    </div>
</body>
</html>
        """

    @classmethod
    def from_orca(cls, atoms, hess_file: str, skip_imaginary: bool = False, **kwargs):
        """
        Create NormalViewer from ORCA .hess file.

        Args:
            atoms: ASE Atoms object of equilibrium structure
            hess_file: Path to ORCA .hess file
            skip_imaginary: If True, skip translational/rotational modes (freq < 10 cm^-1). Default False.
            **kwargs: Additional viewer settings

        Returns:
            NormalViewer instance

        Example:
            >>> from ase.io import read
            >>> atoms = read("molecule.xyz")
            >>> viewer = NormalViewer.from_orca(atoms, "orca.hess")
            >>> viewer.show()
        """
        from .hessian_parsers import parse_orca_hess, reshape_modes_to_atoms, get_real_vibrations

        # Parse ORCA .hess file
        frequencies, normal_modes, n_atoms_hess = parse_orca_hess(hess_file)

        # Verify atom count matches
        if isinstance(atoms, Atoms):
            n_atoms = len(atoms)
        else:
            n_atoms = len(atoms.get('symbols', []))

        if n_atoms != n_atoms_hess:
            raise ValueError(
                f"Atom count mismatch: structure has {n_atoms} atoms, "
                f"but Hessian file has {n_atoms_hess} atoms"
            )

        # Filter out translations/rotations if requested
        if skip_imaginary:
            frequencies, normal_modes = get_real_vibrations(frequencies, normal_modes)

        # Reshape modes to per-atom format
        mode_vectors = reshape_modes_to_atoms(normal_modes, n_atoms)

        return cls(atoms, mode_vectors=mode_vectors, frequencies=frequencies.tolist(), **kwargs)


class OverlayViewer(BaseViewer):
    """Overlay viewer for comparing multiple molecules simultaneously."""
    
    def __init__(self, data: Union[Atoms, List[Atoms], Dict[str, Any], str, List], **kwargs):
        """
        Initialize the overlay viewer with one or more molecules.
        
        Args:
            data: Can be a single Atoms object, a list of Atoms objects,
                  or any format supported by BaseViewer
            **kwargs: Additional settings for visualization
        """
        super().__init__(data)
        
        # Ensure data is always a list for overlay viewer
        if not isinstance(self.data, list):
            self.data = [self.data]
        
        self.settings = {
            "bondThreshold": 1.2,  # Scale factor for covalent radii sum
            "bondThickness": 0.1,
            "atomSize": 0.4,
            "backgroundColor": "#1f2937",
            "style": "cartoon",
            "showCell": True,
            "showBond": True,
            "showShadow": False,
            "viewMode": "Perspective",
            "rotationMode": "TrackBall",
            "colorBy": "Atom",
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the overlay viewer."""
        js_data = self._get_js_data()
        
        # Load overlay_viewer template
        try:
            template_path = os.path.join(os.path.dirname(__file__), "templates", "overlay_viewer.html")
            with open(template_path, 'r', encoding='utf-8') as f:
                html = f.read()
        except FileNotFoundError:
            return self._generate_simple_html("Overlay Viewer")
        
        # Inline JavaScript dependencies
        vendor_dir = os.path.join(os.path.dirname(__file__), "static", "js", "vendor")
        
        # Inline Three.js
        three_path = os.path.join(vendor_dir, "three.min.js")
        if os.path.exists(three_path):
            with open(three_path, 'r', encoding='utf-8') as f:
                three_js = f.read()
            html = html.replace(
                '<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>',
                f'<script>{three_js}</script>'
            )
        
        # Inline OrbitControls
        orbit_path = os.path.join(vendor_dir, "OrbitControls.js")
        if os.path.exists(orbit_path):
            with open(orbit_path, 'r', encoding='utf-8') as f:
                orbit_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"></script>',
                f'<script>{orbit_js}</script>'
            )
        
        # Inline TrackballControls
        trackball_path = os.path.join(vendor_dir, "TrackballControls.js")
        if os.path.exists(trackball_path):
            with open(trackball_path, 'r', encoding='utf-8') as f:
                trackball_js = f.read()
            html = html.replace(
                '<script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>',
                f'<script>{trackball_js}</script>'
            )
        
        # Inline styles.js
        styles_path = os.path.join(os.path.dirname(__file__), "static", "js", "styles.js")
        if os.path.exists(styles_path):
            with open(styles_path, 'r', encoding='utf-8') as f:
                styles_js = f.read()
            # Replace external script tag with inline version
            external_tag = '<script src="/static/js/styles.js"></script>'
            if external_tag in html:
                html = html.replace(external_tag, f'<script>\n{styles_js}\n</script>')
        
        # Inject molecular data and settings
        settings_json = json.dumps(self.settings)
        script = f"""
<script>
    (function() {{
        var MAX_RETRIES = 100;
        var retryCount = 0;

        function initData() {{
            if (typeof settings === 'undefined') {{
                window.settings = {{}};
            }}
            const pythonSettings = {settings_json};
            Object.assign(settings, pythonSettings);

            function trySetData() {{
                retryCount++;
                if (typeof renderer !== 'undefined' && renderer !== null) {{
                    if (typeof setMolecularData === 'function') {{
                        setMolecularData({js_data});
                    }}
                }} else if (retryCount < MAX_RETRIES) {{
                    setTimeout(trySetData, 50);
                }}
            }}

            trySetData();
        }}

        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', initData);
        }} else {{
            initData();
        }}
    }})();
</script>
</body>
"""
        return html.replace('</body>', script)
    
    def _generate_simple_html(self, title: str) -> str:
        """Generate a simple HTML fallback."""
        num_molecules = len(self.data) if isinstance(self.data, list) else 1
        return f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <meta charset="UTF-8">
    <style>
        body {{ font-family: Arial, sans-serif; background: #111827; color: #f9fafb; 
                display: flex; justify-content: center; align-items: center; height: 100vh; margin: 0; }}
        .container {{ text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <p>Loaded {num_molecules} molecule(s) for overlay comparison</p>
    </div>
</body>
</html>
"""
