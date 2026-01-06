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
            import base64
            
            html = self._generate_html()
            
            # Use base64 data URI for better Jupyter compatibility
            html_bytes = html.encode('utf-8')
            html_b64 = base64.b64encode(html_bytes).decode('ascii')
            
            # Create iframe with data URI
            iframe_html = f'''
<div style="width: {width if isinstance(width, str) else str(width) + 'px'}; height: {height}px; position: relative;">
    <iframe 
        src="data:text/html;base64,{html_b64}" 
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
            "bondThreshold": 1.7,
            "bondThickness": 0.1,
            "atomSize": 0.4,
            "animationSpeed": 30,
            "forceScale": 1.0,
            "backgroundColor": "#1f2937",
            "style": "Cartoon",
            "showCell": True,
            "showBond": True,
            "showShadow": False,
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
        function initData() {{
            console.log('Initializing molecular data...');
            
            // Apply settings from Python
            const pythonSettings = {settings_json};
            Object.assign(settings, pythonSettings);
            console.log('Applied settings:', settings);
            
            // Wait for renderer to be initialized
            function trySetData() {{
                if (typeof renderer !== 'undefined' && renderer !== null) {{
                    console.log('Renderer ready, setting data');
                    if (typeof setMolecularData === 'function') {{
                        console.log('setMolecularData found, calling with data');
                        setMolecularData({js_data});
                    }} else {{
                        console.error('setMolecularData function not found!');
                    }}
                }} else {{
                    console.log('Renderer not ready, waiting...');
                    setTimeout(trySetData, 50);
                }}
            }}
            
            trySetData();
        }}
        
        // Try immediate execution (for iframes where DOM is already ready)
        if (document.readyState === 'loading') {{
            console.log('DOM still loading, waiting for DOMContentLoaded');
            document.addEventListener('DOMContentLoaded', initData);
        }} else {{
            console.log('DOM already ready, initializing immediately');
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
    
    def __init__(self, data: Union[Atoms, Dict[str, Any], str, List], **kwargs):
        super().__init__(data)
        self.settings = {
            "bondThreshold": 1.7,
            "bondThickness": 0.1,
            "atomSize": 0.4,
            "animationSpeed": 30,
            "forceScale": 1.0,
            "backgroundColor": "#1f2937",
            "style": "default",
            "showCell": True,
            "showBond": True,
            "showShadow": False,
            "displacementAmplitude": 1.0,
            "showModeVector": False,
            **kwargs
        }
    
    def _generate_html(self) -> str:
        """Generate the HTML content for the normal mode viewer."""
        js_data = self._get_js_data()
        
        # Load normal_viewer template (now with THREE.js support)
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
            # Replace external script tag with inline version
            external_tag = '<script src="/static/js/styles.js"></script>'
            if external_tag in html:
                html = html.replace(external_tag, f'<script>\n{styles_js}\n</script>')
        
        # Inject molecular data and settings
        settings_json = json.dumps(self.settings)
        script = f"""
<script>
    (function() {{
        function initData() {{
            console.log('Initializing normal mode data...');
            
            // Apply settings from Python
            const pythonSettings = {settings_json};
            Object.assign(settings, pythonSettings);
            console.log('Applied settings:', settings);
            
            // Wait for renderer to be initialized
            function trySetData() {{
                if (typeof renderer !== 'undefined' && renderer !== null) {{
                    console.log('Renderer ready, setting data');
                    if (typeof setMolecularData === 'function') {{
                        setMolecularData({js_data});
                    }} else {{
                        console.error('setMolecularData function not found!');
                    }}
                }} else {{
                    console.log('Renderer not ready, waiting...');
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
            "bondThreshold": 1.7,
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
        
        # Inject molecular data and settings before closing script tag
        settings_json = json.dumps(self.settings)
        
        # Insert initialization code before the last closing </script> tag
        init_script = f"""
// Initialize data from Python
(function() {{
    function initData() {{
        console.log('Initializing overlay data...');
        
        // Apply settings from Python
        const pythonSettings = {settings_json};
        Object.assign(settings, pythonSettings);
        console.log('Applied settings:', settings);
        
        // Wait for renderer to be initialized
        function trySetData() {{
            if (typeof renderer !== 'undefined' && renderer !== null) {{
                console.log('Renderer ready, setting data');
                if (typeof setMolecularData === 'function') {{
                    setMolecularData({js_data});
                }} else {{
                    console.error('setMolecularData function not found!');
                }}
            }} else {{
                console.log('Renderer not ready, waiting...');
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
        # Replace the last occurrence of </script></body>
        last_script_close = html.rfind('</script>')
        if last_script_close != -1:
            html = html[:last_script_close] + init_script + html[last_script_close+9:]
        
        return html
    
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
