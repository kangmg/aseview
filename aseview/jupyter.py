import ipywidgets as widgets
from IPython.display import display, HTML
import json
from typing import Union, Dict, Any
from ase import Atoms
from .viewers import BaseViewer, MolecularData


class JupyterViewer(widgets.VBox):
    """Jupyter widget wrapper for molecular viewers."""
    
    def __init__(self, viewer: BaseViewer):
        """
        Initialize the Jupyter widget.
        
        Args:
            viewer: A BaseViewer instance (MolecularViewer, NormalViewer, or OverlayViewer)
        """
        self.viewer = viewer
        
        # Create the main HTML output widget
        self.html_output = widgets.HTML(
            value=self.viewer.get_html(),
            placeholder='Molecular Viewer',
            description='',
        )
        
        # Create control widgets
        self.sidebar_toggle = widgets.ToggleButton(
            value=True,
            description='Toggle Sidebar',
            disabled=False,
            button_style='info',
            tooltip='Toggle sidebar visibility',
        )
        
        self.save_button = widgets.Button(
            description='Save HTML',
            disabled=False,
            button_style='success',
            tooltip='Save as HTML file',
        )
        
        self.filename_input = widgets.Text(
            value='molecule_viewer.html',
            placeholder='Filename',
            description='Filename:',
            disabled=False
        )
        
        # Set up event handlers
        self.sidebar_toggle.observe(self._on_sidebar_toggle, names='value')
        self.save_button.on_click(self._on_save_click)
        
        # Create controls box
        controls = widgets.HBox([self.sidebar_toggle, self.filename_input, self.save_button])
        
        # Initialize the parent VBox
        super().__init__(children=[controls, self.html_output])
        
        # Display the widget
        display(self)
    
    def _on_sidebar_toggle(self, change):
        """Handle sidebar toggle button click."""
        # In a full implementation, this would send a message to the JavaScript
        # to toggle the sidebar. For now, we'll just update the HTML.
        pass
    
    def _on_save_click(self, button):
        """Handle save button click."""
        filename = self.filename_input.value
        html_content = self.viewer.get_html()
        
        # In a Jupyter notebook, we can't directly save files to the filesystem
        # for security reasons. Instead, we'll show how to save it.
        js_code = f"""
        var element = document.createElement('a');
        element.setAttribute('href', 'data:text/html;charset=utf-8,' + encodeURIComponent(`{html_content}`));
        element.setAttribute('download', '{filename}');
        element.style.display = 'none';
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
        """
        
        display(HTML(f'<script>{js_code}</script>'))
        print(f"Saving {filename}...")


def view_molecule(data: Union[Atoms, Dict[str, Any], str], viewer_type: str = "molecular", **kwargs):
    """
    Create and display a molecular viewer in Jupyter.
    
    Args:
        data: ASE Atoms object, dictionary with molecular data, or path to JSON file
        viewer_type: Type of viewer ("molecular", "normal", or "overlay")
        **kwargs: Additional settings for the viewer
        
    Returns:
        JupyterViewer widget
    """
    # Import here to avoid circular imports
    from .viewers import MolecularViewer, NormalViewer, OverlayViewer
    
    # Create the appropriate viewer
    if viewer_type == "molecular":
        viewer = MolecularViewer(data, **kwargs)
    elif viewer_type == "normal":
        viewer = NormalViewer(data, **kwargs)
    elif viewer_type == "overlay":
        viewer = OverlayViewer(data, **kwargs)
    else:
        raise ValueError(f"Unknown viewer type: {viewer_type}")
    
    # Create and return the Jupyter widget
    return JupyterViewer(viewer)