import ipywidgets as widgets
from IPython.display import display, HTML
import json
import base64
from typing import Union, Dict, Any, List
from ase import Atoms
from .wrapper import BaseViewer, MolecularData


class JupyterViewer(widgets.VBox):
    """Jupyter widget wrapper for molecular viewers."""

    def __init__(self, viewer: BaseViewer, height: int = 600):
        """
        Initialize the Jupyter widget.

        Args:
            viewer: A BaseViewer instance (MolecularViewer, NormalViewer, or OverlayViewer)
            height: Height of the viewer in pixels
        """
        self.viewer = viewer
        self._height = height

        # Create the main HTML output widget with iframe
        html_content = self.viewer.get_html()

        # Escape HTML for srcdoc attribute to handle large molecules
        escaped_html = (html_content
            .replace('&', '&amp;')
            .replace('"', '&quot;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
        )

        iframe_html = f'''
<div style="width: 100%; height: {height}px; position: relative;">
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
        self.html_output = widgets.HTML(value=iframe_html)

        # Create control widgets
        self.save_button = widgets.Button(
            description='Save HTML',
            disabled=False,
            button_style='success',
            tooltip='Save as HTML file',
            icon='download'
        )

        self.filename_input = widgets.Text(
            value='molecule_viewer.html',
            placeholder='Filename',
            description='Filename:',
            disabled=False,
            layout=widgets.Layout(width='250px')
        )

        # Set up event handlers
        self.save_button.on_click(self._on_save_click)

        # Create controls box
        controls = widgets.HBox([self.filename_input, self.save_button])

        # Initialize the parent VBox
        super().__init__(children=[controls, self.html_output])

    def _on_save_click(self, button):
        """Handle save button click."""
        filename = self.filename_input.value
        html_content = self.viewer.get_html()

        # Use base64 to safely encode the HTML content
        html_b64 = base64.b64encode(html_content.encode('utf-8')).decode('ascii')

        # Create download link using data URI
        js_code = f"""
        (function() {{
            var link = document.createElement('a');
            link.href = 'data:text/html;base64,{html_b64}';
            link.download = '{filename}';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }})();
        """

        display(HTML(f'<script>{js_code}</script>'))
        print(f"Downloading {filename}...")


def view_molecule(
    data: Union[Atoms, Dict[str, Any], str, List],
    viewer_type: str = "molecular",
    height: int = 600,
    **kwargs
):
    """
    Create and display a molecular viewer in Jupyter.

    Args:
        data: ASE Atoms object, dictionary with molecular data, path to JSON file,
              or list of any of these
        viewer_type: Type of viewer ("molecular", "normal", or "overlay")
        height: Height of the viewer in pixels
        **kwargs: Additional settings for the viewer

    Returns:
        JupyterViewer widget
    """
    # Import here to avoid circular imports
    from .wrapper import MolecularViewer, NormalViewer, OverlayViewer

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
    return JupyterViewer(viewer, height=height)
