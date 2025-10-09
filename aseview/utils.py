"""
Utility functions for molecular data handling
"""
import json
from typing import Dict, Any, List, Union
from ase import Atoms


def load_ase_data(filepath: str) -> Union[Atoms, List[Atoms]]:
    """
    Load molecular data from various file formats using ASE.
    
    Args:
        filepath: Path to the file containing molecular data
        
    Returns:
        ASE Atoms object or list of ASE Atoms objects
    """
    try:
        from ase.io import read
        data = read(filepath, index=':')
        if len(data) == 1:
            return data[0]
        return data
    except Exception as e:
        raise ValueError(f"Could not load file {filepath}: {str(e)}")


def convert_to_json(data: Union[Atoms, List[Atoms], Dict, List[Dict]], filepath: str = None) -> str:
    """
    Convert molecular data to JSON format.
    
    Args:
        data: Molecular data in various formats
        filepath: Optional path to save the JSON file
        
    Returns:
        JSON string representation of the data
    """
    from .wrapper import MolecularData
    
    if isinstance(data, Atoms):
        json_data = MolecularData.from_atoms(data)
    elif isinstance(data, list) and all(isinstance(item, Atoms) for item in data):
        json_data = [MolecularData.from_atoms(item) for item in data]
    elif isinstance(data, dict) or (isinstance(data, list) and all(isinstance(item, dict) for item in data)):
        json_data = data
    else:
        raise ValueError("Unsupported data type for conversion")
    
    json_string = json.dumps(json_data, indent=2)
    
    if filepath:
        with open(filepath, 'w') as f:
            f.write(json_string)
    
    return json_string


def load_json_data(filepath: str) -> Union[Dict, List[Dict]]:
    """
    Load molecular data from a JSON file.
    
    Args:
        filepath: Path to the JSON file
        
    Returns:
        Dictionary or list of dictionaries containing molecular data
    """
    with open(filepath, 'r') as f:
        return json.load(f)