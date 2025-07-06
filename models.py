"""
This module contains the data models for the application.
Since we're primarily using in-memory storage for this prototype,
these classes serve as data structures rather than database models.
"""

from typing import Dict, List, Optional
from datetime import datetime

class Molecule:
    """Represents a molecule structure with its properties"""
    
    def __init__(
        self, 
        mol_id: Optional[str] = None, 
        smiles: Optional[str] = None, 
        name: Optional[str] = None, 
        mol_file: Optional[str] = None
    ):
        self.id = mol_id
        self.smiles = smiles  # SMILES string representation
        self.name = name
        self.mol_file = mol_file  # Could be an SDF, PDB, or other format
        self.svg = None  # SVG representation for 2D visualization
        self.conformers: List[Dict] = []
        self.properties: Dict = {}
        self.pdb_data: Optional[str] = None
        self.trajectory_frames: List[Dict] = []
        self.structural_analysis: Dict = {}
        
    def to_dict(self) -> Dict:
        """Convert molecule to dictionary for JSON serialization"""
        result = {
            'id': self.id,
            'name': self.name,
            'smiles': self.smiles,
            'properties': self.properties,
            'conformers_count': len(self.conformers),
            'pdb_data': self.pdb_data,
            'trajectory_frames': self.trajectory_frames,
            'structural_analysis': self.structural_analysis
        }
        if self.svg:
            result['svg'] = self.svg
        return result


class BindingSite:
    """Represents an enzyme binding site"""
    
    def __init__(
        self, 
        site_id: Optional[str] = None, 
        protein_id: Optional[str] = None, 
        name: Optional[str] = None, 
        coordinates: Optional[List[float]] = None
    ):
        self.id = site_id
        self.protein_id = protein_id
        self.name = name
        self.coordinates = coordinates or []  # 3D coordinates of the binding site
        self.properties: Dict = {}
    
    def to_dict(self) -> Dict:
        """Convert binding site to dictionary for JSON serialization"""
        return {
            'id': self.id,
            'protein_id': self.protein_id,
            'name': self.name,
            'coordinates': self.coordinates,
            'properties': self.properties
        }


class DockingResult:
    """Represents the result of a docking prediction"""
    
    def __init__(
        self, 
        molecule: Optional[Molecule] = None, 
        binding_site: Optional[BindingSite] = None
    ):
        self.id: Optional[str] = None  # Unique identifier
        self.molecule = molecule
        self.binding_site = binding_site
        self.score: Optional[float] = None
        self.conformer_id: Optional[str] = None
        self.pose_data: Optional[List[float]] = None  # 3D coordinates of the docked pose
        self.interactions: List[Dict] = []  # Key interactions between molecule and binding site
        self.binding_affinity: Optional[float] = None  # Binding affinity in kcal/mol
        self.toxicity: Dict = {}  # Toxicity predictions
        self.timestamp: Optional[datetime] = datetime.now()
    
    def to_dict(self) -> Dict:
        """Convert docking result to dictionary for JSON serialization"""
        return {
            'id': self.id,
            'molecule': self.molecule.to_dict() if self.molecule else None,
            'binding_site': self.binding_site.to_dict() if self.binding_site else None,
            'score': self.score,
            'conformer_id': self.conformer_id,
            'interactions': self.interactions,
            'binding_affinity': self.binding_affinity,
            'toxicity': self.toxicity,
            'pose_data': self.pose_data,
            'timestamp': self.timestamp.isoformat() if self.timestamp else None
        }
