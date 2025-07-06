"""
Molecule service - handles molecule rendering and processing
This module focuses on the 2D rendering functionality
"""

import logging
import io
import base64
import math
from typing import Optional, Dict, Any, Tuple, List
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from models import Molecule
from routes.api import molecules  # Using the in-memory storage

logger = logging.getLogger(__name__)

class MoleculeRenderingError(Exception):
    """Exception raised for errors in the molecule rendering process."""
    pass

def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Parse a SMILES string into an RDKit molecule object
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        RDKit Mol object or None if parsing failed
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return None
        return mol
    except Exception as e:
        logger.error(f"Error parsing SMILES {smiles}: {str(e)}")
        return None

def parse_mol_file(mol_file_content: str) -> Optional[Chem.Mol]:
    """
    Parse a MOL file content into an RDKit molecule object
    
    Args:
        mol_file_content: String content of a MOL file
        
    Returns:
        RDKit Mol object or None if parsing failed
    """
    try:
        mol = Chem.MolFromMolBlock(mol_file_content)
        if mol is None:
            logger.warning("Failed to parse MOL file")
            return None
        return mol
    except Exception as e:
        logger.error(f"Error parsing MOL file: {str(e)}")
        return None

def prepare_molecule_for_rendering(mol: Chem.Mol) -> Chem.Mol:
    """
    Prepare the molecule for rendering by adding hydrogens and computing 2D coordinates
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Prepared RDKit Mol object
    """
    try:
        # Make a copy of the molecule to avoid modifying the original
        mol_copy = Chem.Mol(mol)
        
        # Clean up the molecule
        Chem.SanitizeMol(mol_copy)
        
        # Compute 2D coordinates if they don't exist
        if not mol_copy.GetNumConformers():
            AllChem.Compute2DCoords(mol_copy)
        
        return mol_copy
    except Exception as e:
        logger.error(f"Error preparing molecule for rendering: {str(e)}")
        raise MoleculeRenderingError(f"Failed to prepare molecule: {str(e)}")

def generate_molecule_svg(mol: Chem.Mol, width: int = 300, height: int = 200) -> str:
    """
    Generate an SVG representation of a molecule
    
    Args:
        mol: RDKit Mol object
        width: Width of the SVG image
        height: Height of the SVG image
        
    Returns:
        SVG string representation of the molecule
    """
    try:
        # Prepare molecule for rendering
        mol = prepare_molecule_for_rendering(mol)
        
        # Get 2D coordinates
        conf = mol.GetConformer()
        
        # Get atom and bond information
        atoms = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                'symbol': symbol,
                'x': pos.x, 
                'y': pos.y,
                'z': pos.z,
                'idx': atom.GetIdx(),
                'charge': atom.GetFormalCharge(),
                'is_aromatic': atom.GetIsAromatic()
            })
        
        bonds = []
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            is_aromatic = bond.GetIsAromatic()
            bonds.append({
                'begin': begin_idx,
                'end': end_idx,
                'type': str(bond_type),
                'is_aromatic': is_aromatic
            })
        
        # Scale coordinates to fit in the given dimensions
        x_coords = [a['x'] for a in atoms]
        y_coords = [a['y'] for a in atoms]
        
        if x_coords and y_coords:  # Check if lists are not empty
            min_x, max_x = min(x_coords), max(x_coords)
            min_y, max_y = min(y_coords), max(y_coords)
            
            padding = 20  # Padding in pixels
            available_width = width - 2 * padding
            available_height = height - 2 * padding
            
            x_range = max_x - min_x if max_x > min_x else 1.0
            y_range = max_y - min_y if max_y > min_y else 1.0
            
            scale_x = available_width / x_range
            scale_y = available_height / y_range
            
            # Use the smaller scale to maintain aspect ratio
            scale = min(scale_x, scale_y)
            
            # Transform coordinates
            for atom in atoms:
                atom['svg_x'] = padding + (atom['x'] - min_x) * scale
                atom['svg_y'] = height - (padding + (atom['y'] - min_y) * scale)  # SVG y-axis is flipped
        else:
            for atom in atoms:
                atom['svg_x'] = width / 2
                atom['svg_y'] = height / 2
        
        # Generate SVG
        svg_parts = [
            f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg" version="1.1">',
            f'<rect width="100%" height="100%" fill="#222" opacity="0" />'
        ]
        
        # Draw bonds
        for bond in bonds:
            begin = atoms[bond['begin']]
            end = atoms[bond['end']]
            x1, y1 = begin['svg_x'], begin['svg_y']
            x2, y2 = end['svg_x'], end['svg_y']
            
            # Determine styling based on bond type
            stroke_width = 2
            if bond['is_aromatic']:
                stroke = "#3ea5bf"  # Cyan-ish for aromatic bonds
            else:
                stroke = "#88b"  # Bluish for normal bonds
            
            # Draw the bond line
            svg_parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{stroke_width}" />')
            
            # Add second line for double or triple bonds
            if bond['type'] in ['DOUBLE', 'TRIPLE']:
                # Calculate perpendicular offset
                dx = x2 - x1
                dy = y2 - y1
                length = math.sqrt(dx*dx + dy*dy)
                
                if length > 0:
                    dx, dy = dx/length, dy/length
                    offset = 3  # Pixel offset for parallel line
                    
                    # Perpendicular vector
                    px, py = -dy*offset, dx*offset
                    
                    # Draw parallel line
                    svg_parts.append(f'<line x1="{x1+px}" y1="{y1+py}" x2="{x2+px}" y2="{y2+py}" stroke="{stroke}" stroke-width="{stroke_width}" />')
                    
                    # Add a third line for triple bonds
                    if bond['type'] == 'TRIPLE':
                        svg_parts.append(f'<line x1="{x1-px}" y1="{y1-py}" x2="{x2-px}" y2="{y2-py}" stroke="{stroke}" stroke-width="{stroke_width}" />')
        
        # Draw atoms
        for atom in atoms:
            symbol = atom['symbol']
            x, y = atom['svg_x'], atom['svg_y']
            
            # Skip drawing carbon atoms in most cases
            if symbol == 'C' and any(bond['begin'] == atom['idx'] or bond['end'] == atom['idx'] for bond in bonds):
                # Only draw carbons that have bonds
                continue
                
            # Determine text color based on atom type
            atom_colors = {
                'C': '#88a',  # Carbon - blue-gray
                'N': '#77d',  # Nitrogen - blue
                'O': '#e55',  # Oxygen - red
                'S': '#db0',  # Sulfur - yellow
                'P': '#e70',  # Phosphorus - orange
                'F': '#6cd',  # Fluorine - light blue
                'Cl': '#0c0', # Chlorine - green
                'Br': '#a60', # Bromine - brown
                'I': '#909',  # Iodine - purple
                'H': '#aaa'   # Hydrogen - light gray
            }
            fill_color = atom_colors.get(symbol, '#ddd')  # Default to light gray
            
            # Add charge if non-zero
            display_text = symbol
            if atom['charge'] > 0:
                display_text += f'<tspan dy="-5" font-size="10">+{atom["charge"] if atom["charge"] > 1 else ""}</tspan>'
            elif atom['charge'] < 0:
                display_text += f'<tspan dy="-5" font-size="10">-{abs(atom["charge"]) if atom["charge"] < -1 else ""}</tspan>'
            
            # Draw the atom symbol with a slight background for readability
            svg_parts.append(f'<circle cx="{x}" cy="{y}" r="10" fill="#222" fill-opacity="0.7" />')
            svg_parts.append(f'<text x="{x}" y="{y}" text-anchor="middle" dominant-baseline="middle" fill="{fill_color}" font-family="Arial" font-size="14">{display_text}</text>')
        
        svg_parts.append('</svg>')
        return ''.join(svg_parts)
    except Exception as e:
        logger.error(f"Error generating SVG: {str(e)}")
        # Return a minimal fallback SVG with error message
        return f'''<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
            <rect width="100%" height="100%" fill="#343a40" />
            <text x="50%" y="50%" text-anchor="middle" dominant-baseline="middle" fill="#f8f9fa" font-family="Arial" font-size="14">
                Rendering Error
            </text>
        </svg>'''

def process_molecule(molecule: Molecule) -> Molecule:
    """
    Process a molecule to add 2D visualization
    
    Args:
        molecule: Molecule object
        
    Returns:
        Updated Molecule object with SVG visualization
    """
    try:
        mol = None
        
        # Try to parse from SMILES if available
        if molecule.smiles:
            mol = parse_smiles(molecule.smiles)
        
        # If SMILES parsing failed or SMILES not available, try MOL file
        if mol is None and molecule.mol_file:
            mol = parse_mol_file(molecule.mol_file)
        
        # If both methods failed, we can't visualize
        if mol is None:
            logger.warning(f"Could not create RDKit molecule for {molecule.name or molecule.id or 'unknown'}")
            molecule.svg = None
            return molecule
        
        # Generate SVG
        svg = generate_molecule_svg(mol)
        molecule.svg = svg
        
        # Initialize properties dictionary if it doesn't exist
        if not hasattr(molecule, 'properties') or not molecule.properties:
            molecule.properties = {}
        
        # Add basic properties
        try:
            # Calculate simple properties without using descriptors that might be missing
            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
            
            # Create a simple formula string
            formula_parts = []
            for symbol in sorted(atom_counts.keys()):
                count = atom_counts[symbol]
                if count == 1:
                    formula_parts.append(symbol)
                else:
                    formula_parts.append(f"{symbol}{count}")
            
            molecule.properties['formula'] = ''.join(formula_parts)
            molecule.properties['num_atoms'] = mol.GetNumAtoms()
            molecule.properties['num_bonds'] = mol.GetNumBonds()
            
            # Count rings using a simple method (may not be as accurate as RDKit's method)
            # Euler formula for planar graphs: V - E + F = 2
            # For molecules: F (faces) = R (rings) + 1 (outer face)
            # So: R = E - V + 1
            # Where V = number of atoms, E = number of bonds
            try:
                num_rings = max(0, molecule.properties['num_bonds'] - molecule.properties['num_atoms'] + 1)
                molecule.properties['num_rings'] = num_rings
            except:
                # If there's any issue, just skip this property
                pass
            
            # Try to import the descriptors module for more advanced properties
            try:
                from rdkit.Chem import Descriptors, rdMolDescriptors
                
                # Add more detailed properties if the modules are available
                try:
                    molecule.properties['formula'] = rdMolDescriptors.CalcMolFormula(mol)
                except:
                    pass
                
                try:
                    molecule.properties['weight'] = Descriptors.MolWt(mol)
                except:
                    pass
                
                try:
                    molecule.properties['logp'] = Descriptors.MolLogP(mol)
                except:
                    pass
                
                try:
                    molecule.properties['num_rings'] = rdMolDescriptors.CalcNumRings(mol)
                except:
                    pass
                    
            except ImportError:
                logger.warning("RDKit descriptors module not available, using simplified properties")
        
        except Exception as e:
            logger.warning(f"Could not calculate all properties: {str(e)}")
        
        return molecule
    except Exception as e:
        logger.error(f"Error processing molecule: {str(e)}")
        # Don't lose the original molecule data, just fail the visualization
        molecule.svg = None
        return molecule

def generate_conformers(mol: Chem.Mol, n_conformers: int = 10, max_attempts: int = 1000) -> Tuple[bool, Chem.Mol]:
    """
    Generate multiple conformers for a molecule
    
    Args:
        mol: RDKit Mol object
        n_conformers: Number of conformers to generate
        max_attempts: Maximum number of attempts for conformer generation
        
    Returns:
        Tuple of (success, molecule with conformers)
    """
    try:
        # Make a copy to avoid modifying the original
        mol_copy = Chem.Mol(mol)
        
        # Add hydrogens 
        mol_copy = Chem.AddHs(mol_copy)
        
        # Generate conformers
        result = AllChem.EmbedMultipleConfs(
            mol_copy, 
            numConfs=n_conformers,
            maxAttempts=max_attempts,
            randomSeed=42,  # For reproducibility
            pruneRmsThresh=0.5  # Remove similar conformers
        )
        
        if result == -1 or mol_copy.GetNumConformers() == 0:
            logger.warning("Failed to generate any conformers")
            return False, mol
            
        # Optimize conformers
        for conf_id in range(mol_copy.GetNumConformers()):
            try:
                AllChem.UFFOptimizeMolecule(mol_copy, confId=conf_id)
            except Exception as e:
                logger.warning(f"Failed to optimize conformer {conf_id}: {str(e)}")
        
        return True, mol_copy
    except Exception as e:
        logger.error(f"Error in conformer generation: {str(e)}")
        return False, mol


class MoleculeService:
    def __init__(self):
        self.molecules = molecules  # Using the existing in-memory storage

    def get_molecule(self, mol_id: str) -> Optional[Molecule]:
        """Retrieve a molecule by its ID"""
        return self.molecules.get(mol_id)

    def save_molecule(self, molecule: Molecule) -> Molecule:
        """Save a molecule to storage"""
        self.molecules[molecule.id] = molecule
        return molecule

    def get_all_molecules(self) -> Dict[str, Molecule]:
        """Get all molecules"""
        return self.molecules

    def delete_molecule(self, mol_id: str) -> bool:
        """Delete a molecule by its ID"""
        if mol_id in self.molecules:
            del self.molecules[mol_id]
            return True
        return False

    def update_molecule(self, mol_id: str, updates: Dict) -> Optional[Molecule]:
        """Update molecule properties"""
        molecule = self.get_molecule(mol_id)
        if molecule:
            for key, value in updates.items():
                if hasattr(molecule, key):
                    setattr(molecule, key, value)
            return molecule
        return None
