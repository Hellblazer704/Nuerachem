from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np
from typing import Dict, List, Optional
from models import Molecule, BindingSite, DockingResult

class VisualizationService:
    def generate_3d_structure(self, molecule: Molecule) -> str:
        mol = Chem.MolFromSmiles(molecule.smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToPDBBlock(mol)

    def analyze_structure(self, molecule: Molecule) -> Dict:
        mol = Chem.MolFromSmiles(molecule.smiles)
        return {
            'molecular_weight': Chem.Descriptors.ExactMolWt(mol),
            'rotatable_bonds': Chem.Descriptors.NumRotatableBonds(mol),
            'h_bond_donors': Chem.Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Chem.Descriptors.NumHAcceptors(mol),
            'rings': Chem.Descriptors.RingCount(mol)
        }

    def generate_docking_frames(self, docking_result: DockingResult) -> List[Dict]:
        start_coords = np.array(docking_result.molecule.conformers[0]['coordinates'])
        end_coords = np.array(docking_result.pose_data)
        
        frames = []
        for i in range(30):  # Generate 30 frames
            t = i / 29.0
            interpolated_coords = start_coords * (1-t) + end_coords * t
            
            frame = {
                'frame_id': i,
                'coordinates': interpolated_coords.tolist(),
                'progress': t,
                'score': docking_result.score * t
            }
            frames.append(frame)
            
        return frames