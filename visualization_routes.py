from flask import Blueprint, jsonify
from services.visualization_service import VisualizationService
from services.molecule_service import MoleculeService
from services.docking_service import DockingService
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

viz_bp = Blueprint('visualization', __name__)
viz_service = VisualizationService()
mol_service = MoleculeService()
dock_service = DockingService()

@viz_bp.route('/api/molecules/<mol_id>/3d', methods=['GET'])
def get_3d_structure(mol_id):
    try:
        molecule = mol_service.get_molecule(mol_id)
        if not molecule or not molecule.smiles:
            return jsonify({'error': 'Molecule not found'}), 404

        # Generate 3D structure
        mol = Chem.MolFromSmiles(molecule.smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Convert to PDB format
        pdb_data = Chem.MolToPDBBlock(mol)
        return jsonify({'pdb_data': pdb_data})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@viz_bp.route('/api/molecules/<mol_id>/analysis', methods=['GET'])
def analyze_molecule(mol_id):
    try:
        molecule = mol_service.get_molecule(mol_id)
        if not molecule or not molecule.smiles:
            return jsonify({'error': 'Molecule not found'}), 404

        mol = Chem.MolFromSmiles(molecule.smiles)
        analysis = {
            'molecular_weight': Descriptors.ExactMolWt(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'h_bond_donors': Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
            'rings': Descriptors.RingCount(mol)
        }
        return jsonify(analysis)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@viz_bp.route('/api/docking/<docking_id>/animation', methods=['GET'])
def get_docking_animation(docking_id):
    try:
        docking_result = dock_service.get_result(docking_id)
        if not docking_result:
            return jsonify({'error': 'Docking result not found'}), 404
            
        frames = dock_service.generate_animation_frames(docking_result)
        return jsonify({'frames': frames})
    except Exception as e:
        return jsonify({'error': str(e)}), 500