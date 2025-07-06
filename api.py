from flask import Blueprint, jsonify, request
from models import Molecule, BindingSite, DockingResult

# Create Blueprint
api_bp = Blueprint('api', __name__)

# In-memory storage (moved here to avoid circular imports)
molecules = {}
binding_sites = {}
docking_results = {}

@api_bp.route('/api/molecules', methods=['GET'])
def get_molecules():
    return jsonify([mol.to_dict() for mol in molecules.values()])

@api_bp.route('/api/molecules/<mol_id>', methods=['GET'])
def get_molecule(mol_id):
    molecule = molecules.get(mol_id)
    if molecule:
        return jsonify(molecule.to_dict())
    return jsonify({'error': 'Molecule not found'}), 404

@api_bp.route('/api/binding-sites', methods=['GET'])
def get_binding_sites():
    return jsonify([site.to_dict() for site in binding_sites.values()])

@api_bp.route('/api/binding-sites/<site_id>', methods=['GET'])
def get_binding_site(site_id):
    site = binding_sites.get(site_id)
    if site:
        return jsonify(site.to_dict())
    return jsonify({'error': 'Binding site not found'}), 404

@api_bp.route('/api/docking/<docking_id>', methods=['GET'])
def get_docking_result(docking_id):
    result = docking_results.get(docking_id)
    if result:
        return jsonify(result.to_dict())
    return jsonify({'error': 'Docking result not found'}), 404
