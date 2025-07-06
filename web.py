"""
Web routes for the Neural Chem AI application
"""

import logging
from flask import Blueprint, render_template, redirect, url_for, request, flash
from models import Molecule
import services.molecule_service as molecule_service
import uuid

# Create a blueprint for the web routes
from flask import Blueprint

web_bp = Blueprint('web', __name__)
logger = logging.getLogger(__name__)

# Import in-memory storage from API
from routes.api import molecules, binding_sites, docking_results
from models import BindingSite

@web_bp.route('/')
def index():
    """Render the homepage"""
    try:
        return render_template('index.html')
    except Exception as e:
        logger.error(f"Error rendering index: {str(e)}", exc_info=True)
        return "An error occurred while loading the page", 500

@web_bp.route('/molecules')
def list_molecules():
    """List all molecules"""
    try:
        return render_template('molecules.html', molecules=list(molecules.values()))
    except Exception as e:
        logger.error(f"Error listing molecules: {str(e)}", exc_info=True)
        return "An error occurred while loading the molecules", 500

@web_bp.route('/molecules/add', methods=['GET', 'POST'])
def add_molecule():
    """Add a new molecule"""
    try:
        if request.method == 'POST':
            # Generate a unique ID
            mol_id = str(uuid.uuid4())
            
            # Create the molecule object
            molecule = Molecule(
                mol_id=mol_id,
                smiles=request.form.get('smiles'),
                name=request.form.get('name'),
                mol_file=request.form.get('mol_file')
            )
            
            # Process the molecule (generate SVG, etc.)
            molecule = molecule_service.process_molecule(molecule)
            
            # Store the molecule
            molecules[mol_id] = molecule
            
            flash('Molecule added successfully', 'success')
            return redirect(url_for('web.list_molecules'))
        
        return render_template('add_molecule.html')
    except Exception as e:
        logger.error(f"Error adding molecule: {str(e)}", exc_info=True)
        flash(f"Error adding molecule: {str(e)}", 'danger')
        return redirect(url_for('web.list_molecules'))

@web_bp.route('/molecules/<molecule_id>')
def view_molecule(molecule_id):
    """View a specific molecule"""
    try:
        molecule = molecules.get(molecule_id)
        if not molecule:
            flash('Molecule not found', 'danger')
            return redirect(url_for('web.list_molecules'))
            
        # Add ML predictions to molecule
        import services.ml_service as ml_service
        molecule = ml_service.enrich_molecule_with_predictions(molecule)
            
        # Get toxicity prediction
        import services.binding_site_service as binding_site_service
        toxicity = binding_site_service.calculate_toxicity(molecule)
            
        return render_template('molecule_detail.html', molecule=molecule, toxicity=toxicity)
    except Exception as e:
        logger.error(f"Error viewing molecule: {str(e)}", exc_info=True)
        flash(f"Error viewing molecule: {str(e)}", 'danger')
        return redirect(url_for('web.list_molecules'))

@web_bp.route('/binding-sites')
def list_binding_sites():
    """List all binding sites"""
    try:
        return render_template('binding_sites.html', binding_sites=list(binding_sites.values()))
    except Exception as e:
        logger.error(f"Error listing binding sites: {str(e)}", exc_info=True)
        return "An error occurred while loading the binding sites", 500

@web_bp.route('/binding-sites/add', methods=['GET', 'POST'])
def add_binding_site():
    """Add a new binding site"""
    try:
        if request.method == 'POST':
            # Generate a unique ID
            site_id = str(uuid.uuid4())
            
            # Parse coordinates
            coords_str = request.form.get('coordinates', '0,0,0')
            try:
                coordinates = [float(x.strip()) for x in coords_str.split(',') if x.strip()]
                if len(coordinates) < 3:
                    coordinates.extend([0.0] * (3 - len(coordinates)))
            except:
                coordinates = [0.0, 0.0, 0.0]
            
            # Create the binding site object
            binding_site = BindingSite(
                site_id=site_id,
                protein_id=request.form.get('protein_id'),
                name=request.form.get('name'),
                coordinates=coordinates
            )
            
            # Process the binding site
            import services.binding_site_service as binding_site_service
            binding_site = binding_site_service.process_binding_site(binding_site)
            
            # Store the binding site
            binding_sites[site_id] = binding_site
            
            flash('Binding site added successfully', 'success')
            return redirect(url_for('web.list_binding_sites'))
        
        return render_template('add_binding_site.html')
    except Exception as e:
        logger.error(f"Error adding binding site: {str(e)}", exc_info=True)
        flash(f"Error adding binding site: {str(e)}", 'danger')
        return redirect(url_for('web.list_binding_sites'))

@web_bp.route('/binding-sites/<binding_site_id>')
def view_binding_site(binding_site_id):
    """View a specific binding site"""
    try:
        binding_site = binding_sites.get(binding_site_id)
        if not binding_site:
            flash('Binding site not found', 'danger')
            return redirect(url_for('web.list_binding_sites'))
            
        return render_template('binding_site_detail.html', binding_site=binding_site)
    except Exception as e:
        logger.error(f"Error viewing binding site: {str(e)}", exc_info=True)
        flash(f"Error viewing binding site: {str(e)}", 'danger')
        return redirect(url_for('web.list_binding_sites'))

@web_bp.route('/docking')
def docking_page():
    """Show docking interface"""
    try:
        return render_template('docking.html', 
                              molecules=list(molecules.values()), 
                              binding_sites=list(binding_sites.values()),
                              docking_results=list(docking_results.values()))
    except Exception as e:
        logger.error(f"Error loading docking page: {str(e)}", exc_info=True)
        return "An error occurred while loading the docking page", 500

@web_bp.route('/docking/predict', methods=['GET', 'POST'])
def perform_docking():
    """Perform molecular docking prediction"""
    try:
        if request.method == 'POST':
            molecule_id = request.form.get('molecule_id')
            binding_site_id = request.form.get('binding_site_id')
            
            if not molecule_id or not binding_site_id:
                flash('Both molecule and binding site are required', 'danger')
                return redirect(url_for('web.docking_page'))
                
            # Retrieve molecule and binding site
            molecule = molecules.get(molecule_id)
            binding_site = binding_sites.get(binding_site_id)
            
            if not molecule or not binding_site:
                flash('Molecule or binding site not found', 'danger')
                return redirect(url_for('web.docking_page'))
                
            # Predict docking
            import services.binding_site_service as binding_site_service
            import services.ml_service as ml_service
            
            # Get docking result
            docking_result = binding_site_service.predict_docking(molecule, binding_site)
            
            # Add binding affinity prediction
            docking_result.binding_affinity = ml_service.binding_affinity_prediction(molecule, binding_site)
            
            # Generate a unique ID for the result
            result_id = str(uuid.uuid4())
            docking_result.id = result_id  # Set the ID in the result object
            docking_results[result_id] = docking_result
            
            return redirect(url_for('web.view_docking_result', result_id=result_id))
            
        return redirect(url_for('web.docking_page'))
    except Exception as e:
        logger.error(f"Error performing docking: {str(e)}", exc_info=True)
        flash(f"Error performing docking: {str(e)}", 'danger')
        return redirect(url_for('web.docking_page'))

@web_bp.route('/docking/results/<result_id>')
def view_docking_result(result_id):
    """View a specific docking result"""
    try:
        docking_result = docking_results.get(result_id)
        if not docking_result:
            flash('Docking result not found', 'danger')
            return redirect(url_for('web.docking_page'))
            
        return render_template('docking_result.html', result=docking_result)
    except Exception as e:
        logger.error(f"Error viewing docking result: {str(e)}", exc_info=True)
        flash(f"Error viewing docking result: {str(e)}", 'danger')
        return redirect(url_for('web.docking_page'))
