"""
Sample data initialization for the Neural Chem AI application.
This module contains functions to create sample molecules and binding sites
for testing and demonstration purposes.
"""

import logging
import uuid
from models import Molecule, BindingSite
from services import molecule_service

logger = logging.getLogger(__name__)

# Well-known SMILES strings for common molecules
SAMPLE_MOLECULES = [
    {
        "name": "Aspirin",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "description": "Common pain reliever and anti-inflammatory drug"
    },
    {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "description": "Stimulant found in coffee and tea"
    },
    {
        "name": "Ibuprofen",
        "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "description": "Non-steroidal anti-inflammatory drug (NSAID)"
    }
]

# Sample binding sites
SAMPLE_BINDING_SITES = [
    {
        "name": "Cyclooxygenase-1 (COX-1)",
        "protein_id": "P23219",
        "coordinates": [23.4, 45.6, 12.3],
        "description": "Enzyme involved in prostaglandin synthesis"
    },
    {
        "name": "Adenosine A2A Receptor",
        "protein_id": "P29274",
        "coordinates": [10.1, 22.5, 33.7],
        "description": "G protein-coupled receptor that binds adenosine"
    }
]

from models import Molecule, BindingSite, DockingResult
from storage import molecules, binding_sites, docking_results

def create_sample_molecules():
    molecules = {}
    
    # Add Aspirin as a sample molecule
    aspirin = Molecule(
        mol_id="mol1",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        name="Aspirin"
    )
    molecules[aspirin.id] = aspirin
    
    # Add Ibuprofen as another sample
    ibuprofen = Molecule(
        mol_id="mol2",
        smiles="CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        name="Ibuprofen"
    )
    molecules[ibuprofen.id] = ibuprofen
    
    return molecules

def create_sample_binding_sites():
    binding_sites = {}
    
    # Add sample binding site
    site = BindingSite(
        site_id="site1",
        protein_id="1ABC",
        name="Active Site",
        coordinates=[10.0, 15.0, 20.0]
    )
    binding_sites[site.id] = site
    
    return binding_sites

def create_sample_molecules():
    """Create sample molecules for demonstration"""
    molecules = {}
    
    for mol_data in SAMPLE_MOLECULES:
        try:
            # Generate a unique ID
            mol_id = str(uuid.uuid4())
            
            # Create the molecule object
            molecule = Molecule(
                mol_id=mol_id,
                smiles=mol_data["smiles"],
                name=mol_data["name"]
            )
            
            # Process the molecule (generate SVG, properties)
            molecule = molecule_service.process_molecule(molecule)
            
            # Add description to properties
            if not hasattr(molecule, 'properties') or not molecule.properties:
                molecule.properties = {}
            
            molecule.properties['description'] = mol_data.get("description", "")
            
            # Store the molecule
            molecules[mol_id] = molecule
            logger.info(f"Created sample molecule: {molecule.name}")
            
        except Exception as e:
            logger.error(f"Error creating sample molecule {mol_data['name']}: {str(e)}")
    
    return molecules

def create_sample_binding_sites():
    """Create sample binding sites for demonstration"""
    binding_sites = {}
    
    for site_data in SAMPLE_BINDING_SITES:
        try:
            # Generate a unique ID
            site_id = str(uuid.uuid4())
            
            # Create the binding site object
            binding_site = BindingSite(
                site_id=site_id,
                protein_id=site_data["protein_id"],
                name=site_data["name"],
                coordinates=site_data["coordinates"]
            )
            
            # Add description to properties
            binding_site.properties['description'] = site_data.get("description", "")
            
            # Store the binding site
            binding_sites[site_id] = binding_site
            logger.info(f"Created sample binding site: {binding_site.name}")
            
        except Exception as e:
            logger.error(f"Error creating sample binding site {site_data['name']}: {str(e)}")
    
    return binding_sites


def create_sample_docking_results():
    docking_results_dict = {}
    
    # Create a sample docking result
    if molecules and binding_sites:
        first_molecule = next(iter(molecules.values()))
        first_site = next(iter(binding_sites.values()))
        
        result = DockingResult(first_molecule, first_site)
        result.id = "test_dock1"
        result.score = -8.5
        result.binding_affinity = -9.2
        result.pose_data = [1.0, 2.0, 3.0]
        
        docking_results_dict[result.id] = result
    
    return docking_results_dict