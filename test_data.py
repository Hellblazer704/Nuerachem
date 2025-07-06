from models import Molecule, BindingSite, DockingResult

def create_sample_data():
    # Create a sample molecule (Aspirin)
    molecule = Molecule(
        mol_id="test_mol1",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        name="Aspirin"
    )

    # Create a sample binding site
    binding_site = BindingSite(
        site_id="test_site1",
        protein_id="1ABC",
        name="Active Site",
        coordinates=[10.0, 15.0, 20.0]
    )

    # Create a sample docking result
    docking_result = DockingResult(molecule, binding_site)
    docking_result.id = "test_dock1"
    docking_result.score = -8.5
    docking_result.pose_data = [12.0, 17.0, 22.0]
    
    return molecule, binding_site, docking_result