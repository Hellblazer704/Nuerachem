class MoleculeViewer {
    constructor(containerId) {
        this.viewer = $3Dmol.createViewer(document.getElementById(containerId), {
            backgroundColor: 'white'
        });
        this.currentAnimation = null;
    }

    async loadMolecule(moleculeId) {
        try {
            const response = await fetch(`/api/molecules/${moleculeId}/3d`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            const data = await response.json();
            if (!data.pdb_data) {
                throw new Error('No structure data received');
            }
            this.displayStructure(data.pdb_data);
        } catch (error) {
            console.error('Error loading molecule:', error);
            throw error;
        }
    }

    displayStructure(pdbData) {
        try {
            this.viewer.clear();
            this.viewer.addModel(pdbData, "pdb");
            this.viewer.setStyle({}, {
                stick: {},
                cartoon: { color: 'spectrum' }
            });
            this.viewer.zoomTo();
            this.viewer.render();
        } catch (error) {
            console.error('Error displaying structure:', error);
            throw error;
        }
    }
}