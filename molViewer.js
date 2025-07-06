/**
 * Molecule Viewer - Handles molecule visualization in the frontend
 */

document.addEventListener('DOMContentLoaded', function() {
    initializeMoleculeViewers();
});

/**
 * Initialize molecule viewers for all molecule containers on the page
 */
function initializeMoleculeViewers() {
    const moleculeContainers = document.querySelectorAll('.molecule-container');
    
    moleculeContainers.forEach(container => {
        // Get the SVG element inside the container
        const svgElement = container.querySelector('.molecule-svg');
        
        if (svgElement) {
            // Set up zoom and pan behavior if needed
            setupSvgInteractions(svgElement);
            
            // Ensure SVG is visible and properly sized
            svgElement.style.width = '100%';
            svgElement.style.height = 'auto';
            svgElement.style.display = 'block';
            
            // Add any error handling/fallback for missing SVG content
            if (!svgElement.innerHTML.trim()) {
                svgElement.innerHTML = `
                    <svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
                        <rect width="100%" height="100%" fill="#f8f9fa" />
                        <text x="50%" y="50%" text-anchor="middle" dominant-baseline="middle" fill="#6c757d">
                            Visualization unavailable
                        </text>
                    </svg>
                `;
            }
        }
    });
}

/**
 * Set up interactive behaviors for SVG molecule visualizations
 * @param {HTMLElement} svgElement - The SVG element to add interactions to
 */
function setupSvgInteractions(svgElement) {
    let isDragging = false;
    let startX, startY;
    let currentTranslateX = 0;
    let currentTranslateY = 0;
    let currentScale = 1;
    
    // Create a wrapper around the SVG content for transformations
    const svgContent = svgElement.querySelector('svg');
    
    if (!svgContent) return;
    
    // Setup mousewheel zoom
    svgElement.addEventListener('wheel', function(e) {
        e.preventDefault();
        
        const scaleAmount = e.deltaY > 0 ? 0.9 : 1.1;
        currentScale *= scaleAmount;
        
        // Limit scaling
        currentScale = Math.max(0.5, Math.min(currentScale, 5));
        
        applyTransform();
    });
    
    // Setup mouse drag for panning
    svgElement.addEventListener('mousedown', function(e) {
        isDragging = true;
        startX = e.clientX - currentTranslateX;
        startY = e.clientY - currentTranslateY;
        svgElement.style.cursor = 'grabbing';
    });
    
    document.addEventListener('mousemove', function(e) {
        if (!isDragging) return;
        
        currentTranslateX = e.clientX - startX;
        currentTranslateY = e.clientY - startY;
        applyTransform();
    });
    
    document.addEventListener('mouseup', function() {
        isDragging = false;
        svgElement.style.cursor = 'grab';
    });
    
    // Double-click to reset view
    svgElement.addEventListener('dblclick', function() {
        currentScale = 1;
        currentTranslateX = 0;
        currentTranslateY = 0;
        applyTransform();
    });
    
    // Apply all transformations
    function applyTransform() {
        svgContent.style.transform = `translate(${currentTranslateX}px, ${currentTranslateY}px) scale(${currentScale})`;
        svgContent.style.transformOrigin = 'center';
    }
    
    // Set initial style
    svgContent.style.transition = 'transform 0.1s ease-out';
    svgElement.style.cursor = 'grab';
}

/**
 * Load molecule data from the API
 * @param {string} moleculeId - ID of the molecule to load
 * @returns {Promise} Promise resolving to molecule data
 */
function loadMoleculeData(moleculeId) {
    return fetch(`/api/molecules/${moleculeId}`)
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load molecule data');
            }
            return response.json();
        })
        .then(data => {
            if (data.status === 'success') {
                return data.data;
            } else {
                throw new Error(data.message || 'Failed to load molecule data');
            }
        });
}

/**
 * Display molecule properties in a container
 * @param {Object} molecule - Molecule data object
 * @param {HTMLElement} container - Container element to display properties
 */
function displayMoleculeProperties(molecule, container) {
    if (!molecule || !container) return;
    
    const propertiesContainer = container.querySelector('.molecule-properties') || container;
    
    // Clear existing content
    propertiesContainer.innerHTML = '';
    
    // Create a table for properties
    const table = document.createElement('table');
    table.className = 'table table-sm';
    
    // Add molecule name and basic info
    let html = `
        <thead>
            <tr>
                <th colspan="2">${molecule.name || 'Unnamed Molecule'}</th>
            </tr>
        </thead>
        <tbody>
    `;
    
    // Add SMILES if available
    if (molecule.smiles) {
        html += `
            <tr>
                <td>SMILES</td>
                <td><code>${molecule.smiles}</code></td>
            </tr>
        `;
    }
    
    // Add properties if available
    if (molecule.properties) {
        for (const [key, value] of Object.entries(molecule.properties)) {
            const formattedKey = key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
            html += `
                <tr>
                    <td>${formattedKey}</td>
                    <td>${value}</td>
                </tr>
            `;
        }
    }
    
    html += '</tbody>';
    table.innerHTML = html;
    propertiesContainer.appendChild(table);
}
