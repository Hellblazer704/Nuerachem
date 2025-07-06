"""
Binding site service - handling binding site operations and docking predictions
"""

import logging
import uuid
import math
import random
import numpy as np
from typing import Dict, List, Tuple, Optional
from models import BindingSite, Molecule, DockingResult

logger = logging.getLogger(__name__)

class BindingSiteError(Exception):
    """Exception raised for errors in binding site operations."""
    pass

def calculate_interaction_score(mol_properties: Dict, binding_site: BindingSite) -> float:
    """
    Calculate a docking score between a molecule and binding site
    
    Args:
        mol_properties: Dictionary of molecule properties
        binding_site: BindingSite object
        
    Returns:
        Interaction score (0-1 where 1 is best fit)
    """
    try:
        # Extract key properties for scoring
        # This is a simplified scoring system - real docking would be much more complex
        if not mol_properties:
            return random.uniform(0.1, 0.4)  # Lower score for molecules with no properties
            
        # Get molecular weight if available, otherwise estimate
        weight = mol_properties.get('weight', 0)
        if not weight and 'formula' in mol_properties:
            # Rough estimation based on formula if no weight is directly available
            formula = mol_properties.get('formula', '')
            weight = len(formula) * 10  # Very rough estimation
            
        # Binding pocket compatibility is partially based on size
        binding_pocket_size = 400  # Average binding pocket size in cubic angstroms
        size_compatibility = max(0, 1 - abs(weight - binding_pocket_size) / binding_pocket_size)
        
        # Charge and hydrophobicity also matter
        logp = mol_properties.get('logp', 0)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 0
                
        # Assume optimal logP is around 2 for binding
        logp_factor = max(0, 1 - abs(logp - 2) / 5)
        
        # Number of rings affects rigidity
        num_rings = int(mol_properties.get('num_rings', 0))
        if isinstance(num_rings, str):
            try:
                num_rings = int(num_rings)
            except:
                num_rings = 0
        
        # Slightly prefer molecules with 1-3 rings
        ring_factor = 0.5
        if 1 <= num_rings <= 3:
            ring_factor = 0.9
        
        # Calculate binding site specific factors from properties
        site_properties = binding_site.properties or {}
        binding_preference = site_properties.get('polarity', 0.5)  # Default moderate polarity
        
        # Calculate polarity match
        polarity_match = 1 - abs(binding_preference - (0.5 - logp/10))
        
        # Combine factors - weighted sum
        score = (
            0.3 * size_compatibility + 
            0.3 * logp_factor + 
            0.2 * ring_factor + 
            0.2 * polarity_match
        )
        
        # Add a small random factor to simulate the imprecision of docking predictions
        score = min(1.0, max(0.1, score + random.uniform(-0.1, 0.1)))
        
        return score
        
    except Exception as e:
        logger.error(f"Error calculating interaction score: {str(e)}")
        # Return a moderate score if calculation fails
        return 0.5

def predict_interactions(molecule: Molecule, binding_site: BindingSite) -> List[Dict]:
    """
    Predict the interactions between a molecule and binding site
    
    Args:
        molecule: Molecule object
        binding_site: BindingSite object
        
    Returns:
        List of interactions with type and strength
    """
    try:
        # Get molecule properties
        properties = molecule.properties or {}
        
        # If no properties available, return basic interactions
        if not properties:
            return [
                {'type': 'unknown', 'strength': 0.5}
            ]
            
        # Otherwise try to identify most likely interactions
        interactions = []
        
        # Check for hydrogen bond potential
        if any(atom in properties.get('formula', '') for atom in ['O', 'N']):
            interactions.append({
                'type': 'hydrogen_bond', 
                'strength': random.uniform(0.6, 0.9)
            })
            
        # Check for hydrophobic interaction potential based on logP
        logp = properties.get('logp', 0)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 0
                
        if logp > 2:
            interactions.append({
                'type': 'hydrophobic', 
                'strength': min(0.9, 0.5 + logp/10)
            })
            
        # Check for ionic interactions based on charged groups
        if any(group in properties.get('formula', '') for group in ['N', 'P', 'S']):
            interactions.append({
                'type': 'ionic', 
                'strength': random.uniform(0.5, 0.8)
            })
            
        # Add at least one interaction if none were found
        if not interactions:
            interactions.append({
                'type': 'van_der_waals', 
                'strength': random.uniform(0.3, 0.6)
            })
            
        return interactions
            
    except Exception as e:
        logger.error(f"Error predicting interactions: {str(e)}")
        # Return basic interaction if prediction fails
        return [{'type': 'interaction', 'strength': 0.5}]

def calculate_binding_pose(molecule: Molecule, binding_site: BindingSite) -> Optional[List[float]]:
    """
    Calculate a simulated binding pose for visualization
    
    Args:
        molecule: Molecule object
        binding_site: BindingSite object
        
    Returns:
        List of coordinates representing the pose or None if calculation fails
    """
    try:
        # This is a simplified simulation that returns a random pose
        # Real docking would calculate actual 3D coordinates
        
        # Use binding site coordinates as center
        center = binding_site.coordinates or [0, 0, 0]
        
        # Generate a random offset within a reasonable radius (5 Angstroms)
        radius = 5.0
        
        # Generate a point on a sphere with random radius up to 'radius'
        phi = random.uniform(0, 2 * math.pi)
        costheta = random.uniform(-1, 1)
        theta = math.acos(costheta)
        r = random.uniform(0, radius)
        
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        
        pose = [
            center[0] + x,
            center[1] + y,
            center[2] + z
        ]
        
        return pose
            
    except Exception as e:
        logger.error(f"Error calculating binding pose: {str(e)}")
        return None

def calculate_toxicity(molecule: Molecule) -> Dict:
    """
    Calculate toxicity predictions for a molecule
    
    Args:
        molecule: Molecule object
        
    Returns:
        Dictionary of toxicity predictions
    """
    try:
        # This is a simplified toxicity prediction model
        # Real toxicity prediction would use ML models and chemical databases
        
        properties = molecule.properties or {}
        
        # Get molecular weight and logP
        weight = properties.get('weight', 0)
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except:
                weight = 300  # Default molecular weight
                
        logp = properties.get('logp', 0)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 2  # Default logP
                
        # Lipinski's Rule of Five violations
        # MW < 500, logP < 5, H-bond donors < 5, H-bond acceptors < 10
        violations = 0
        if weight > 500:
            violations += 1
        if logp > 5:
            violations += 1
            
        # Simplified model for toxicity based on common rules
        # Higher logP and molecular weight often correlate with increased toxicity
        
        # Base prediction - scale from 0 to 1 where higher is more toxic
        base_toxicity = (logp / 10) * 0.5 + (weight / 1000) * 0.5
        
        # Adjust and add random variation to simulate uncertainty
        toxicity = min(1.0, max(0, base_toxicity + random.uniform(-0.1, 0.1)))
        
        # Define categories
        category = "Low"
        if toxicity > 0.7:
            category = "High"
        elif toxicity > 0.4:
            category = "Moderate"
            
        # Return toxicity predictions
        return {
            'toxicity_score': round(toxicity, 2),
            'toxicity_category': category,
            'lipinski_violations': violations,
            'confidence': 'Low' if violations > 1 else 'Moderate'
        }
            
    except Exception as e:
        logger.error(f"Error calculating toxicity: {str(e)}")
        # Return default values if calculation fails
        return {
            'toxicity_score': 0.5,
            'toxicity_category': 'Moderate',
            'lipinski_violations': 0,
            'confidence': 'Low'
        }

def process_binding_site(binding_site: BindingSite) -> BindingSite:
    """
    Process a binding site to add additional properties and metadata
    
    Args:
        binding_site: BindingSite object
        
    Returns:
        Updated BindingSite object
    """
    try:
        # Initialize properties dictionary if it doesn't exist
        if not hasattr(binding_site, 'properties') or not binding_site.properties:
            binding_site.properties = {}
            
        # Add default coordinates if missing
        if not binding_site.coordinates or len(binding_site.coordinates) < 3:
            binding_site.coordinates = [0.0, 0.0, 0.0]
            
        # Add some example properties useful for docking
        if 'volume' not in binding_site.properties:
            binding_site.properties['volume'] = random.uniform(100, 1000)  # in cubic angstroms
            
        if 'polarity' not in binding_site.properties:
            binding_site.properties['polarity'] = random.uniform(0.2, 0.8)  # 0 to 1
            
        if 'flexibility' not in binding_site.properties:
            binding_site.properties['flexibility'] = random.uniform(0.1, 0.9)  # 0 to 1
            
        if 'charge' not in binding_site.properties:
            binding_site.properties['charge'] = random.uniform(-2, 2)  # net charge
            
        return binding_site
            
    except Exception as e:
        logger.error(f"Error processing binding site: {str(e)}")
        return binding_site

def predict_docking(molecule: Molecule, binding_site: BindingSite) -> DockingResult:
    """
    Predict docking between a molecule and binding site
    
    Args:
        molecule: Molecule object
        binding_site: BindingSite object
        
    Returns:
        DockingResult object with predictions
    """
    try:
        # Create new docking result
        result = DockingResult(molecule=molecule, binding_site=binding_site)
        
        # Calculate interaction score
        score = calculate_interaction_score(molecule.properties, binding_site)
        result.score = score
        
        # Predict binding pose
        result.pose_data = calculate_binding_pose(molecule, binding_site)
        
        # Predict interactions
        result.interactions = predict_interactions(molecule, binding_site)
        
        # Add toxicity prediction to results
        result.toxicity = calculate_toxicity(molecule)
        
        return result
            
    except Exception as e:
        logger.error(f"Error predicting docking: {str(e)}")
        # Return partial result if prediction fails
        result = DockingResult(molecule=molecule, binding_site=binding_site)
        result.score = 0.5
        return result