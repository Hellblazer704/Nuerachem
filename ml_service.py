"""
Machine Learning service - handles ML-based predictions and model management
"""

import logging
import os
import numpy as np
import random
from typing import Dict, List, Optional, Tuple, Any
from models import Molecule, BindingSite
import services.molecule_service as molecule_service

logger = logging.getLogger(__name__)

# Mock feature extraction
def extract_features_from_molecule(molecule: Molecule) -> np.ndarray:
    """
    Extract numeric features from a molecule for ML models
    
    Args:
        molecule: Molecule object
        
    Returns:
        Numpy array of features
    """
    try:
        # In a real implementation, this would use sophisticated chemistry-aware featurization
        # Examples include Morgan fingerprints, molecular descriptors, etc.
        
        # Get properties
        properties = molecule.properties or {}
        
        # Extract basic features (with defaults if missing)
        features = []
        
        # Molecular weight
        weight = properties.get('weight', 0)
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except:
                weight = 300
        features.append(weight)
        
        # LogP (octanol-water partition coefficient)
        logp = properties.get('logp', 0)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 0
        features.append(logp)
        
        # Number of atoms
        num_atoms = properties.get('num_atoms', 0)
        if isinstance(num_atoms, str):
            try:
                num_atoms = int(num_atoms)
            except:
                num_atoms = 20
        features.append(num_atoms)
        
        # Number of rings
        num_rings = properties.get('num_rings', 0)
        if isinstance(num_rings, str):
            try:
                num_rings = int(num_rings)
            except:
                num_rings = 0
        features.append(num_rings)
        
        # Normalize feature vector to unit length
        features = np.array(features, dtype=np.float32)
        
        # Pad with zeros to fixed length (for demonstration)
        padded = np.zeros(10, dtype=np.float32)
        padded[:len(features)] = features
        
        return padded
        
    except Exception as e:
        logger.error(f"Error extracting features: {str(e)}")
        # Return zero vector if extraction fails
        return np.zeros(10, dtype=np.float32)

def property_prediction(molecule: Molecule) -> Dict[str, Any]:
    """
    Predict additional molecular properties using ML
    
    Args:
        molecule: Molecule object
        
    Returns:
        Dictionary of predicted properties
    """
    try:
        # Extract features
        features = extract_features_from_molecule(molecule)
        
        # In a real implementation, we would load pretrained models
        # For this demo, we'll generate reasonable values
        
        # Get existing properties
        properties = molecule.properties or {}
        
        # Get molecular weight and logP for predictions
        weight = properties.get('weight', 300)
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except:
                weight = 300
                
        logp = properties.get('logp', 2)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 2
        
        # Calculate predictions
        predictions = {}
        
        # Solubility prediction (logS)
        # Rule of thumb: logS ~ 0.5 - logP
        predictions['water_solubility'] = max(-7, min(1, 0.5 - logp + random.uniform(-0.5, 0.5)))
        
        # Bioavailability prediction
        # Rule of thumb: MW < 500 and logP < 5 tends to have better bioavailability
        base_bioavailability = 1.0
        if weight > 500:
            base_bioavailability -= (weight - 500) / 500
        if logp > 5:
            base_bioavailability -= (logp - 5) / 5
        predictions['bioavailability'] = max(0, min(1, base_bioavailability + random.uniform(-0.1, 0.1)))
        
        # Blood-brain barrier penetration
        # Rule of thumb: small molecules with moderate logP tend to cross better
        predictions['bbb_penetration'] = (1 if (weight < 400 and 1 < logp < 4) else 0.2) + random.uniform(-0.1, 0.1)
        
        # Half-life prediction (hours)
        # Complex property influenced by many factors
        predictions['half_life'] = 4 + weight/100 + random.uniform(0, 5)
        
        # Format values for display
        for key, value in predictions.items():
            if isinstance(value, float):
                if -1 < value < 1:  # Small values
                    predictions[key] = round(value, 2)
                else:  # Larger values
                    predictions[key] = round(value, 1)
        
        return predictions
        
    except Exception as e:
        logger.error(f"Error in property prediction: {str(e)}")
        # Return default predictions if calculation fails
        return {
            'water_solubility': -2.0,
            'bioavailability': 0.5,
            'bbb_penetration': 0.3,
            'half_life': 5.0
        }

def binding_affinity_prediction(molecule: Molecule, binding_site: BindingSite) -> float:
    """
    Predict binding affinity between a molecule and binding site using ML
    
    Args:
        molecule: Molecule object
        binding_site: BindingSite object
        
    Returns:
        Predicted binding affinity in kcal/mol
    """
    try:
        # Extract molecule features
        mol_features = extract_features_from_molecule(molecule)
        
        # Get site properties for basic features
        site_properties = binding_site.properties or {}
        
        # In a real implementation, this would use a sophisticated ML model
        # For this demo, we'll simulate predictions
        
        # Base binding affinity - negative values indicate stronger binding
        # Typical range: 0 to -15 kcal/mol
        
        # Get key parameters
        weight = molecule.properties.get('weight', 300)
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except:
                weight = 300
                
        logp = molecule.properties.get('logp', 2)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 2
                
        # Size matching component
        # Optimal binding occurs at specific size ratios
        site_volume = site_properties.get('volume', 500)
        size_match = -5 * (1 - min(1, abs(weight/3 - site_volume/10) / 100))
        
        # Polarity matching component
        site_polarity = site_properties.get('polarity', 0.5)
        mol_polarity = max(0, min(1, 0.5 - logp/10))
        polarity_match = -4 * (1 - abs(site_polarity - mol_polarity))
        
        # Flexibility penalty
        num_rings = molecule.properties.get('num_rings', 0)
        if isinstance(num_rings, str):
            try:
                num_rings = int(num_rings)
            except:
                num_rings = 0
        flexibility_penalty = -1 if num_rings > 2 else 0
        
        # Combine components and add noise
        binding_affinity = size_match + polarity_match + flexibility_penalty + random.uniform(-1, 1)
        
        # Ensure reasonable range
        binding_affinity = max(-15, min(0, binding_affinity))
        
        return round(binding_affinity, 1)
        
    except Exception as e:
        logger.error(f"Error in binding affinity prediction: {str(e)}")
        # Return moderate affinity if prediction fails
        return -7.5

def drug_likeness_score(molecule: Molecule) -> Dict[str, Any]:
    """
    Calculate drug-likeness score and related properties
    
    Args:
        molecule: Molecule object
        
    Returns:
        Dictionary with drug-likeness scores and properties
    """
    try:
        # Get properties
        properties = molecule.properties or {}
        
        # Extract key properties
        weight = properties.get('weight', 300)
        if isinstance(weight, str):
            try:
                weight = float(weight)
            except:
                weight = 300
                
        logp = properties.get('logp', 2)
        if isinstance(logp, str):
            try:
                logp = float(logp)
            except:
                logp = 2
        
        # Check Lipinski's Rule of Five
        # 1. MW < 500
        # 2. logP < 5
        # 3. H-bond donors < 5 (not calculated in our simplified model)
        # 4. H-bond acceptors < 10 (not calculated in our simplified model)
        lipinski_violations = 0
        if weight > 500:
            lipinski_violations += 1
        if logp > 5:
            lipinski_violations += 1
            
        # QED (Quantitative Estimate of Drug-likeness)
        # Simplified calculation - real QED uses multiple descriptors
        # QED ranges from 0 to 1, where 1 is most drug-like
        
        # Desirability functions for key properties
        d_weight = 1.0 if 200 <= weight <= 400 else max(0, 1 - abs(weight - 300) / 300)
        d_logp = 1.0 if 1 <= logp <= 4 else max(0, 1 - abs(logp - 2.5) / 5)
        
        # Combine and add noise
        qed = (d_weight * d_logp) ** 0.5 + random.uniform(-0.1, 0.1)
        qed = max(0, min(1, qed))
        
        # Drug Safety Score (invented metric)
        # Higher score means safer
        safety_score = 1.0 - logp/10 - weight/1000 + random.uniform(-0.1, 0.1)
        safety_score = max(0, min(1, safety_score))
        
        # Format results
        result = {
            'lipinski_violations': lipinski_violations,
            'lipinski_pass': 'Yes' if lipinski_violations <= 1 else 'No',
            'qed_score': round(qed, 2),
            'drug_likeness': round(qed * 10, 1),  # Scale to 0-10
            'safety_score': round(safety_score * 10, 1)  # Scale to 0-10
        }
        
        return result
        
    except Exception as e:
        logger.error(f"Error in drug-likeness calculation: {str(e)}")
        # Return default values if calculation fails
        return {
            'lipinski_violations': 0,
            'lipinski_pass': 'Yes',
            'qed_score': 0.5,
            'drug_likeness': 5.0,
            'safety_score': 6.0
        }

def structure_similarity_search(query_molecule: Molecule, molecules: List[Molecule]) -> List[Tuple[Molecule, float]]:
    """
    Search for structurally similar molecules using ML
    
    Args:
        query_molecule: Molecule to use as query
        molecules: List of molecules to search through
        
    Returns:
        List of (molecule, similarity_score) pairs, sorted by similarity
    """
    try:
        # Extract features from query molecule
        query_features = extract_features_from_molecule(query_molecule)
        
        # Calculate similarity for each molecule
        similarities = []
        for mol in molecules:
            # Skip the query molecule itself
            if mol.id == query_molecule.id:
                continue
                
            # Extract features
            mol_features = extract_features_from_molecule(mol)
            
            # Calculate cosine similarity
            dot_product = np.dot(query_features, mol_features)
            norm_query = np.linalg.norm(query_features)
            norm_mol = np.linalg.norm(mol_features)
            
            if norm_query > 0 and norm_mol > 0:
                similarity = dot_product / (norm_query * norm_mol)
            else:
                similarity = 0
                
            # Add some noise for variety
            similarity = min(1.0, max(0, similarity + random.uniform(-0.1, 0.1)))
            
            similarities.append((mol, similarity))
            
        # Sort by similarity score (highest first)
        similarities.sort(key=lambda x: x[1], reverse=True)
        
        return similarities
        
    except Exception as e:
        logger.error(f"Error in similarity search: {str(e)}")
        # Return random ordering if calculation fails
        return [(mol, random.random()) for mol in molecules if mol.id != query_molecule.id]

def enrich_molecule_with_predictions(molecule: Molecule) -> Molecule:
    """
    Enrich a molecule with ML-based predictions for various properties
    
    Args:
        molecule: Molecule object
        
    Returns:
        Enriched Molecule object
    """
    try:
        # Skip if molecule has no properties
        if not hasattr(molecule, 'properties') or not molecule.properties:
            return molecule
            
        # Make predictions
        properties = property_prediction(molecule)
        drug_likeness = drug_likeness_score(molecule)
        
        # Add predictions to properties
        for key, value in properties.items():
            molecule.properties[f'pred_{key}'] = value
            
        for key, value in drug_likeness.items():
            molecule.properties[f'pred_{key}'] = value
            
        return molecule
        
    except Exception as e:
        logger.error(f"Error enriching molecule: {str(e)}")
        return molecule