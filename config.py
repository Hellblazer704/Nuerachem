"""
Configuration settings for the Neural Chem AI application
"""

import os

# Flask Configuration
DEBUG = os.environ.get('FLASK_DEBUG', 'True').lower() in ('true', '1', 'yes')
TESTING = os.environ.get('FLASK_TESTING', 'False').lower() in ('true', '1', 'yes')
SECRET_KEY = os.environ.get('SESSION_SECRET', 'neural-chem-ai-default-key')

# RDKit Configuration
MAX_CONFORMERS = int(os.environ.get('MAX_CONFORMERS', 50))  # Maximum number of conformers to generate per molecule
CONFORMER_ENERGY_THRESHOLD = float(os.environ.get('CONFORMER_ENERGY_THRESHOLD', 10.0))  # Energy threshold for conformer filtering (kcal/mol)

# Machine Learning Configuration
MODEL_WEIGHTS_PATH = os.environ.get('MODEL_WEIGHTS_PATH', 'ml_models/weights')

# API Configuration
MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB max upload size

# Visualization Configuration
VISUALIZATION_MODE = os.environ.get('VISUALIZATION_MODE', 'interactive')  # Options: 'static', 'interactive'
