import logging
from app import app

# Configure logging with more detailed format
logging.basicConfig(
    level=logging.DEBUG, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    try:
        # Running the Flask application on host 0.0.0.0 and port 5000
        logger.info("Starting Neural Chem AI application")
        app.run(host="0.0.0.0", port=5000, debug=True)
    except Exception as e:
        logger.error(f"Failed to start application: {e}", exc_info=True)
