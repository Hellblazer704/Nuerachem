from typing import Dict, List, Optional
from models import DockingResult
from storage import docking_results

class DockingService:
    def __init__(self):
        self.docking_results = docking_results

    def get_result(self, docking_id: str) -> Optional[DockingResult]:
        return self.docking_results.get(docking_id)

    def generate_animation_frames(self, docking_result: DockingResult) -> List[Dict]:
        if not docking_result or not docking_result.molecule:
            return []
            
        frames = []
        # Generate 30 frames for smooth animation
        for i in range(30):
            frame = {
                'frame_id': i,
                'pdb_data': self._generate_frame_data(docking_result, i/29),
                'progress': i/29,
                'score': docking_result.score * (i/29) if docking_result.score else 0
            }
            frames.append(frame)
        return frames

    def _generate_frame_data(self, docking_result: DockingResult, progress: float) -> str:
        # For testing, return a simple structure
        return """
ATOM      1  C   LIG     1       0.000   0.000   0.000
ATOM      2  C   LIG     1       1.500   0.000   0.000
END
"""