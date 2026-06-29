import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# Make utils_sav_io importable when pytest runs from the repo root
sys.path.insert(0, str(Path(__file__).parent))
