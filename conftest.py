import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parent))
