import sys
from pathlib import Path

# Ensure the package is importable whether or not it is installed editable.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
