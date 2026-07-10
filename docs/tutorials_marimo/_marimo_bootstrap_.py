import sys
from pathlib import Path


def add_repo_root(levels: int = 2):
    root = Path(__file__).resolve().parents[levels]
    sys.path.insert(0, str(root))
