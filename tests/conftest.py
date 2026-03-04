import os
import sys

# ensure the maws package in workspace root is importable
root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if root not in sys.path:
    sys.path.insert(0, root)
