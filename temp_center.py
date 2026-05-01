
from schrodinger import structure
import numpy as np
import sys
try:
    st = next(structure.StructureReader("1XKK_prepared.maegz"))
    parts = "A:797".split(':')
    chain = parts[0].strip()
    resnum = parts[1].strip()
    atoms = [a for a in st.atom if a.chain.strip() == chain.strip() and a.resnum == int(resnum)]
    if not atoms: sys.exit(1)
    center = np.mean([a.xyz for a in atoms], axis=0)
    print(f"{center[0]},{center[1]},{center[2]}")
except Exception:
    sys.exit(1)
