"""
pubchem-gaussian-pipeline
=========================

Automated pipeline for generating Gaussian input files and SLURM
submission scripts from molecule names, PubChem CIDs, or SMILES strings.

Pipeline steps
--------------
1. Resolve molecule names → PubChem CIDs + properties  (pubchem)
2. Download 3D SDF files from PubChem                  (pubchem)
3. Convert SDF → XYZ via Open Babel                    (geometry)
4. Write Gaussian .com input files (opt + freq)        (gaussian)
5. Generate SLURM .sh submission scripts               (slurm)
"""

from .pubchem import (
    build_molecule_table,
    download_sdfs,
    resolve_pubchem_record,
    resolve_with_fallback,
    score_candidate,
)
from .geometry import sdf_to_xyz, convert_sdfs_to_xyz
from .gaussian import xyz_to_gaussian_coords, write_gaussian_com, write_gaussian_coms
from .slurm import write_slurm_script, write_slurm_scripts
from .utils import sanitize_basename, ensure_dir, normalize_cid
