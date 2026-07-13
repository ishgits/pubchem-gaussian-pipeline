"""
pubchem-gaussian-pipeline
=========================

Automated pipeline for generating Gaussian input files and SLURM
submission scripts from molecule names, PubChem CIDs, or SMILES strings.

Pipeline steps
--------------
1. Resolve molecule names → PubChem CIDs + properties  (pubchem)
2. Conformer search: RDKit ETKDGv3 + MMFF94 rank,
   top-3 distinct conformers per molecule               (conformers)
3. Write Gaussian .com input files (opt + freq),
   one per conformer                                     (gaussian)
4. Generate SLURM .sh submission scripts                (slurm)

The v1.1 path (PubChem 3D SDF → Open Babel XYZ, via ``download_sdfs`` /
``convert_sdfs_to_xyz``) remains available for single-geometry use.
"""

# Pipeline version for provenance (M-06). Bump manually when the code that
# produces scientific outputs changes. The v2 conformer pathway is the released
# default at 2.0.0 (MOD-01). Recorded in conformer_log.csv alongside a
# best-effort git commit so two runs from different revisions are
# distinguishable; changing this string deliberately changes the
# pipeline_version stamped into provenance logs going forward.
__version__ = "2.0.0"

from .pubchem import (
    build_molecule_table,
    download_sdfs,
    resolve_pubchem_record,
    resolve_with_fallback,
    score_candidate,
)
from .geometry import sdf_to_xyz, convert_sdfs_to_xyz
from .conformers import (
    UNCONVERGED_FF_SEED,
    check_conformer_eligibility,
    generate_conformers,
    search_conformers,
    select_converged_top_n,
    select_top_n,
)
from .gaussian import (
    xyz_to_gaussian_coords,
    write_gaussian_com,
    write_gaussian_coms,
    write_gaussian_coms_from_conformers,
)
from .slurm import write_slurm_script, write_slurm_scripts
from .utils import (
    sanitize_basename,
    ensure_dir,
    normalize_cid,
    git_short_sha,
    pipeline_provenance,
)
