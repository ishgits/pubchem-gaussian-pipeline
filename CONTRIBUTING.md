# Contributing

## Scientific invariants

These are load-bearing. Do not change one silently — if a change needs to
alter any of them, say so explicitly in the PR description.

- Route lines are chemistry. Do not alter the functional, basis set, solvent,
  temperature, or optimization/frequency semantics without an explicit,
  recorded scientific decision.
- Charge and multiplicity are physical. Defaults are neutral singlet, and any
  change in handling must be explicit.
- Coordinates are Ångström, conformer force-field energies are kcal/mol, and
  DFT energies are Hartree. Never mix or silently convert them.
- Never place illustrative, mocked, random, or placeholder science in a
  non-test pipeline path.
- The Link1 frequency step must read geometry and wavefunction from the
  optimization checkpoint through `Geom=AllChk Guess=Read`.
- RDKit/MMFF or UFF geometries are DFT starting points, not validated minima.
- PubChem's stereo-bearing SMILES is the source of record. Molecules with
  undefined stereochemistry take the provisional path (one loudly flagged
  arbitrated structure, `dE=NA`) — they are never silently guessed and never
  presented as the defined stereoisomer.
- UFF is a recorded fallback only. Relative force-field energies are
  meaningful only within one molecule and one force-field method.
- The conformer stage does not infer electronic multiplicity. Gaussian charge
  and multiplicity remain explicit inputs.

The provenance guarantees in [`docs/provenance.md`](docs/provenance.md) are
part of the contract: the manifest is authoritative, artifacts carry one
`artifact_id` back-pointer, runs are immutable with no resume/append, and
collisions or hash disagreements fail before any output is written.

## Development rules

- Work on a branch; do not commit directly to `main`.
- Keep `tests/` offline and independent of Gaussian, a cluster, and
  Open Babel.
- Every generated artifact must be traceable through the manifest lineage:
  molecule identity → conformer → XYZ → COM → SH.
- Validate before mutation. A failed stage must not partially rewrite
  existing logs or outputs.
- One input record maps to exactly one unique output record and path.
- Do not claim validation because execution succeeded. Validation requires an
  expected value, invariant, reference, or explicit mapping check.

## Required checks

Before opening a PR:

```bash
pytest tests/ -q
python scripts/check_invariants.py
git diff --check
test -z "$(git ls-files -ci --exclude-standard)"
```

For a release candidate, repeat the same checks from a clean `git archive`
and note whether the pinned dependency stack in `requirements-lock.txt` was
actually used.
