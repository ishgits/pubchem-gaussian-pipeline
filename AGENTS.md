# AGENTS.md — pubchem-gaussian-pipeline

Universal operating instructions for any coding or reviewing agent working in
this repository. This file is the source of truth. The exact frozen v2.0 release
boundary is defined in `docs/release-contract-v2.0.md` and is normative by
reference.

## 1. Project objective

Generate Gaussian input files and SLURM submission scripts from molecule names.
The primary v2 pathway is:

```text
PubChem stereo SMILES -> RDKit conformer search -> ranked XYZ artifacts
-> Gaussian Link1 opt/freq inputs -> SLURM scripts
```

The pipeline produces inputs; it does not run Gaussian. The Open Babel v1.1
single-geometry pathway is deprecated compatibility functionality and is not the
default.

## 2. Scientific invariants

These are load-bearing. A silent change is prohibited.

- Route lines are chemistry. Do not alter the functional, basis set, solvent,
  temperature, or optimization/frequency semantics without an explicit approved
  deviation.
- Charge and multiplicity are physical. Defaults are neutral singlet, and any
  change in handling must be explicit.
- Coordinates are Ångström, conformer force-field energies are kcal/mol, and DFT
  energies are Hartree. Never mix or silently convert them.
- Never place illustrative, mocked, random, or placeholder science in a
  non-test pipeline path.
- The Link1 frequency step must read geometry and wavefunction from the
  optimization checkpoint through `Geom=AllChk Guess=Read` unless an approved
  deviation says otherwise.
- RDKit/MMFF or UFF geometries are DFT starting points, not validated minima.
- PubChem stereo-bearing SMILES is the source of record. Undefined stereo is
  skipped and logged rather than guessed.
- UFF is a recorded fallback only. Relative force-field energies are meaningful
  only within one molecule and one force-field method.
- The conformer stage does not infer electronic multiplicity. Gaussian charge
  and multiplicity remain explicit inputs.

## 3. Frozen v2.0 provenance and artifact contract

Follow `docs/release-contract-v2.0.md` exactly.

- `run_manifest.json` is the authoritative complete provenance record.
- XYZ and COM are scientific artifacts. They carry `run_id`, `artifact_id`,
  `config_hash`, and the exact per-artifact fields in the release contract.
- SLURM scripts are operational artifacts and identify their exact source COM.
- The supported transfer/archive unit is the complete run package, including the
  manifest, stage logs, XYZ, COM, and SH directories.
- Full conformer-search configuration belongs in the manifest. Do not require or
  add duplicate copies of every knob to each XYZ or COM unless the frozen
  contract is deliberately revised by Ish.
- Flat directories remain for v2.0. The `runs/<study>/<run_id>/` redesign is
  deferred to v2.1.
- All source and destination path collisions, duplicate artifact IDs, missing or
  zero-byte required files, and manifest/hash disagreements fail before output
  mutation.
- Resume and append reuse require the same clean, nonblank source commit. Dirty
  or unavailable git provenance disables reuse and forces regeneration.
- The legacy Open Babel v1.1 path is explicitly exempt from the strict v2
  manifest contract.

## 4. Development rules

- Work on a branch; never commit directly to `main`.
- Keep `tests/` offline and independent of Gaussian, a cluster, and Open Babel.
- Every generated v2 artifact must be traceable through the manifest lineage:
  molecule identity -> conformer -> XYZ -> COM -> SH.
- Validate before mutation. A failed stage must not partially rewrite existing
  logs or outputs.
- One input record must map to exactly one unique output record and path.
- Do not claim validation because execution succeeded. Validation requires an
  expected value, invariant, reference, or explicit mapping check.
- Maintain `docs/implementation-status.md` honestly. Planned behavior is not
  implemented behavior.

## 5. Required checks

The objective floor must be green before review:

```bash
pytest tests/ -q
python scripts/check_invariants.py
git diff --check
test -z "$(git ls-files -ci --exclude-standard)"
```

A final release candidate must also pass the same checks from a clean
`git archive` and must record whether the pinned dependency stack was actually
used.

## 6. Definition of done

A change is done only when:

1. required checks are green;
2. the implementation matches `docs/architecture.md`,
   `docs/implementation-plan.md`, and `docs/release-contract-v2.0.md`;
3. `docs/implementation-status.md` accurately distinguishes implemented,
   pending, rejected, and deferred work;
4. all contract-level Blocker and Major findings are resolved or rejected with
   an approved contract-based justification;
5. one final holistic base-to-head re-review finds no contract-level Blocker or
   Major regression;
6. Ish gives final merge approval.

## 7. Review guidelines

Review the actual base-to-head diff and code. Never trust the status document
without verification.

### Holistic audit dimensions

Check all supported public entry points and every stage handoff for:

```text
scientific invariants
molecular identity
complete-group integrity
manifest and artifact linkage
file hashes
one-to-one path mapping
collision handling
failure-before-mutation
resume and append behavior
zero-job behavior
clean-archive reproducibility
```

Test normal, blank, missing, zero-byte, stale, duplicate, colliding, damaged,
and no-git inputs where relevant.

### Severity boundary

A release-blocking finding must demonstrate:

- wrong chemistry or a scientific-invariant violation;
- silent data loss or overwrite;
- identity corruption;
- violation of the frozen release contract;
- a supported workflow producing invalid jobs;
- required-check regression.

Treat silent route/unit/charge/Link1 changes and placeholder science as Blockers.
Treat broken manifest lineage, missing required contract metadata, unverifiable
reuse, or non-runnable supported jobs as Majors.

Do **not** report absence of every duplicated search knob inside an XYZ or COM as
a finding when the artifact satisfies the exact per-artifact matrix and links to
a complete valid manifest. Suggestions that broaden the frozen contract belong
in a deferred v2.1 section unless Ish explicitly approves the expansion.

### Finding format and re-review

- Classify findings as Blocker, Major, Moderate, Minor, or Verified Strength.
- Give each finding a stable non-reused ID.
- State location, evidence, consequence, required remediation, and verification.
- On re-review classify prior findings as Resolved, Partially resolved,
  Unresolved, Rejected with justification, or Regressed.
- State explicitly whether each new finding violates the frozen contract or is a
  proposed contract expansion.

### Stop rule

After the implementation is aligned to the frozen contract, conduct one
holistic audit and one final re-review. Recommendations that do not demonstrate a
contract violation are recorded for v2.1 and do not keep v2.0 in an open-ended
review loop.
