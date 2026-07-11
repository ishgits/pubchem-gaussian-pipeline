# AGENTS.md — pubchem-gaussian-pipeline

Universal operating instructions for any coding or reviewing agent (Claude Code,
Codex, Cowork) working in this repository. This file is the single source of
truth. `CLAUDE.md` and the `prompts/` briefs point back here.

## 1. Project objective

Generate Gaussian input files (`.com`) and SLURM submission scripts (`.sh`) from
molecule names, via PubChem lookup → 3D SDF → Open Babel XYZ → Gaussian opt+freq
(Link1 pattern). The pipeline produces *inputs* to quantum-chemistry jobs; it
does not run Gaussian. Correctness means the generated inputs faithfully encode
the intended chemistry and are fully reproducible from recorded configuration.

## 2. Scientific invariants (never violate silently)

These are load-bearing. Changing any of them alters the science and MUST be
recorded as an explicit deviation in `docs/implementation-status.md`, never made
quietly.

- **Route lines are chemistry, not strings.** The functional, basis set,
  solvation model (`scrf=(iefpcm,solvent=...)`), and `temperature=` in the
  opt/freq route lines determine the computed energies and thermochemistry.
  Never edit, reformat, or "clean up" a route line without recording why.
- **Charge and multiplicity are physical.** Defaults are neutral singlet
  (`0 1`). Do not change per-molecule charge/multiplicity handling without a
  recorded deviation.
- **Units are never changed silently.** Coordinates are Ångström, energies
  Hartree, temperature Kelvin. Any conversion must be explicit and labeled.
- **No placeholder science.** Never replace a real calculation, coordinate, or
  API result with a hardcoded, illustrative, random, or mocked value in a
  non-test code path. Tests may use fixtures; `pipeline/` may not.
- **The Link1 opt→freq contract holds.** The frequency job must read geometry
  from the optimization checkpoint (`Geom=AllChk Guess=Read`). Do not decouple
  opt and freq without a recorded deviation.
- **Starting geometries are approximate, and the code must say so.** PubChem +
  Open Babel geometries are DFT-*starting* points, not optimized minima. Do not
  add claims of geometric optimality.

## 3. Development rules

- Never silently change units.
- Never replace physical calculations or real API results with illustrative
  placeholders.
- Every generated scientific output must record the configuration that produced
  it (route lines, charge/mult, nproc) and the relevant software versions
  (pipeline version, Open Babel version).
- Generated outputs must be traceable to their inputs (molecule name → CID →
  SDF → XYZ → .com, via the log CSVs).
- Do not claim validation based only on successful execution. "It ran" ≠ "it is
  correct." Validation requires checking the output against an expected value or
  reference.
- Work on a branch, never commit directly to `main`.
- Keep `tests/` runnable with no network and no Open Babel / Gaussian
  (pure-function tests only).

## 4. Required checks (the objective floor)

CI (`.github/workflows/review-readiness.yml`) must be green before any human or
agent review begins:

- `pytest tests/ -q` passes.
- `python scripts/check_invariants.py` passes (placeholder-sentinel scan +
  route-line physics-token check).

If the floor is red, the change is not review-ready. Do not request review.

## 5. Definition of done

A change is done when, and only when:

1. Required checks are green.
2. `docs/implementation-status.md` is updated (what was/wasn't implemented,
   deviations, tests added, known limitations, questions requiring scientific
   judgment).
3. Any deviation from `docs/architecture.md` is recorded, not silent.
4. A reviewer's Blocker and Major findings are all Resolved or Rejected-with-
   justification.
5. Ish has given the final merge approval (human gate).

## 6. Review guidelines (read by Codex)

When reviewing a pull request in this repository:

- Inspect the actual diff and code. Do **not** trust
  `docs/implementation-status.md`; verify its claims against the code.
- Review for: correctness, scientific validity (§2 invariants),
  reproducibility, failure handling, testing, documentation, and scope
  compliance (does the PR do what the plan said, and nothing silently extra).
- Treat any silent change to a route line, unit, charge/multiplicity, or the
  Link1 contract as at least a **Blocker**.
- Treat placeholder/mocked values in `pipeline/` as a **Blocker**.
- Treat "validated because it ran" claims as at least a **Major**.
- Treat missing provenance (config/version not recorded in output) as a
  **Major**.
- Classify every finding as **Blocker / Major / Moderate / Minor / Verified
  Strength**. Every finding must state: (1) evidence — file and line range,
  (2) consequence, (3) required remediation.
- Give every finding a stable ID: `B-01`, `M-01`, `MOD-01`, `MIN-01`. IDs must
  not be reused or renumbered across rounds.
- On re-review, classify each prior finding as **Resolved / Partially resolved
  / Unresolved / Rejected (with justification) / Regressed**, and check the fix
  for new regressions.
