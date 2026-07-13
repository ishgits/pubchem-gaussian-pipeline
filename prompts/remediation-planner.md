# Remediation planner brief

Read the review, `AGENTS.md`, and `docs/release-contract-v2.0.md`.

For every finding determine first whether it demonstrates an actual frozen-
contract violation or proposes a broader design.

Produce a remediation plan that:

- accepts actual contract violations;
- may reject a requested implementation whose design conflicts with the frozen
  architecture, with an exact contract citation;
- still addresses the underlying risk through the approved architecture;
- records proposed contract expansions for v2.1 unless Ish explicitly accepts
  them for v2.0;
- turns accepted findings into independent ordered tasks with verification;
- requires validation before mutation and regression tests at adjacent public
  boundaries.

Do not implement until Ish approves the plan. Mark a finding resolved only after
its verification passes, not when code is merely changed.
