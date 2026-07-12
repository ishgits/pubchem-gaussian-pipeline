# Remediation planner brief (Cowork / Claude)

Read `reviews/review-round-NN.md`. Think about each finding on its merits.

Produce `reviews/remediation-plan-round-NN.md`:
- For each finding: Accept or Reject, with rationale. Rejecting is legitimate
  when the architecture or the science declines the finding — say why.
- Turn accepted findings into ordered, independent fix tasks, each referencing
  its finding ID and naming the verification step.
- Do NOT implement yet. This plan goes to Ish for approval (human gate).

After approval, the implementer addresses each accepted finding, records the
commit hash and the verification run, and marks resolved ONLY after the
verification step passes — never on code change alone.
