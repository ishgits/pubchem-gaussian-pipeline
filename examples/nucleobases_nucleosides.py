"""
Example molecule configuration: Nucleobases and Nucleosides
===========================================================

This is a real-world example from a study comparing canonical and
non-canonical nucleobases and their corresponding nucleosides (base + ribose).

Copy this into the "Define Your Molecules" cell in run_pipeline.ipynb and
modify to suit your project.
"""

# ── Molecule list ──────────────────────────────────────────────────────────
# Each entry is a human-readable label. For multi-component labels like
# "Adenine + Ribose", you MUST provide an alias (below) that maps to a
# single PubChem-resolvable name.

molecules = [
    "Ribose",
    "Water",
    # Canonical nucleobases
    "Adenine",
    "Thymine",
    "Cytosine",
    "Guanine",
    "Uracil",
    # Non-canonical / alternative nucleobases
    "Xanthine",
    "Hypoxanthine",
    "Barbituric acid",
    "2,6-Diaminopurine",
    "Purine",
    "Pyrimidine",
    "Isocytosine",
    "Imidazole",
    # Nucleosides (base + ribose)
    "Adenine + Ribose",
    "Thymine + Ribose",
    "Cytosine + Ribose",
    "Guanine + Ribose",
    "Uracil + Ribose",
    "Xanthine + Ribose",
    "Hypoxanthine + Ribose",
    "Barbituric acid + Ribose",
    "2,6-Diaminopurine + Ribose",
    "Purine + Ribose",
    "Pyrimidine + Ribose",
    "Isocytosine + Ribose",
    "Imidazole + Ribose",
]

# ── Alias map ──────────────────────────────────────────────────────────────
# Maps your label → the PubChem query string that returns the right record.
# Especially important for nucleosides, which have standard names in PubChem
# that differ from the "Base + Ribose" convention.

alias = {
    "Ribose": "beta-D-ribofuranose",
    "Water": "water",

    # Standard nucleosides
    "Adenine + Ribose": "adenosine",
    "Guanine + Ribose": "guanosine",
    "Cytosine + Ribose": "cytidine",
    "Uracil + Ribose": "uridine",
    "Hypoxanthine + Ribose": "inosine",
    "Xanthine + Ribose": "xanthosine",
    "Isocytosine + Ribose": "isocytidine",

    # Less standard nucleosides
    "Thymine + Ribose": "5-methyluridine",      # aka ribothymidine
    "Purine + Ribose": "nebularine",             # purine riboside
    "Pyrimidine + Ribose": "pyrimidine riboside",
    "Imidazole + Ribose": "imidazole riboside",
    "2,6-Diaminopurine + Ribose": "2,6-diaminopurine riboside",
    "Barbituric acid + Ribose": "barbituric acid riboside",
}

# ── Fallback queries ──────────────────────────────────────────────────────
# If the primary alias fails, these alternative names are tried in order.
# Useful for obscure molecules where PubChem naming is inconsistent.

fallback_queries = {
    "Imidazole + Ribose": [
        "1-(beta-D-ribofuranosyl)imidazole",
        "1-β-D-ribofuranosylimidazole",
        "N-ribosylimidazole",
    ],
    "Pyrimidine + Ribose": [
        "1-(beta-D-ribofuranosyl)pyrimidine",
        "1-β-D-ribofuranosylpyrimidine",
        "pyrimidine N-riboside",
    ],
    "Barbituric acid + Ribose": [
        "1-(beta-D-ribofuranosyl)barbituric acid",
        "1-(beta-D-ribofuranosyl)pyrimidine-2,4,6-trione",
        "N-ribosylbarbituric acid",
        "ribosylbarbituric acid",
    ],
}

# ── Expected formulas (optional) ──────────────────────────────────────────
# If provided, the scoring heuristic will prefer PubChem records whose
# molecular formula matches. Leave empty if you don't need this check.

expected_formulas = {}
