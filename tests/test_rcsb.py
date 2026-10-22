"""Offline test for the mutation confirmation in rcsb.py, using entries copied from
the REAL RCSB GraphQL responses (2026-10-20): pdbx_mutation null, gene glued to
the variant in the title ("KRASG12D"), isoform name "KRAS4B", other-protein entity.
Run: uv run python tests/test_rcsb.py
"""
from pathlib import Path

from protein_mcp.sources import rcsb

KRAS = {"rcsb_polymer_entity": {"pdbx_description": "Isoform 2B of GTPase KRas", "pdbx_mutation": None},
        "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P01116"]}}


def entry(pdb_id, title, entities=(KRAS,)):
    return {"rcsb_id": pdb_id, "struct": {"title": title}, "polymer_entities": list(entities)}


def parsed(pdb_id, title, entities=(KRAS,)):
    return rcsb._parse_pdb_entry(entry(pdb_id, title, entities), "P01116")


def main():
    # 9TBR: the known false negative (gene glued to the variant, empty mutation field)
    s = parsed("9TBR", "From KRASG12D to pan-KRAS inhibitors")
    assert s["mutations"] == [] and s["mutations_from_title"] == ["G12D"], s
    assert rcsb._mutation_matches(s, "G12D")
    assert rcsb._mutation_matches(s, "G12")                          # bare position
    assert not rcsb._mutation_matches(s, "G12C")

    # other real glued forms
    assert parsed("7OK3", "Crystal Structure of KRasG13C in Complex with Nucleotide-based "
                          "Covalent Inhibitor edaGDP")["mutations_from_title"] == ["G13C"]
    assert parsed("5UFQ", "K-RasG12D(GNP)/R11.1.6 complex")["mutations_from_title"] == ["G12D"]

    # 9E3S: G12N must still be rejected for G12D (the previous fix must hold)
    cypa = {"rcsb_polymer_entity": {"pdbx_description": "Peptidyl-prolyl cis-trans isomerase A",
                                    "pdbx_mutation": None},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P62937"]}}
    s = parsed("9E3S", "Tricomplex of RMC-9945, KRAS G12N, and CypA", (cypa, KRAS))
    assert s["mutations_from_title"] == ["G12N"] and s["partner_chains"]
    assert not rcsb._mutation_matches(s, "G12D")

    # isoform names are not mutations: "KRAS4B" must not yield "S4B"
    s = parsed("7T1F", "Crystal structure of GDP-bound T50I mutant of human KRAS4B")
    assert s["mutations_from_title"] == ["T50I"], s
    s = parsed("7RSC", 'NMR-driven structure of the KRAS4B-G12D "alpha-alpha" dimer on a lipid bilayer nanodisc')
    assert s["mutations_from_title"] == ["G12D"], s

    # a digit before the variant still blocks it (not a gene prefix)
    assert parsed("XXXX", "Complex 1G12D")["mutations_from_title"] == []
    print("ALL OK")


main()
