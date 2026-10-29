"""Offline test for analyze_interfaces (protein_mcp/files/interfaces.py, sources/pisa.py,
rcsb.fetch_assemblies) and its report section. Fixtures are REAL replies (2026-10-28):
3WD5 (TNF + adalimumab Fab; I 21 3, TNF trimer by crystal symmetry, RCSB assembly = asymmetric unit)
and 2W2N (PCSK9 prodomain + catalytic domain + LDLR EGF-A with His tags; lattice-shifted contacts):
mmCIFs, RCSB chains and assemblies GraphQL replies, PDBe PISA interfaces/multimers XML.
Run: uv run python tests/test_interfaces.py
"""
import asyncio
import gzip
import hashlib
import json
import tempfile
from pathlib import Path

import httpx

from protein_mcp import dossier, report
from protein_mcp.files import choice, interfaces
from protein_mcp.sources import pisa, rcsb

FIX = Path(__file__).parent / "fixtures" / "prep"
CHAINS = json.loads((FIX / "rcsb_chains_3WD5_2W2N.json").read_text(encoding="utf-8"))
ASSEMBLIES = (FIX / "rcsb_assemblies_3WD5_2W2N.json").read_bytes()
PISA = {(k, p): gzip.decompress((FIX / f"pisa_{k}_{p}.xml.gz").read_bytes())
        for k in ("interfaces", "multimers") for p in ("3wd5", "2w2n")}
ACC = {"3WD5": "P01375", "2W2N": "Q8NBP7"}
STATE = {"pisa_down": False, "multimers_down": False, "rcsb_down": False, "rewrite": None}


def handler(request: httpx.Request) -> httpx.Response:
    url = str(request.url)
    if url == rcsb.RCSB_GRAPHQL:
        return httpx.Response(500, text="down") if STATE["rcsb_down"] else httpx.Response(200, content=ASSEMBLIES)
    if "pisa/cgi-bin" in url:
        if STATE["rewrite"]:                         # prepare_structures runs meanwhile
            STATE["rewrite"]()
        if STATE["pisa_down"]:
            return httpx.Response(503, text="Service Unavailable")
        kind = "interfaces" if "interfaces.pisa" in url else "multimers"
        if STATE["multimers_down"] and kind == "multimers":       # real answer for cryo-EM 5A1A
            return httpx.Response(200, content=b"<pisa_multimers><status>Ok</status><pdb_entry><pdb_code>x</pdb_code>"
                                               b"<status>No symmetry operations</status></pdb_entry></pisa_multimers>")
        return httpx.Response(200, content=PISA[(kind, url.rsplit("?", 1)[1])])
    return httpx.Response(500, text=f"unexpected {url}")


def run(*args):
    original = httpx.AsyncClient
    httpx.AsyncClient = lambda **kw: original(transport=httpx.MockTransport(handler), **kw)
    try:
        return asyncio.run(interfaces.analyze_interfaces(*args))
    finally:
        httpx.AsyncClient = original


def make_manifest(pdb: str) -> Path:
    """A manifest record as prepare_structures writes it, with the real mmCIF and its SHA-256."""
    acc = ACC[pdb]
    prep = dossier.DOSSIER_DIR / acc / "prep"
    (prep / "raw").mkdir(parents=True, exist_ok=True)
    data = gzip.decompress((FIX / f"{pdb}.cif.gz").read_bytes())
    (prep / "raw" / f"{pdb}.cif").write_bytes(data)
    entry = next(e for e in CHAINS["data"]["entries"] if e["rcsb_id"] == pdb)
    rec = {"pdb_id": pdb, "file": {"ok": True, "path": f"raw/{pdb}.cif", "sha256": hashlib.sha256(data).hexdigest()},
           **rcsb._parse_chains(entry, acc), "links": {"rcsb": f"https://www.rcsb.org/structure/{pdb}"}}
    m = {"manifest_version": 1, "uniprot_id": acc, "created_at": "2026-10-28T10:00:00Z", "generated_by": "prepare_structures",
         "source_dossier": {"file": "structures.json"}, "structures": [rec], "warnings": [], "decisions": [], "pending": [],
         "user_choice": {"structure": pdb, "reason": "kept through the write"}}
    choice.write_manifest(prep / "manifest.json", m)
    return prep / "manifest.json"


def test_parsers():
    p = pisa.parse_interfaces(PISA[("interfaces", "3wd5")])
    assert p["ok"] and {i["id"] for i in p["interfaces"]} == set(range(1, 13))       # 12 protein-protein interfaces
    two = next(i for i in p["interfaces"] if i["id"] == 2)
    assert two["css"] == 1.0 and [(m["chain"], m["symop"]) for m in two["molecules"]] == [("A", "x,y,z"), ("A", "Z,X,Y")]
    # every set keeps ALL its assemblies (reviewer 2026-10-28: a deduplication emptied sets 2 and 3)
    mm = pisa.parse_multimers(PISA[("multimers", "3wd5")])
    assert [[a["composition"] for a in s["assemblies"]] for s in mm["sets"]] == [["A[3]", "LH"], ["A[3]", "L", "H"], ["LH", "A"]]
    assert mm["total_assemblies"] == 3 and mm["sets"][0]["assemblies"][0]["verdict"] == "This assembly appears to be stable in solution."
    assert not pisa.parse_interfaces(b"<pisa_interfaces><status>Ok</status><pdb_entry><status>Entry not found</status>"
                                     b"</pdb_entry></pisa_interfaces>")["ok"]
    # an interface without <id> is skipped and counted, not a crash (reviewer 2026-10-28)
    odd = pisa.parse_interfaces(b"<pisa_interfaces><status>Ok</status><pdb_entry><status>Ok</status><interface>"
                                b"<molecule><chain_id>A</chain_id><class>Protein</class><symop>x,y,z</symop></molecule>"
                                b"<molecule><chain_id>B</chain_id><class>Protein</class><symop>x,y,z</symop></molecule>"
                                b"</interface></pdb_entry></pisa_interfaces>")
    assert odd == {"ok": True, "interfaces": [], "unreadable_interfaces": 1}
    a = rcsb._parse_assemblies(next(e for e in json.loads(ASSEMBLIES)["data"]["entries"] if e["rcsb_id"] == "3WD5"))
    assert a["assemblies"] == [{"id": "3WD5-1", "defined_by": "author_defined_assembly", "method": None,
                                "oligomeric_count": 3, "oligomeric_details": "trimeric", "polymer_chains": ["A", "H", "L"],
                                "operators": ["1"], "equals_asymmetric_unit": True}]
    assert interfaces.sequence_flags("MKHHHHHHPMSXD")[1:] == [
        "1 non-standard or unknown residue (X)",
        "6 consecutive histidines at residues 3-8 (often an affinity tag; not stated by RCSB)"]


def test_3wd5():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    assert "run prepare_structures first" in run("P01375")["error"]
    path = make_manifest("3WD5")
    r = run("P01375", ["3wd5", "1ABC"])
    assert r["ok"], r
    assert r["selection"] == "1 of 1 structures in the manifest (requested: 3WD5, 1ABC)" and r["pdb_ids_not_in_manifest"] == ["1ABC"]
    m = json.loads(path.read_text(encoding="utf-8"))
    assert m["user_choice"]["reason"] == "kept through the write"
    b = m["structures"][0]["interfaces"]
    assert b["raw_sha256"] == m["structures"][0]["file"]["sha256"] and b["spacegroup"] == "I 21 3"
    by = {(c["partner_chain"], c["symmetry"]["code"]): c for c in b["contacts"]}
    # the TNF trimer is made by crystal symmetry: two 3-fold copies, both PISA interface 2 (css 1)
    aa = by[("A", "5_555")]
    assert aa["symmetry"]["operator"] == "z,x,y" and aa["n_target_residues"] == 31 and not aa["in_asymmetric_unit"]
    assert aa["partner_role"] == "target" and aa["pisa_interface"]["id"] == 2 and aa["pisa_interface"]["css"] == 1.0
    assert by[("A", "9_555")]["n_target_residues"] == 29
    # same chain, inverse operator: the agreement count uses the target-side molecule (reviewer: was "5 of 29")
    assert by[("A", "9_555")]["pisa_interface"]["target_residues_in_pisa"] == "28 of 29"
    assert aa["pisa_interface"]["target_residues_in_pisa"] == "31 of 31"
    # the Fab heavy chain touches its own TNF and the 3-fold neighbour (PISA 7 and 4, both css 0)
    assert (by[("H", "1_555")]["n_target_residues"], by[("H", "1_555")]["pisa_interface"]["id"]) == (7, 7)
    assert (by[("H", "9_555")]["n_target_residues"], by[("H", "9_555")]["pisa_interface"]["id"]) == (14, 4)
    lc = by[("L", "1_555")]
    assert lc["in_asymmetric_unit"] and lc["n_target_residues"] == 14 and lc["pisa_interface"]["css"] == 0.0
    assert lc["pisa_interface"]["target_residues_in_pisa"] == "14 of 14"
    # UniProt numbering through RCSB's alignment (TNF chain = soluble form; auth 6 -> UniProt 82)
    assert aa["target_residues"][0] == {"auth_seq_id": "6", "resname": "ARG", "uniprot": 82}
    assert all(set(c) >= {"target_residues", "symmetry"} and "_partner_ids" not in c for c in b["contacts"])
    assert b["rcsb_assemblies"][0]["equals_asymmetric_unit"] is True
    assert [a["composition"] for a in b["pisa"]["multimer_sets"][0]["assemblies"]] == ["A[3]", "LH"]
    seqs = {s["entity_id"]: s for s in b["sequences"]}
    assert seqs["1"]["role"] == "target" and seqs["1"]["length"] == 157
    assert seqs["3"]["description"] == "Adalimumab Heavy Chain" and seqs["3"]["uniprot_ids"] == []
    assert seqs["3"]["sequence"].startswith("EVQLVESGG") and seqs["3"]["chains"] == [{"label": "C", "auth": "H"}]
    # chat view: one line per contact, the sources quoted
    c = r["structures"][0]
    assert any(l.startswith("A-H (partner) 9_555 [y,z,x] symmetry copy: 14 target residues") for l in c["contacts"])
    assert c["rcsb_assemblies"] == ["3WD5-1: author_defined_assembly, polymer chains A, H, L, operators 1; same as the asymmetric unit"]
    assert c["pisa_multimer_sets"][0] == ("set 1: A[3] (This assembly appears to be stable in solution.) + "
                                          "LH (This assembly appears to be stable in solution.)")
    assert len(c["pisa_multimer_sets"]) == 3 and c["pisa_errors"] is None

    # report: section 6 subsection from the manifest only
    m["created_at"] = "2026-10-28T10:00:00Z"
    doc = report.Doc()
    report._interfaces(doc, m)
    md = report.render_markdown(doc)
    assert "1 of 1 structures analysed" in md and "PISA multimer set 1 (one alternative; its assemblies coexist): A[3]: " in md
    assert "| A | H (partner) | 9_555<br>y,z,x<br>symmetry copy | 14: 148, 149" in md and "#4, css 0, 445 Å², 14 of 14" in md
    assert "RCSB assembly 3WD5-1: author_defined_assembly" in md and "same as the asymmetric unit" in md
    empty = report.Doc()
    report._interfaces(empty, {"structures": [{"pdb_id": "X"}]})
    assert "Not analysed" in report.render_markdown(empty)

    # PISA down: said, contacts still computed; RCSB down: assemblies/sequences unknown (None, not [])
    STATE["pisa_down"] = True
    r = run("P01375")
    STATE["pisa_down"] = False
    assert r["ok"] and any("PISA interfaces not available" in w for w in r["warnings"]) and len(r["structures"][0]["pisa_errors"]) == 2
    # not "no matching interface": nothing is known from PISA (reviewer 2026-10-28)
    assert all(l.endswith("; PISA not available") for l in r["structures"][0]["contacts"])
    doc = report.Doc()
    report._interfaces(doc, json.loads(path.read_text(encoding="utf-8")))
    assert "| PISA not available |" in report.render_markdown(doc)
    # one PISA service fails alone (cryo-EM 5A1A): the other one is kept
    STATE["multimers_down"] = True
    r = run("P01375")
    STATE["multimers_down"] = False
    b = json.loads(path.read_text(encoding="utf-8"))["structures"][0]["interfaces"]
    assert b["pisa"]["interfaces"] and b["pisa"]["multimer_sets"] is None
    assert b["pisa"]["multimers_error"] == "PISA status: No symmetry operations"
    assert all(c["pisa_interface"] for c in b["contacts"])
    STATE["rcsb_down"] = True
    r = run("P01375")
    STATE["rcsb_down"] = False
    b = json.loads(path.read_text(encoding="utf-8"))["structures"][0]["interfaces"]
    assert r["ok"] and b["sequences"] is None and b["rcsb_assemblies"] is None and b["contacts"]
    assert any("RCSB assemblies and sequences unknown" in w for w in r["warnings"])

    # prepare_structures rewrites the manifest during the analysis: nothing written, said
    def rewrite():
        mm = json.loads(path.read_text(encoding="utf-8"))
        mm["created_at"] = "2026-10-29T00:00:00Z"
        choice.write_manifest(path, mm)
        STATE["rewrite"] = None
    STATE["rewrite"] = rewrite
    assert "rewrote the manifest" in run("P01375")["error"]

    # the raw file changed after the manifest: skipped, never analysed
    raw = path.parent / "raw" / "3WD5.cif"
    raw.write_bytes(raw.read_bytes() + b"\n")
    r = run("P01375")
    assert not r["ok"] and "SHA-256" in r["skipped"][0]["reason"]


def test_2w2n():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    path = make_manifest("2W2N")
    r = run("Q8NBP7")
    assert r["ok"], r
    b = json.loads(path.read_text(encoding="utf-8"))["structures"][0]["interfaces"]
    by = {(c["target_chain"], c["partner_chain"], c["symmetry"]["code"]): c for c in b["contacts"]}
    # PCSK9 in two target entities (A catalytic, P prodomain): their contact is target-target
    ap = by[("A", "P", "1_555")]
    assert ap["partner_role"] == "target" and ap["n_target_residues"] == 40 and ap["pisa_interface"]["id"] == 1
    # lattice shift folded into the operator as PISA writes it (-Y+1/2,X-1/2,Z-1/4)
    shifted = by[("A", "P", "2_544")]
    assert shifted["symmetry"] == {"code": "2_544", "operator": "-y+1/2,x-1/2,z-1/4", "lattice_shift": [0, -1, -1]}
    assert shifted["pisa_interface"]["id"] == 4
    ae = by[("A", "E", "1_555")]
    assert ae["partner_role"] == "partner" and ae["pisa_interface"]["target_residues_in_pisa"] == "15 of 16"
    # matching needs the same operator (reviewer 2026-10-28, 4Y6O: a 1-residue lattice contact got css 1 by
    # sharing residues with the large interface). The same residues at another operator match nothing.
    pis = pisa.parse_interfaces(PISA[("interfaces", "2w2n")])["interfaces"]
    moved = dict(shifted, symmetry={"code": "1_655", "operator": "x+1,y,z", "lattice_shift": [1, 0, 0]})
    assert interfaces.match_pisa(moved, pis) is None
    assert interfaces.match_pisa(dict(shifted, symmetry=dict(ap["symmetry"])), pis)["id"] == 1     # by operator, not residues
    # every contact found here has its PISA interface at the same operator
    assert all(c["pisa_interface"] for c in b["contacts"]), [c["symmetry"] for c in b["contacts"] if not c["pisa_interface"]]
    flags = {s["entity_id"]: s["flags"] for s in b["sequences"]}
    assert any(f.startswith("6 consecutive histidines at residues 3-8") for f in flags["2"])
    assert any(f.startswith("6 consecutive histidines at residues 307-312") for f in flags["1"])
    lld = next(s for s in b["sequences"] if s["entity_id"] == "2")
    assert lld["role"] == "partner" and lld["description"] == "LOW-DENSITY LIPOPROTEIN RECEPTOR" and lld["uniprot_ids"] == ["P01130"]


def test_carry_over():
    sha = "a" * 64
    old = {"structures": [{"pdb_id": "3WD5", "interfaces": {"raw_sha256": sha}},
                          {"pdb_id": "1TNF", "interfaces": {"raw_sha256": "b" * 64}},
                          {"pdb_id": "2AZ5", "interfaces": {"raw_sha256": "c" * 64}}]}
    new = {"structures": [{"pdb_id": "3WD5", "file": {"ok": True, "sha256": sha}},
                          {"pdb_id": "1TNF", "file": {"ok": True, "sha256": "d" * 64}}]}
    st = interfaces.carry_over(old, new)
    assert st == {"kept": ["3WD5"], "dropped": ["1TNF: structure file changed or not downloaded", "2AZ5: not in the new manifest"]}
    assert new["structures"][0]["interfaces"]["raw_sha256"] == sha and "interfaces" not in new["structures"][1]
    assert interfaces.carry_over(None, {"structures": []}) == {"kept": [], "dropped": []}


def main():
    test_parsers()
    test_3wd5()
    test_2w2n()
    test_carry_over()
    print("ALL OK")


main()
