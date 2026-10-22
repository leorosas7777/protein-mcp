"""Offline test for report.py. The dossier is built from REAL data: the profile from the
UniProt/Open Targets fixtures (through target_profile + attach_dossier, as server.py
does) and three real dossier files saved on 2026-10-20 (structures G12D, ligands G12C,
literature 'inhibitor'; one literature title contains 'CD8⁺', which Arial lacks).
Key check (anti-hallucination): every link in the report exists in the dossier or is
built from an ID that does.
Run: uv run python tests/test_report.py
"""
import asyncio
import json
import re
import shutil
import tempfile
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parent.parent  # fixtures live in tests/fixtures
from protein_mcp import dossier
from protein_mcp import report
from protein_mcp import target_profile

FIX = ROOT / "tests" / "fixtures"
ACC = "P01116"


async def save_profile():
    uniprot = json.loads((FIX / "uniprot_P01116_profile.json").read_text(encoding="utf-8"))
    ot = json.loads((FIX / "opentargets_ENSG00000133703.json").read_text(encoding="utf-8"))

    def handler(req):
        return httpx.Response(200, json=uniprot if req.url.host == "rest.uniprot.org" else ot)

    async with httpx.AsyncClient(transport=httpx.MockTransport(handler)) as c:
        full = await target_profile.get_protein_profile(ACC, _client=c)
    dossier.attach_dossier(full, ACC, "profile", "get_protein_profile", {"max_diseases": 10})


def links_are_sourced(md: str, folder: Path):
    blob = "".join(p.read_text(encoding="utf-8") for p in folder.glob("*.json") if not p.name.startswith("report"))
    for url in re.findall(r"\]\((https?://[^)\s]+)\)", md):
        if url in blob:
            continue
        m = re.fullmatch(r"https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/", url) \
            or re.fullmatch(r"https://www\.ebi\.ac\.uk/chembl/target_report_card/(CHEMBL\d+)/", url)
        assert m and f'"{m.group(1)}"' in blob, f"link not backed by the dossier: {url}"


def main():
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())
    folder = dossier.DOSSIER_DIR / ACC

    # no dossier -> clear error, nothing written
    r = report.build_report(ACC)
    assert r["ok"] is False and "No dossier" in r["error"]
    assert report.build_report("not an id")["ok"] is False

    # only the profile -> report with the other sections listed as missing
    asyncio.run(save_profile())
    r = report.build_report(ACC)
    assert r["ok"] and set(r["missing"]) == {"get_structures", "get_ligands", "search_literature"}, r
    md = Path(r["markdown_path"]).read_text(encoding="utf-8")
    assert 'Run get_ligands("P01116")' in md and "## 4. Ligands" not in md

    # full dossier
    for f in (FIX / "dossier_P01116").glob("*.json"):
        shutil.copy(f, folder / f.name)
    (folder / "notes.json").write_text('{"not": "a dossier file"}', encoding="utf-8")   # ignored
    r = report.build_report(ACC)
    assert r["ok"] and r["missing"] == [] and r["pages"] >= 4, r
    assert r["included"]["get_ligands"] == ["ligands_g12c.json"]
    pdf = Path(r["pdf_path"]).read_bytes()
    assert pdf.startswith(b"%PDF") and len(pdf) > 20_000
    md = Path(r["markdown_path"]).read_text(encoding="utf-8")

    links_are_sourced(md, folder)
    # values quoted from the records
    assert "| G12C |" in md and "KIKK182-185→RIRR" in md and "G10→GG" in md
    assert md.count("Unnamed in UniProt. Note:") == 2                    # never given a name
    assert "[Noonan syndrome]" in md and "2262" in md and "[non-small cell lung carcinoma]" in md
    assert "ADAGRASIB" in md and "4 (approved 2022)" in md and "no name in ChEMBL" in md
    assert "4 of 127 entries shown" in md and "9E3S" in md                # selection + discarded
    assert "5 of " in md and "CD8⁺" in md                                 # Markdown keeps the original
    assert "variant filter 'G12D'" in md and "AF-P01116-F1" in md
    # the PDF font cannot draw '⁺': it is mapped, not dropped
    assert report._printable("CD8⁺ T") == "CD8+ T"

    # an open PDF (Windows lock) -> clear error instead of a crash
    real = report.render_pdf
    report.render_pdf = lambda *a, **k: (_ for _ in ()).throw(PermissionError(13, "Permission denied"))
    try:
        r = report.build_report(ACC)
    finally:
        report.render_pdf = real
    assert r["ok"] is False and "close it" in r["error"]

    # mechanism compounds without IC50/Ki/Kd (real PCSK9 shape): a table with name and
    # phase; dossiers written before 2026-10-21 hold bare ids and still render
    lig_file = folder / "ligands_g12c.json"
    saved = json.loads(lig_file.read_text(encoding="utf-8"))
    saved["result"]["mechanism_compounds_without_qualifying_activity"] = [{
        "chembl_id": "CHEMBL2364655", "name": "EVOLOCUMAB",
        "link": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL2364655/",
        "molecule_type": "Antibody", "max_phase": 4.0, "first_approval": 2015,
        "mechanism_of_action": "Subtilisin/kexin type 9 inhibitor", "action_type": "INHIBITOR"}]
    lig_file.write_text(json.dumps(saved), encoding="utf-8")
    md2 = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    assert "EVOLOCUMAB" in md2 and "| Antibody | 4 (approved 2015) |" in md2, md2[-1500:]
    links_are_sourced(md2, folder)
    saved["result"]["mechanism_compounds_without_qualifying_activity"] = ["CHEMBL2364655"]
    lig_file.write_text(json.dumps(saved), encoding="utf-8")
    md3 = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    assert "no qualifying activity: CHEMBL2364655" in md3

    # other-variants-only and unchecked compounds get their own tables; a failed
    # molecule lookup is "not retrieved"; ChEMBL max_phase -1 is spelled out
    base = {"link": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL553/", "chembl_id": "CHEMBL553",
            "molecule_type": "Small molecule", "first_approval": None, "action_type": "INHIBITOR",
            "mechanism_of_action": "Epidermal growth factor receptor erbB1 inhibitor"}
    saved["result"]["mechanism_compounds_without_qualifying_activity"] = []
    saved["result"]["mechanism_compounds_other_variants_only"] = [
        {**base, "name": "ERLOTINIB", "molecule_record_retrieved": True, "max_phase": 4.0}]
    saved["result"]["mechanism_compounds_unchecked"] = [
        {**base, "name": None, "molecule_record_retrieved": False, "max_phase": -1.0}]
    saved["result"]["confidence"]["min_applied"] = None
    lig_file.write_text(json.dumps(saved), encoding="utf-8")
    md4 = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    assert "only for variants other than 'G12C' (1;" in md4 and "ERLOTINIB" in md4, md4[-2000:]
    assert "could not be checked because their activity download failed (1;" in md4
    assert "molecule record not retrieved" in md4 and "clinical, phase unknown (ChEMBL -1)" in md4
    assert "Minimum assay confidence: not applied." in md4 and "confidence: None" not in md4
    links_are_sourced(md4, folder)
    assert "mechanism lookup failed" not in md4   # old dossiers have no "mechanisms" block

    # failed mechanism lookup: the lists are None and section 4 warns (not only Methods)
    for field in ("mechanism_compounds_without_qualifying_activity", "mechanism_compounds_other_variants_only",
                  "mechanism_compounds_unchecked"):
        saved["result"][field] = None
    saved["result"]["mechanisms"] = {"retrieved": False, "records_read": None, "records_total": None}
    lig_file.write_text(json.dumps(saved), encoding="utf-8")
    md5 = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
    sec4 = md5.split("4. Ligands")[1].split("5. Literature")[0]
    assert "Warning: the ChEMBL mechanism lookup failed" in sec4, sec4[:1500]

    # truncated download: the Activities line follows mechanism_compounds_download,
    # never "plus every activity" when the extra download failed (unchecked may be empty)
    saved["result"]["mechanisms"] = {"retrieved": True, "records_read": 1, "records_total": 1}
    saved["result"]["mechanism_compounds_unchecked"] = []
    act_block = saved["result"]["activities"]
    act_block.update({"truncated": True, "ranked_download_cutoff_pchembl": 7.94})
    expected = {"failed": "the extra download for mechanism compounds failed (0 unchecked;",
                "complete": "plus every activity of mechanism compounds",
                "mechanism_lookup_failed": "mechanism lookup failed: no extra download for drugs"}
    for state, text in expected.items():
        act_block["mechanism_compounds_download"] = state
        lig_file.write_text(json.dumps(saved), encoding="utf-8")
        line = next(l for l in Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8")
                    .splitlines() if "fetched (truncated" in l)
        assert text in line, (state, line)
        assert state == "complete" or "plus every activity" not in line, (state, line)
    # several Ensembl genes (real TNF case): genes skipped are listed whether a later
    # gene was found, none was, or a query failed after them
    prof_file = folder / "profile.json"
    prof = json.loads(prof_file.read_text(encoding="utf-8"))
    ot = prof["result"]["open_targets"]
    variants = {
        "found": {**ot, "ensembl_genes_without_target": ["ENSG00000204490"],
                  "other_ensembl_genes_not_queried": ["ENSG00000232810"]},
        "none": {"ok": True, "found": False, "ensembl_gene_id": None, "data_version": "26.09",
                 "ensembl_genes_without_target": ["ENSG00000204490"],
                 "note": "sin registro en Open Targets para ninguno de los 1 genes Ensembl de UniProt."},
        "error": {"ok": False, "error": "Open Targets request failed: timeout",
                  "ensembl_genes_without_target": ["ENSG00000204490"]},
    }
    for case, value in variants.items():
        prof["result"]["open_targets"] = value
        prof_file.write_text(json.dumps(prof), encoding="utf-8")
        sec2 = Path(report.build_report(ACC)["markdown_path"]).read_text(encoding="utf-8") \
            .split("## 2.")[1].split("## 3.")[0]
        assert "ENSG00000204490" in sec2, (case, sec2[:800])
        assert case != "found" or ("not queried: ENSG00000232810" in sec2 and "[Noonan syndrome]" in sec2)
        assert case != "error" or "could not be queried" in sec2
    prof["result"]["open_targets"] = ot
    prof_file.write_text(json.dumps(prof), encoding="utf-8")

    # False/0 parameters are shown, so a non-default call never reads as the default one
    assert report._params_text({"topic": "", "include_preprints": False, "max_models": 0}) == \
        "include_preprints=False, max_models=0"
    assert report._params_text({"topic": ""}) == "defaults"
    print(json.dumps({"pages": report.build_report(ACC)["pages"], "md_kb": round(len(md) / 1024, 1)}))
    print("ALL OK")


main()
