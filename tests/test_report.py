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
import sys
import tempfile
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))  # project root
import common  # noqa: E402
import report  # noqa: E402
import target_profile  # noqa: E402

FIX = ROOT / "tests" / "fixtures"
ACC = "P01116"


async def save_profile():
    uniprot = json.loads((FIX / "uniprot_P01116_profile.json").read_text(encoding="utf-8"))
    ot = json.loads((FIX / "opentargets_ENSG00000133703.json").read_text(encoding="utf-8"))

    def handler(req):
        return httpx.Response(200, json=uniprot if req.url.host == "rest.uniprot.org" else ot)

    async with httpx.AsyncClient(transport=httpx.MockTransport(handler)) as c:
        full = await target_profile.get_protein_profile(ACC, _client=c)
    common.attach_dossier(full, ACC, "profile", "get_protein_profile", {"max_diseases": 10})


def links_are_sourced(md: str, folder: Path):
    blob = "".join(p.read_text(encoding="utf-8") for p in folder.glob("*.json") if not p.name.startswith("report"))
    for url in re.findall(r"\]\((https?://[^)\s]+)\)", md):
        if url in blob:
            continue
        m = re.fullmatch(r"https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/", url) \
            or re.fullmatch(r"https://www\.ebi\.ac\.uk/chembl/target_report_card/(CHEMBL\d+)/", url)
        assert m and f'"{m.group(1)}"' in blob, f"link not backed by the dossier: {url}"


def main():
    common.DOSSIER_DIR = Path(tempfile.mkdtemp())
    folder = common.DOSSIER_DIR / ACC

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
    print(json.dumps({"pages": report.build_report(ACC)["pages"], "md_kb": round(len(md) / 1024, 1)}))
    print("ALL OK")


main()
