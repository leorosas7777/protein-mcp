"""Live smoke test: calls the REAL APIs through server.py and checks invariants.

Two kinds of checks:
- EXTERNAL: facts that do not come from this code (approved drugs, known
  diseases). They catch a pipeline that runs fine but tells the analyst
  something false, which offline fixtures cannot.
- REGRESSION: real defects already fixed once (9E3S G12N, salts, antibodies).
Plus anti-hallucination checks on the output: unnamed things stay unnamed,
every link in the report is backed by the dossier.

Dossiers go to a temp folder; dossiers/ is never touched.
Run:  uv run python tests/smoke_live.py              (all, ~2-3 min)
      uv run python tests/smoke_live.py --skip-slow  (without EGFR, ~1 min)
      uv run python tests/smoke_live.py kras pcsk9   (only these groups)
Groups: kras pcsk9 egfr notfound tnf prep (Phase 2 downloads, BRAF V600E).
Exit code 1 if any check fails.
"""
import asyncio
import json
import logging
import re
import sys
import tempfile
import time
from pathlib import Path

from protein_mcp import dossier
from protein_mcp import server

logging.basicConfig(level=logging.WARNING, stream=sys.stderr)
logging.getLogger("httpx").setLevel(logging.WARNING)

RESULTS = []   # (group, kind, name, ok, detail)
GROUP = ""


def check(kind, name, ok, detail=""):
    RESULTS.append((GROUP, kind, name, bool(ok), "" if ok else str(detail)[:300]))


def names(items):   # non-dict items (old bare-id format) count as unnamed
    return {x.get("name") for x in items or [] if isinstance(x, dict)}


def by_name(items, name):
    return next((x for x in items or [] if isinstance(x, dict) and x.get("name") == name), {})


def links_backed_by_dossier(md, folder):
    """Every URL in the report is in the dossier, or built from an id that is."""
    blob = "".join(p.read_text(encoding="utf-8") for p in folder.glob("*.json"))
    bad = []
    for url in re.findall(r"\]\((https?://[^)\s]+)\)", md):
        if url in blob:
            continue
        m = re.fullmatch(r"https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/", url) \
            or re.fullmatch(r"https://www\.ebi\.ac\.uk/chembl/target_report_card/(CHEMBL\d+)/", url)
        if not (m and f'"{m.group(1)}"' in blob):
            bad.append(url)
    return bad


# ------------------------------------------------------------------ groups
async def kras():
    r = await server.resolve_protein("KRAS")
    check("EXTERNAL", "resolve KRAS -> P01116 first",
          r["candidates"] and r["candidates"][0]["uniprot_id"] == "P01116", r.get("candidates"))

    lig = await server.get_ligands("P01116", mutation="G12C")
    ligs = lig.get("ligands") or []
    for drug in ("SOTORASIB", "ADAGRASIB"):
        check("EXTERNAL", f"G12C ligands include approved {drug}",
              (by_name(ligs, drug).get("max_phase") or 0) >= 4, [(l["name"], l["max_phase"]) for l in ligs[:6]])
    check("REGRESSION", "G12C label/text conflicts are flagged (23 real rows labelled G12C, text G12D)",
          (lig.get("activities") or {}).get("excluded", {}).get("mutation_label_conflicts_description", 0) > 0,
          lig.get("activities"))
    check("ANTI-HALLUCINATION", "unnamed research compounds stay None",
          any(l.get("name") is None for l in ligs), "every ligand has a name: names may be invented")
    check("ANTI-HALLUCINATION", "each ligand link is built from its own ChEMBL id",
          all(l["chembl_id"] in (l.get("link") or "") for l in ligs))

    st = await server.get_structures("P01116", mutation="G12D", max_models=4)
    exp = st.get("experimental") or {}
    returned = exp.get("structures") or []
    check("REGRESSION", "G12D filter: every returned entry names G12D",
          returned and all("G12D" in (s.get("mutations") or []) + (s.get("mutations_from_title") or [])
                           for s in returned), [(s.get("pdb_id"), s.get("mutations")) for s in returned])
    check("REGRESSION", "9E3S (G12N) is never returned for G12D",
          all(s.get("pdb_id") != "9E3S" for s in returned))
    af = st.get("alphafold") or {}
    check("REGRESSION", "AlphaFold keeps the canonical model (not isoform P01116-2)",
          af.get("is_canonical_sequence") is True, af.get("entry_id"))

    pr = await server.get_protein_profile("P01116", max_diseases=5)
    dis = ((pr.get("uniprot") or {}).get("diseases") or {}).get("items") or []
    check("EXTERNAL", "UniProt diseases include Noonan syndrome 3", "Noonan syndrome 3" in names(dis), names(dis))
    check("ANTI-HALLUCINATION", "diseases without a UniProt name stay None and carry the note",
          all(d.get("note") for d in dis if d.get("name") is None)
          and sum(d.get("name") is None for d in dis) >= 1, dis)
    ot = (pr.get("open_targets") or {})
    onc = ((ot.get("associations") or {}).get("top_oncology") or [{}])[0]
    check("EXTERNAL", "Open Targets top oncology = non-small cell lung carcinoma",
          onc.get("name") == "non-small cell lung carcinoma", onc.get("name"))
    check("EXTERNAL", "Open Targets: KRAS small molecule has an approved drug",
          "Approved Drug" in ((ot.get("tractability") or {}).get("small molecule") or {}).get("yes", []))

    lit = await server.search_literature("P01116", max_results=5)
    papers = lit.get("papers") or []
    check("REGRESSION", "literature returns papers", len(papers) > 0, lit.get("error"))
    bad = []
    for p in papers:
        for url in re.findall(r"https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/", json.dumps(p)):
            if url != p.get("pmid"):
                bad.append((p.get("pmid"), url))
    check("ANTI-HALLUCINATION", "every PubMed link uses the paper's own PMID", not bad, bad)

    rep = await server.build_report("P01116")
    if rep.get("ok"):
        md = Path(rep["markdown_path"]).read_text(encoding="utf-8")
        bad = links_backed_by_dossier(md, Path(rep["markdown_path"]).parent)
        check("ANTI-HALLUCINATION", "every report link is backed by the dossier", not bad, bad[:5])
    else:
        check("REGRESSION", "build_report works", False, rep.get("error"))


async def pcsk9():
    lig = await server.get_ligands("Q8NBP7")
    mech = lig.get("mechanism_compounds_without_qualifying_activity") or []
    for drug in ("EVOLOCUMAB", "ALIROCUMAB"):
        x = by_name(mech, drug)
        check("EXTERNAL", f"PCSK9: {drug} listed as an approved antibody",
              x.get("molecule_type") == "Antibody" and (x.get("max_phase") or 0) >= 4, mech[:3])
    check("REGRESSION", "PCSK9: mechanism compounds are objects, not bare ids",
          mech and all(isinstance(x, dict) for x in mech), mech[:3])
    pr = await server.get_protein_profile("Q8NBP7", max_diseases=5)
    dis = names(((pr.get("uniprot") or {}).get("diseases") or {}).get("items"))
    check("EXTERNAL", "PCSK9: UniProt disease is familial hypercholesterolemia",
          any("ypercholesterolemia" in (n or "") for n in dis), dis)


async def egfr():
    t = time.time()
    lig = await server.get_ligands("P00533")
    a = lig.get("activities") or {}
    check("REGRESSION", "EGFR: truncation declared with its pChEMBL cutoff",
          a.get("truncated") and a.get("ranked_download_cutoff_pchembl") is not None, a)
    ligs = lig.get("ligands") or []
    for drug in ("OSIMERTINIB", "ERLOTINIB", "GEFITINIB"):
        check("EXTERNAL", f"EGFR: approved {drug} is ranked", (by_name(ligs, drug).get("max_phase") or 0) >= 4,
              sorted(n for n in names(ligs) if n))
    check("REGRESSION", "EGFR: erlotinib found via its parent id (mechanism recorded on the salt)",
          by_name(ligs, "ERLOTINIB").get("chembl_id") == "CHEMBL553")
    mech = names(lig.get("mechanism_compounds_without_qualifying_activity"))
    for drug in ("CETUXIMAB", "PANITUMUMAB"):
        check("EXTERNAL", f"EGFR: antibody {drug} listed", drug in mech, mech)
    conf = lig.get("confidence") or {}
    check("REGRESSION", "EGFR: every assay record read (5966 > old 5000 cap; a third went unfiltered)",
          conf.get("assay_records_total") and conf.get("assay_records_read") == conf.get("assay_records_total"), conf)
    check("INFO", f"EGFR get_ligands time {time.time() - t:.0f}s, chat {len(json.dumps(lig)) / 1024:.0f} KB", True)

    mut = await server.get_ligands("P00533", mutation="C797S")
    other = names(mut.get("mechanism_compounds_other_variants_only"))
    none_ = names(mut.get("mechanism_compounds_without_qualifying_activity"))
    check("REGRESSION", "EGFR C797S: lapatinib (IC50 on other variants) is 'other variants only', not 'no IC50'",
          "LAPATINIB" in other and "LAPATINIB" not in none_, {"other": sorted(n for n in other if n),
                                                              "none": sorted(n for n in none_ if n)})
    check("EXTERNAL", "EGFR C797S: antibody cetuximab still listed as without IC50/Ki/Kd", "CETUXIMAB" in none_, none_)


async def notfound():
    pr = await server.get_protein_profile("P0DZZ9")
    check("REGRESSION", "nonexistent accession -> clear error in get_protein_profile",
          pr.get("ok") is False and "not found" in (pr.get("error") or ""), pr)
    pr = await server.get_protein_profile("Q0ZZZ9")
    check("REGRESSION", "inactive accession (HTTP 200 Inactive) -> error",
          pr.get("ok") is False, pr)
    mouse = await server.get_protein_profile("P32883")
    check("ANTI-HALLUCINATION", "mouse Kras: Open Targets says 'sin registro' (human only)",
          (mouse.get("open_targets") or {}).get("found") is False, mouse.get("open_targets"))


async def tnf():
    pr = await server.get_protein_profile("P01375")
    ot = pr.get("open_targets") or {}
    check("REGRESSION", "TNF: 7 Ensembl genes in UniProt, Open Targets found on a later one (not 'sin registro')",
          ot.get("found") is True and ot.get("symbol") == "TNF" and ot.get("uniprot_crosscheck") is True,
          {k: ot.get(k) for k in ("found", "ensembl_gene_id", "ensembl_genes_without_target")})
    check("EXTERNAL", "TNF: at least one gene was skipped before the one with a target",
          len(ot.get("ensembl_genes_without_target") or []) >= 1, ot.get("ensembl_genes_without_target"))


async def prep():
    """Phase 2 steps 1-2: prepare_structures on BRAF V600E (defects found in the 8.2 validation and Desktop tests)."""
    st = await server.get_structures("P15056", mutation="V600E")
    have = {s["pdb_id"] for s in (st.get("experimental") or {}).get("structures") or []}
    wanted = ["3OG7", "5JSM", "7ZR0", "8C7Y"]
    check("EXTERNAL", "BRAF V600E structures include 3OG7, 5JSM, 7ZR0 and 8C7Y (needed by this group)",
          set(wanted) <= have, sorted(have))
    rep = await server.prepare_structures("P15056", mutation="V600E", pdb_ids=wanted + ["1UWH"])
    check("REGRESSION", "prepare_structures works", rep.get("ok"), rep.get("error"))
    if not rep.get("ok"):
        return
    check("ANTI-HALLUCINATION", "an ID outside the dossier file is not downloaded",
          rep["pdb_ids_not_in_dossier"] == ["1UWH"] and not any(s["pdb_id"] == "1UWH" for s in rep["structures"]),
          rep["pdb_ids_not_in_dossier"])
    base = Path(rep["prep_dir"])
    m = json.loads(Path(rep["manifest_path"]).read_text(encoding="utf-8"))
    by_id = {s["pdb_id"]: s for s in m["structures"]}
    o = by_id.get("3OG7") or {}
    chains = [c["auth_asym_id"] for t in o.get("target_entities") or [] for c in t["chains"]]
    check("REGRESSION", "3OG7: target chains A, B from RCSB (its mmCIF _struct_ref cites Q5IBP5)",
          chains == ["A", "B"], chains)
    lig = [i for i in o.get("nonpolymer_instances") or [] if i["comp_id"] == "032"]
    check("REGRESSION", "3OG7: one vemurafenib (032) instance, on chain A only",
          len(lig) == 1 and [v["auth_asym_id"] for v in lig[0]["neighbor_target_chains"].values()] == ["A"], lig)
    z = by_id.get("7ZR0") or {}
    zc = [c for t in z.get("target_entities") or [] for c in t["chains"]]
    check("REGRESSION", "7ZR0: BRAF chain is label D, auth K", zc == [{"label_asym_id": "D", "auth_asym_id": "K"}], zc)
    atp = [i for i in z.get("nonpolymer_instances") or [] if i["comp_id"] == "ATP"]
    check("REGRESSION", "7ZR0: ATP touches HSP90 only, no BRAF chain",
          atp and all(not i["neighbor_target_chains"] and i["neighbor_other_chains"] for i in atp), atp)
    six = [i for i in (by_id.get("5JSM") or {}).get("nonpolymer_instances") or [] if i["comp_id"] == "6NB"]
    poses = {(base / i["sdf"]["path"]).read_text().splitlines()[4] for i in six if (i["sdf"] or {}).get("ok")}
    check("REGRESSION", "5JSM: two 6NB instances, two SDFs with different poses", len(six) == 2 and len(poses) == 2,
          [i["label_asym_id"] for i in six])
    recs = [m["sequence"]] + ([m["alphafold"]["file"]] if m["alphafold"] else []) + [s["file"] for s in m["structures"]]         + [i["sdf"] for s in m["structures"] for i in s["nonpolymer_instances"] or [] if i["sdf"]]
    import hashlib
    bad = [r.get("path") or r.get("url") for r in recs if not r.get("ok")
           or hashlib.sha256((base / r["path"]).read_bytes()).hexdigest() != r["sha256"]]
    check("ANTI-HALLUCINATION", "every file in the manifest exists with its SHA-256", not bad, bad)
    src = json.dumps(st)
    urls = re.findall(r'https?://[^"\s]+', json.dumps(m))
    unbacked = [u for u in urls if u not in src and not any(i in u for i in by_id) and "P15056" not in u]
    check("ANTI-HALLUCINATION", "every manifest URL comes from the dossier or a returned ID", not unbacked, unbacked[:5])

    # step 2: cleaning, chain checks, pocket (gemmi)
    c3 = next((s for s in rep["structures"] if s["pdb_id"] == "3OG7"), {})
    checks = [(c["auth_asym_id"], c.get("variant")) for c in c3.get("chain_checks") or []]
    check("REGRESSION", "3OG7: V600E not modelled in chain A (ligand chain), GLU 600 in chain B",
          checks == [("A", "not modelled"), ("B", "variant")], checks)
    pk = (c3.get("pockets") or [{}])[0]
    check("REGRESSION", "3OG7: vemurafenib pocket = 22 residues incl. K483 and DFG, next to missing 597-614",
          pk.get("n_residues") == 22 and all(f"A:{r}" in pk.get("residues", "") for r in ("K483", "D594", "F595", "G596"))
          and {"auth_asym_id": "A", "missing_uniprot_range": [597, 614]} in pk.get("adjacent_missing_ranges", []), pk)
    rec = base / (c3.get("receptor") or "missing")
    het = {line[17:20].strip() for line in rec.read_text().splitlines() if line.startswith("HETATM")} if rec.is_file() else None
    check("REGRESSION", "3OG7: receptor has chains only, no ligand or water", het == set(), het)
    check("EXTERNAL", "AlphaFold is wild type (VAL) at 600",
          ((rep.get("alphafold") or {}).get("variant_position") or {}).get("residue") == "VAL", rep.get("alphafold"))

    # step 2.1 (Desktop test 2026-10-26): ligand names, bivalent ligand, failures in the manifest
    vem = next((x for x in c3.get("ligand_instances") or [] if x["comp_id"] == "032"), {})
    check("EXTERNAL", "3OG7: RCSB (wwPDB) names 032 Vemurafenib / PLX4032, no DrugBank synonyms",
          vem.get("synonyms") == ["Vemurafenib", "PLX4032"], vem)
    c5 = next((s for s in rep["structures"] if s["pdb_id"] == "5JSM"), {})
    p5 = [(p.get("residues_by_chain"), p.get("spans_target_chains")) for p in c5.get("pockets") or []]
    check("REGRESSION", "5JSM: each linked-vemurafenib 6NB pocket spans two chains (24 + 23 residues)",
          p5 == [({"A": 24, "B": 23}, True), ({"C": 24, "D": 23}, True)] and "VEM-3-VEM" in (c5.get("title") or ""), p5)
    check("ANTI-HALLUCINATION", "manifest lists the same failed downloads as the reply",
          m.get("failed_downloads") == rep.get("failed_downloads"), m.get("failed_downloads"))
    sm = rep.get("summary") or {}
    pv = sm.get("pockets_with_variant_in_pocket_chain") or []
    check("ANTI-HALLUCINATION", "summary: V600E modelled in the pocket chain of 8C7Y TXV (A) and 5JSM 6NB, not 3OG7 032",
          "8C7Y TXV U: chain A variant (26 pocket residues); single-chain pocket" in pv
          and any(x.startswith("5JSM 6NB J: chain A variant") for x in pv)
          and any(x.startswith("3OG7 032 C: chain A not modelled") for x in sm.get("pockets_without_variant_in_pocket_chain", [])),
          pv)
    check("ANTI-HALLUCINATION", "summary: 8C7Y chain A has no internal gap; 8C7Y is first by resolution and R-free",
          any(c.startswith("8C7Y:A ") for c in sm.get("chains_without_internal_gaps", [])) and (sm.get("by_resolution") or [""])[0].startswith("8C7Y")
          and (sm.get("by_r_free") or [""])[0].startswith("8C7Y"), {k: sm.get(k) for k in ("chains_without_internal_gaps", "by_resolution", "by_r_free")})

    ks = await server.get_structures("P01116", mutation="G12C")
    nuc = [s["pdb_id"] for s in (ks.get("experimental") or {}).get("structures") or [] if "GDP" in s.get("nucleotides", [])][:1]
    kr = await server.prepare_structures("P01116", mutation="G12C", pdb_ids=nuc, include_alphafold=False) if nuc else {}
    kst = (kr.get("structures") or [{}])[0]
    krec = Path(kr["prep_dir"]) / kst["receptor"] if kst.get("receptor") else None
    khet = {line[17:20].strip() for line in krec.read_text().splitlines() if line.startswith("HETATM")} if krec else None
    check("REGRESSION", f"KRAS G12C {nuc}: receptor keeps GDP and its Mg2+, CYS 12 modelled",
          khet is not None and {"GDP", "MG"} <= khet and [c.get("variant") for c in kst.get("chain_checks", [])][:1] == ["variant"],
          {"het": khet, "checks": kst.get("chain_checks"), "error": kr.get("error")})


GROUPS = {"kras": kras, "pcsk9": pcsk9, "egfr": egfr, "notfound": notfound, "tnf": tnf, "prep": prep}
SLOW = {"egfr"}


async def main(argv):
    global GROUP
    chosen = [a for a in argv if a in GROUPS] or [g for g in GROUPS
                                                  if not ("--skip-slow" in argv and g in SLOW)]
    dossier.DOSSIER_DIR = Path(tempfile.mkdtemp(prefix="smoke_"))
    for g in chosen:
        GROUP = g
        try:
            await GROUPS[g]()
        except Exception as e:   # an API down must show as a failure, not stop the run
            check("ERROR", f"{g} crashed", False, f"{type(e).__name__}: {e}")
    fails = [r for r in RESULTS if not r[3]]
    for g, kind, name, ok, detail in RESULTS:
        print(f"{'PASS' if ok else 'FAIL'}  [{g}/{kind}] {name}" + (f"\n      -> {detail}" if detail else ""))
    print(f"\n{len(RESULTS) - len(fails)} passed, {len(fails)} failed. Dossiers: {dossier.DOSSIER_DIR}")
    print("ALL OK" if not fails else "FAILURES")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main(sys.argv[1:])))
