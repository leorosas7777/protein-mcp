"""build_report: dossier on disk -> Markdown + PDF report.

Option (a), decided 2026-10-20: everything in the report is read from
dossiers/<accession>/*.json (written by the tools through dossier.attach_dossier)
and laid out deterministically. Nothing is generated or inferred here: values
and links come from the saved source records, empty fields say "no record",
and the only fixed text describes how the data were selected. The analyst's
interpretation is added by Claude in the chat, not in this document.

Imports dossier, common, files.choice (manifest path and choice label) and reportlab: it never
calls the source APIs. Section 6 reads prep/manifest.json (Phase 2) when it exists.
"""
import json
import os
import unicodedata
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (ListFlowable, ListItem, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

from protein_mcp import dossier
from protein_mcp.files import choice
from protein_mcp.common import log, parse_accession, utc_now

NO_RECORD = "no record"
LINK_COLOR = "#1a4fa0"

# Order of the data sections and the tool that feeds each one.
SECTIONS = [
    ("get_protein_profile", "Target profile (UniProt) and target-disease evidence (Open Targets)"),
    ("get_structures", "Structures (RCSB PDB, AlphaFold DB)"),
    ("get_ligands", "Ligands and bioactivity (ChEMBL)"),
    ("search_literature", "Literature (Europe PMC)"),
]


# ------------------------------------------------------------ text model
class Link:
    """Inline hyperlink. `url` must come from the dossier (or be built from an ID in it)."""
    def __init__(self, text, url):
        self.text, self.url = str(text), url


class Lines(list):
    """Several lines inside one table cell."""


def _md(x) -> str:
    if x is None:
        return f"*{NO_RECORD}*"
    if isinstance(x, Link):
        text = x.text.replace("|", "\\|").replace("[", "(").replace("]", ")")
        return f"[{text}]({x.url})" if x.url else text
    if isinstance(x, Lines):
        return "<br>".join(_md(i) for i in x)
    if isinstance(x, (list, tuple)):
        return "".join(_md(i) for i in x)
    return str(x).replace("|", "\\|").replace("\n", " ")


_CMAP = None   # code points the PDF font can draw (None = built-in font, no check)


def _printable(text: str) -> str:
    """PDF only: a character the font lacks (e.g. superscript ⁺ in 'CD8⁺') is replaced by
    its compatibility form ('+'), or by '?' if that is missing too. Markdown keeps the original."""
    if not _CMAP:
        return text
    out = []
    for c in text:
        if ord(c) < 128 or ord(c) in _CMAP:
            out.append(c)
        else:
            alt = unicodedata.normalize("NFKC", c)
            out.append(alt if alt != c and all(ord(a) < 128 or ord(a) in _CMAP for a in alt) else "?")
    return "".join(out)


def _pdf(x) -> str:
    if x is None:
        return f"<i>{NO_RECORD}</i>"
    if isinstance(x, Link):
        text = escape(_printable(x.text))
        if not x.url:
            return text
        return f'<link href="{escape(x.url, {chr(34): "&quot;"})}" color="{LINK_COLOR}">{text}</link>'
    if isinstance(x, Lines):
        return "<br/>".join(_pdf(i) for i in x)
    if isinstance(x, (list, tuple)):
        return "".join(_pdf(i) for i in x)
    return escape(_printable(str(x)))


class Doc:
    """Ordered blocks rendered twice: to Markdown and to PDF."""
    def __init__(self):
        self.blocks = []

    def h(self, level, text):
        self.blocks.append(("h", level, text))

    def p(self, *parts, small=False):
        self.blocks.append(("p", list(parts), small))

    def bullets(self, items):
        if items:
            self.blocks.append(("ul", items))

    def table(self, headers, rows, widths):
        if rows:
            self.blocks.append(("table", headers, rows, widths))


# ------------------------------------------------------------ helpers
def _first(entries, tool):
    items = entries.get(tool) or []
    return items[0] if items else None


def _clip(text, n):
    if not text:
        return None
    text = str(text)
    return text if len(text) <= n else text[: n - 1] + "…"


def _num(x, fmt="{:g}"):
    return fmt.format(x) if isinstance(x, (int, float)) else None


def _pubmed(ids, limit=3):
    """Links built from PubMed IDs stored in the dossier."""
    ids = list(ids or [])
    out = []
    for i, pid in enumerate(ids[:limit]):
        if i:
            out.append(", ")
        out.append(Link(f"PMID {pid}", f"https://pubmed.ncbi.nlm.nih.gov/{pid}/"))
    if len(ids) > limit:
        out.append(f" (+{len(ids) - limit} more in dossier)")
    return out or None


def _change_label(v):
    """'G12C' for single-residue substitutions; otherwise built from the raw fields,
    e.g. 'KIKK182-185→RIRR' or 'G10→GG' (insertion)."""
    if v.get("change"):
        return v["change"]
    start, end = v.get("start"), v.get("end")
    span = f"{start}-{end}" if end not in (None, start) else f"{start}"
    return f"{v.get('original') or ''}{span}→{'/'.join(v.get('alternatives') or []) or '?'}"


def _params_text(params):
    # False and 0 are shown: hiding them made two different calls look the same.
    shown = [f"{k}={v!r}" for k, v in (params or {}).items() if v not in ("", None)]
    return ", ".join(shown) or "defaults"


def _phase_text(x):
    """ChEMBL max_phase as text (meanings from ChEMBL's schema documentation)."""
    phase = x.get("max_phase")
    if not isinstance(phase, (int, float)):
        return "none recorded"
    if phase < 0:
        return "clinical, phase unknown (ChEMBL -1)"
    text = "early phase 1 (0.5)" if phase == 0.5 else f"{phase:g}"
    return text + (f" (approved {x['first_approval']})" if x.get("first_approval") else "")


def _compound_cell(x):
    if x.get("molecule_record_retrieved") is False:
        label = "molecule record not retrieved"
    else:
        label = x.get("name") or "no name in ChEMBL"
    return Lines([Link(x.get("chembl_id"), x.get("link")), label])


def _paper_link(p):
    links = p.get("links") or {}
    return links.get("pubmed") or links.get("doi") or links.get("europepmc")


# ------------------------------------------------------------ dossier
def load_dossier(accession: str) -> dict:
    """{tool: [entries]} from dossiers/<accession>/*.json; the default call
    (shortest file name, e.g. ligands.json before ligands_g12c.json) first."""
    folder = dossier.DOSSIER_DIR / accession
    entries = {}
    for path in sorted(folder.glob("*.json")) if folder.is_dir() else []:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            log.warning("Skipping unreadable dossier file %s", path)
            continue
        if not isinstance(data, dict) or "tool" not in data or "result" not in data:
            continue
        data["file"] = path.name
        entries.setdefault(data["tool"], []).append(data)
    for items in entries.values():
        items.sort(key=lambda d: (len(d["file"]), d["file"]))
    return entries


# ------------------------------------------------------------ sections
def _coverage(doc, acc, entries, manifest=None):
    doc.h(2, "Contents of the dossier")
    rows = []
    for tool, label in SECTIONS:
        for d in entries.get(tool) or []:
            rows.append([label, _params_text(d.get("params")), d.get("retrieved_at"), d["file"]])
    doc.table(["Section", "Parameters", "Retrieved (UTC)", "File"], rows, (0.32, 0.25, 0.17, 0.26))
    missing = [f"{label}: not in the dossier. Run {tool}(\"{acc}\") and build the report again."
               for tool, label in SECTIONS if not entries.get(tool)]
    if missing:
        doc.p("Missing sections (left empty instead of filled from memory):")
        doc.bullets(missing)
    if manifest and not manifest.get("unreadable"):
        doc.p(f"Structure preparation: prep/manifest.json, created {manifest.get('created_at')} (section 6).", small=True)


def _uniprot(doc, d):
    up = d["result"]["uniprot"]
    doc.h(2, "1. Target profile (UniProt)")
    doc.p("Source: ", Link(f"UniProtKB {up['uniprot_id']}", up.get("source_url")),
          f" · {'reviewed (Swiss-Prot)' if up.get('reviewed') else 'unreviewed (TrEMBL)'}"
          f" · annotation score {_num(up.get('annotation_score')) or NO_RECORD}/5"
          f" · last annotation update {up.get('last_annotation_update') or NO_RECORD}"
          f" · retrieved {d.get('retrieved_at')}", small=True)

    doc.h(3, "Function")
    if not up.get("function"):
        doc.p(NO_RECORD)
    for f in up.get("function") or []:
        doc.p(f["text"])
        doc.p("Evidence: ", ", ".join(f.get("evidence") or []) or NO_RECORD, " · ",
              _pubmed(f.get("pubmed"), 6) or NO_RECORD, small=True)

    doc.h(3, "Subcellular location")
    locs = []
    for loc in up.get("subcellular_locations") or []:
        text = "; ".join(x for x in (loc.get("location"), loc.get("topology"), loc.get("orientation")) if x)
        locs.append([text + (f" ({loc['molecule']})" if loc.get("molecule") else ""), " · ",
                     _pubmed(loc.get("pubmed"), 3) or "no PubMed evidence listed"])
    doc.bullets(locs) if locs else doc.p(NO_RECORD)

    doc.h(3, f"Disease involvement curated by UniProt ({len(up.get('diseases') or [])})")
    rows = []
    for dis in up.get("diseases") or []:
        links = dis.get("links") or {}
        if dis.get("name"):
            name = Lines([[Link(dis["name"], links.get("uniprot_disease")),
                           f" ({dis['acronym']})" if dis.get("acronym") else ""]])
        else:  # UniProt gives only a note: quoted, never named
            name = Lines(["Unnamed in UniProt. Note: " + (_clip(dis.get("note"), 300) or NO_RECORD)])
        mim = Link(dis["mim"], links.get("omim")) if dis.get("mim") else None
        rows.append([name, mim, _pubmed(dis.get("pubmed"), 3)])
    doc.table(["Disease", "MIM", "Evidence"], rows, (0.55, 0.12, 0.33)) if rows else doc.p(NO_RECORD)

    variants = up.get("natural_variants") or []
    doc.h(3, f"Natural variants ({len(variants)})")
    rows = []
    for v in variants:
        change = _change_label(v)
        dbsnp = (v.get("links") or {}).get("dbsnp") or []
        rows.append([change, _clip(v.get("description"), 180),
                     Lines([Link(r, u) for r, u in zip(v.get("dbsnp") or [], dbsnp)]) or None,
                     _pubmed(v.get("pubmed"), 2)])
    doc.table(["Variant", "UniProt description", "dbSNP", "Evidence"], rows, (0.12, 0.50, 0.15, 0.23)) \
        if rows else doc.p(NO_RECORD)

    muts = up.get("mutagenesis") or []
    if muts:
        doc.h(3, f"Mutagenesis experiments ({len(muts)})")
        doc.table(["Mutation", "Effect reported in UniProt", "Evidence"],
                  [[_change_label(m),
                    _clip(m.get("description"), 200), _pubmed(m.get("pubmed"), 2)] for m in muts],
                  (0.14, 0.62, 0.24))

    xrefs = []
    if up.get("ensembl_gene_ids"):
        xrefs.append("Ensembl gene: " + ", ".join(up["ensembl_gene_ids"]))
    for cid in up.get("chembl_target_ids") or []:
        xrefs.append(["ChEMBL target: ", Link(cid, f"https://www.ebi.ac.uk/chembl/target_report_card/{cid}/")])
    if xrefs:
        doc.h(3, "Cross-references")
        doc.bullets(xrefs)


def _association_rows(rows):
    out = []
    for r in rows:
        scores = sorted((r.get("datatype_scores") or {}).items(), key=lambda kv: -(kv[1] or 0))
        out.append([Link(r.get("name") or r.get("disease_id"), (r.get("links") or {}).get("evidence")),
                    _num(r.get("score"), "{:.3f}"),
                    "; ".join(f"{k.replace('_', ' ')} {v:.2f}" for k, v in scores if v is not None) or None])
    return out


def _open_targets(doc, d):
    ot = d["result"]["open_targets"]
    acc = d["result"]["uniprot_id"]
    n = int((d.get("params") or {}).get("max_diseases") or 10)
    doc.h(2, "2. Target-disease evidence and tractability (Open Targets)")
    if not ot.get("ok") or not ot.get("found"):
        doc.p(ot.get("note") or NO_RECORD if ot.get("ok")
              else f"Open Targets could not be queried: {ot.get('error') or NO_RECORD}")
        # genes queried before the stop (no target, or a failed query) are still reported
        if ot.get("ensembl_genes_without_target"):
            doc.p("Ensembl genes queried without an Open Targets target: "
                  + ", ".join(ot["ensembl_genes_without_target"]), small=True)
        return
    doc.p("Source: ", Link(f"Open Targets Platform, {ot.get('symbol')} ({ot.get('ensembl_gene_id')})",
                           (ot.get("links") or {}).get("target")),
          f" · data release {ot.get('data_version') or NO_RECORD} · retrieved {d.get('retrieved_at')}", small=True)
    if not ot.get("uniprot_crosscheck"):
        doc.p(f"Warning: Open Targets maps this gene to {', '.join(ot.get('swissprot_ids') or []) or NO_RECORD}, "
              f"not to {acc}.")
    if ot.get("ensembl_genes_without_target"):
        doc.p("Other Ensembl genes listed by UniProt, queried first, with no Open Targets target: "
              + ", ".join(ot["ensembl_genes_without_target"]), small=True)
    if ot.get("other_ensembl_genes_not_queried"):
        doc.p("Other Ensembl genes listed by UniProt, not queried: " + ", ".join(ot["other_ensembl_genes_not_queried"]))

    doc.h(3, "Tractability flags")
    names = ot.get("modality_names") or {}
    rows = [[names.get(m, m), ", ".join(x["label"] for x in flags if x.get("value") is True) or "none",
             ", ".join(x["label"] for x in flags if x.get("value") is False) or "none"]
            for m, flags in (ot.get("tractability") or {}).items()]
    doc.table(["Modality", "Flag = yes", "Flag = no"], rows, (0.16, 0.42, 0.42)) if rows else doc.p(NO_RECORD)
    doc.p("Open Targets' own yes/no flags per modality; they are not a druggability verdict.", small=True)

    assoc = ot.get("associations") or {}
    rows = assoc.get("rows") or []
    doc.h(3, f"Top {min(n, len(rows))} associated diseases by overall score")
    doc.p(f"{assoc.get('examined')} of {assoc.get('total')} associations examined "
          f"({assoc.get('selection') or 'ordered by overall score'}). Scores (0-1) aggregate evidence; "
          "they are not proof of causality.", small=True)
    doc.table(["Disease (link to evidence)", "Score", "Evidence types (score)"],
              _association_rows(rows[:n]), (0.36, 0.09, 0.55))

    onc = [r for r in rows if r.get("is_oncology")]
    doc.h(3, f"Top {min(n, len(onc))} oncology associations")
    doc.p(f"{len(onc)} of the {len(rows)} examined associations belong to Open Targets' therapeutic area "
          "'cancer or benign tumor' (its own classification, which includes some syndromes).", small=True)
    doc.table(["Disease (link to evidence)", "Score", "Evidence types (score)"],
              _association_rows(onc[:n]), (0.36, 0.09, 0.55))


def _mutation_cell(s):
    parts = []
    if s.get("mutations"):
        parts.append("field: " + ", ".join(s["mutations"]))
    if s.get("mutations_from_title"):
        parts.append("title: " + ", ".join(s["mutations_from_title"]))
    if s.get("mutations_unparsed"):
        parts.append("other text: " + ", ".join(s["mutations_unparsed"]))
    return Lines(parts) if parts else "none listed (not proof of wild type)"


def _structures(doc, items):
    doc.h(2, "3. Structures")
    for d in items:
        e = d["result"].get("experimental") or {}
        mut = (d.get("params") or {}).get("mutation")
        doc.h(3, "Experimental structures (RCSB PDB)" + (f", variant filter '{mut}'" if mut else ""))
        if not e.get("ok"):
            doc.p(f"RCSB could not be queried: {e.get('error') or NO_RECORD}")
            continue
        doc.p(f"{e.get('returned')} of {e.get('total_entries')} entries shown. {e.get('selection') or ''} "
              f"Retrieved {d.get('retrieved_at')}.", small=True)
        rows = []
        for s in e.get("structures") or []:
            links = s.get("links") or {}
            ligands = [f"{x.get('id')} ({x.get('name')})" for x in s.get("ligands") or []]
            ligands += [f"partner chain: {p}" for p in s.get("partner_chains") or []]
            cite = [Link(k.upper() if k == "doi" else "PubMed", links[k]) for k in ("pubmed", "doi") if links.get(k)]
            rows.append([
                Link(s.get("pdb_id"), links.get("rcsb")),
                _clip(s.get("title"), 110),
                Lines([s.get("method") or NO_RECORD,
                       f"{_num(s.get('resolution_angstrom')) or 'no'} Å" if s.get("resolution_angstrom") else "resolution: no record",
                       s.get("release_date") or ""]),
                _mutation_cell(s),
                Lines(ligands) if ligands else "none",
                ", ".join(s.get("nucleotides") or []) or "none",
                Lines(cite) if cite else None,
            ])
        doc.table(["PDB", "Title", "Method", "Mutations", "Ligands", "Nucl.", "Citation"], rows,
                  (0.07, 0.27, 0.14, 0.14, 0.22, 0.06, 0.10))
        if e.get("discarded_not_confirmed"):
            doc.p("Not confirmed for the variant (sample): " + ", ".join(
                f"{x.get('pdb_id')} ({_clip(x.get('title'), 60)})" for x in e["discarded_not_confirmed"][:5]),
                small=True)

    af = next((d["result"].get("alphafold") for d in items if (d["result"].get("alphafold") or {}).get("ok")), None)
    doc.h(3, "Predicted structure (AlphaFold DB)")
    if not af or not af.get("available", True):
        doc.p((af or {}).get("note") or NO_RECORD)
        return
    links = af.get("links") or {}
    frac = af.get("fraction_of_residues") or {}
    doc.p(Link(af.get("entry_id"), links.get("alphafold")),
          f": mean pLDDT {_num(af.get('mean_plddt')) or NO_RECORD} ({af.get('confidence_band') or NO_RECORD}), "
          f"residues {af.get('residues') or NO_RECORD}, model v{af.get('model_version')}. "
          f"Fraction of residues: very high (>90) {_num(frac.get('very_high_over_90'))}, "
          f"confident (70-90) {_num(frac.get('confident_70_90'))}, low (50-70) {_num(frac.get('low_50_70'))}, "
          f"very low (<50) {_num(frac.get('very_low_under_50'))}. ",
          af.get("selection_reason") or "", " Files: ",
          *[x for k in ("pdb", "cif", "pae_image") if links.get(k) for x in (Link(k.upper(), links[k]), " ")])
    doc.p("A single static prediction without ligands, nucleotides or ions.", small=True)


def _potency_lines(lig):
    lines = []
    for p in (lig.get("potency") or [])[:4]:
        cells = f" [{', '.join(p['cell_types'][:2])}]" if p.get("cell_types") else ""
        lines.append(f"{p.get('mutation') or 'not stated'} · {p.get('standard_type')} · {p.get('assay_class')}: "
                     f"{_num(p.get('best_nM'))} nM (median {_num(p.get('median_nM'))}, n={p.get('n_measurements')}){cells}")
    extra = len(lig.get("potency") or []) - 4
    if extra > 0:
        lines.append(f"+{extra} more groups in dossier")
    return Lines(lines) if lines else None


def _ligand_sources(lig, limit=3):
    seen, out = set(), []
    for p in lig.get("potency") or []:
        for doc in p.get("documents") or []:
            links = doc.get("links") or {}
            url = links.get("pubmed") or links.get("doi") or links.get("chembl_document")
            if url and url not in seen:
                seen.add(url)
                label = f"PMID {doc['pubmed_id']}" if doc.get("pubmed_id") else (doc.get("doi") or doc.get("document_chembl_id"))
                out.append(Link(f"{label} ({doc.get('year')})" if doc.get("year") else label, url))
    return Lines(out[:limit]) if out else None


def _ligands(doc, items):
    doc.h(2, "4. Ligands and bioactivity (ChEMBL)")
    for d in items:
        r = d["result"]
        mut = r.get("mutation_filter")
        doc.h(3, f"Ligands, variant '{mut}'" if mut else "Ligands, all variants")
        t = r.get("target")
        if not t:
            doc.p("No record in ChEMBL for this protein as a single-protein target.")
            continue
        doc.p("Target: ", Link(f"{t.get('chembl_id')} {t.get('name')}", t.get("link")),
              f" ({t.get('type')}, {t.get('organism')}). Retrieved {d.get('retrieved_at')}.", small=True)
        if r.get("other_targets_ignored"):
            doc.p("ChEMBL targets that include this protein but were not used: " + "; ".join(
                f"{x.get('chembl_id')} {x.get('name')} ({x.get('type')})" for x in r["other_targets_ignored"]), small=True)
        a = r.get("activities") or {}
        excluded = "; ".join(f"{k.replace('_', ' ')}: {v}" for k, v in (a.get("excluded") or {}).items())
        cut = a.get("ranked_download_cutoff_pchembl")
        unchecked = r.get("mechanism_compounds_unchecked") or []
        # Dossiers written before 2026-10-21 have no "mechanisms" block: not flagged.
        mech_failed = (r.get("mechanisms") or {}).get("retrieved") is False
        mech_text = {
            "mechanism_lookup_failed": ", mechanism lookup failed: no extra download for drugs",
            "failed": f", the extra download for mechanism compounds failed ({len(unchecked)} unchecked; "
                      "drug statistics only count values above the cutoff)",
            "no_mechanism_compounds": "",
            "partial": ", plus part of the activities of mechanism compounds (some batches cut; see Methods)",
            "complete": ", plus every activity of mechanism compounds",
        }.get(a.get("mechanism_compounds_download"))
        if mech_text is None:   # dossiers written before mechanism_compounds_download existed
            mech_text = (f", the extra download for mechanism compounds failed ({len(unchecked)} unchecked)"
                         if unchecked else ", plus the activities of mechanism compounds (see Methods)")
        trunc = (f" (truncated: most potent first, lowest pChEMBL fetched {_num(cut)}{mech_text})"
                 if cut is not None else " (truncated)") if a.get("truncated") else ""
        conf = (r.get("confidence") or {}).get("min_applied")
        doc.p(f"Activities: {a.get('total_in_chembl')} in ChEMBL, {a.get('fetched')} fetched"
              f"{trunc}, {a.get('used')} used. Excluded: {excluded or 'none'}. "
              f"Minimum assay confidence: {conf if conf is not None else 'not applied'}.", small=True)
        doc.p(r.get("selection") or "")
        if mech_failed:
            doc.p("Warning: the ChEMBL mechanism lookup failed for this file. Drugs without IC50/Ki/Kd "
                  "(e.g. antibodies) and mechanisms of action were not checked: their absence below is not "
                  "evidence that there are none.")
        rows = []
        for lig in r.get("ligands") or []:
            rows.append([
                _compound_cell(lig),
                _phase_text(lig),
                (lig.get("mechanism") or {}).get("mechanism_of_action"),
                _num(lig.get("best_pchembl")),
                _potency_lines(lig),
                _ligand_sources(lig),
            ])
        doc.table(["Compound", "Max phase", "Mechanism", "pChEMBL (best)", "Potency (best per group)", "Sources"],
                  rows, (0.15, 0.10, 0.13, 0.09, 0.35, 0.18))
        items = r.get("mechanism_compounds_without_qualifying_activity") or []
        # Dossiers written before 2026-10-21 hold bare ChEMBL ids here.
        if items and all(isinstance(x, str) for x in items):
            doc.p("Compounds with a recorded mechanism but no qualifying activity: " + ", ".join(items), small=True)
        else:
            _mechanism_table(doc, items, "Drugs and clinical candidates with a recorded mechanism but no qualifying "
                             "IC50/Ki/Kd ({n}; typically antibodies or other proteins, so they are not ranked by potency):")
        _mechanism_table(doc, r.get("mechanism_compounds_other_variants_only") or [],
                         f"Compounds with a recorded mechanism and qualifying IC50/Ki/Kd only for variants other "
                         f"than '{mut}' ({{n}}; no qualifying measurement on this variant in ChEMBL, which is not the same "
                         f"as untested):")
        _mechanism_table(doc, unchecked,
                         "Compounds with a recorded mechanism that could not be checked because their activity "
                         "download failed ({n}; they may have potency data):")


def _mechanism_table(doc, items, intro):
    if not items:
        return
    doc.p(intro.format(n=len(items)))
    rows = [[_compound_cell(x), x.get("molecule_type"), _phase_text(x), x.get("mechanism_of_action")]
            for x in items]
    doc.table(["Compound", "Type", "Max phase", "Mechanism"], rows, (0.25, 0.15, 0.20, 0.40))


def _literature(doc, items):
    doc.h(2, "5. Literature (Europe PMC)")
    for d in items:
        r = d["result"]
        params = d.get("params") or {}
        topic = params.get("topic")
        doc.h(3, "Literature" + (f", topic '{topic}'" if topic else "") + f", sorted by {r.get('sort') or params.get('sort')}")
        doc.p(f"{r.get('returned')} of {r.get('total_hits')} matching records. Query: {r.get('query_used')}. "
              f"Retrieved {d.get('retrieved_at')}.", small=True)
        rows = []
        for i, p in enumerate(r.get("papers") or [], 1):
            types = ", ".join(p.get("pub_types") or []) + (" · PREPRINT" if p.get("is_preprint") else "")
            rows.append([str(i), Lines([Link(p.get("title"), _paper_link(p)), _clip(p.get("authors"), 90) or ""]),
                         f"{p.get('journal') or NO_RECORD} ({p.get('year') or 'n.d.'})", types or None,
                         _num(p.get("cited_by_count"))])
        doc.table(["#", "Title / authors", "Journal (year)", "Type", "Cited by"], rows, (0.04, 0.54, 0.20, 0.14, 0.08))


def load_manifest(accession: str):
    """prep/manifest.json written by prepare_structures (Phase 2), {"unreadable": error} or None."""
    path = choice.manifest_path(accession)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as exc:
        log.warning("Unreadable preparation manifest %s", path)
        return {"unreadable": str(exc)}
    return data if isinstance(data, dict) and data.get("generated_by") == "prepare_structures" else None


def _ranges(ranges):
    return ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in ranges)


def _chain_line(c):
    """One target chain: modelled UniProt range, internal gaps, residue at the variant position."""
    rng = c.get("modelled_uniprot_range")
    name = c["auth_asym_id"] + (f" (label {c['label_asym_id']})"
                                if c.get("label_asym_id") not in (None, c["auth_asym_id"]) else "")
    line = f"{name}: modelled {rng[0]}-{rng[1]}" if rng else f"{name}: modelled range {NO_RECORD}"
    missing = c.get("missing_uniprot_ranges") or []
    # without a UniProt alignment there is no modelled range, so no gap statement either
    line += f"; missing {_ranges(missing)}" if missing else "; no internal gap" if rng else "; gaps: no record"
    # manifest chains carry variant_check; a recorded choice carries variant_status/variant_observed
    vc = c.get("variant_check") or ({"status": c["variant_status"], "observed": c.get("variant_observed")}
                                     if c.get("variant_status") else None)
    if vc:
        line += f"; variant position: {vc.get('status')}" + (f" ({vc['observed']})" if vc.get("observed") else "")
    return line


def _box(box):
    return (" × ".join(_num(x) or "?" for x in box.get("size") or []) + " Å, centre ("
            + ", ".join(_num(x) or "?" for x in box.get("center") or []) + ")")


def _pocket_row(pdb_id, pk, lig):
    ligand = Lines([f"{pk['comp_id']} {pk['label_asym_id']} (auth {pk.get('auth_asym_id')})"]
                   + ([", ".join(lig["synonyms"])] if lig.get("synonyms") else
                      [_clip(lig["name"], 60)] if lig.get("name") else []))
    if pk.get("error"):
        return [pdb_id, ligand, f"error: {pk['error']}", None, None, None]
    n = len(pk.get("residues") or [])
    flags = ((["spans several target chains"] if pk.get("spans_target_chains") else [])
             + (["large box (side > 30 Å)"] if pk.get("large_box") else [])
             + [f"missing {_ranges([a['missing_uniprot_range']])} next to the pocket (chain {a['auth_asym_id']})"
                for a in pk.get("adjacent_missing_ranges") or []])
    af = pk.get("alphafold_plddt") or {}
    return [pdb_id, ligand,
            Lines([f"{n} residues"] + [f"chain {c}: {k}" for c, k in (pk.get("residues_by_chain") or {}).items()]),
            _box(pk.get("box") or {}),
            Lines(flags) if flags else "none",
            # residues_compared counts distinct UniProt positions; n counts residues over all chains (5JSM: 24 vs 47)
            f"{_num(af['mean_plddt'])} (at {af.get('residues_compared')} UniProt positions)"
            if af.get("mean_plddt") is not None else NO_RECORD]


def _interfaces(doc, m):
    """Step 5 (analyze_interfaces): partners, assemblies, PISA and contacts, as stored in the manifest."""
    done = [s for s in m.get("structures") or [] if s.get("interfaces")]
    doc.h(3, "Protein-protein interfaces")
    if not done:
        doc.p("Not analysed (analyze_interfaces not run on this manifest).")
        return
    doc.p(f"{len(done)} of {len(m.get('structures') or [])} structures analysed. Contacts: target residues with a heavy "
          "atom within 5 Å of another polymer chain, including crystal-symmetry copies (code 1_555 = deposited "
          "coordinates). Geometry only: no contact is classed as biological here. RCSB's assembly and PDBe PISA's "
          "opinion are quoted as given; they can disagree.", small=True)
    for s in done:
        b = s["interfaces"]
        doc.p(f"{s['pdb_id']} (analysed {b.get('computed_at')}; space group {b.get('spacegroup') or NO_RECORD})", small=False)
        items = [f"{q['role']} entity {q['entity_id']}: {q['description']}; chains {', '.join(c['auth'] for c in q['chains'])}; "
                 f"UniProt {', '.join(q['uniprot_ids'] or []) or 'none'}; {q['length'] or NO_RECORD} residues"
                 + "".join(f"; {f}" for f in q["flags"][1:]) for q in b.get("sequences") or []]
        if b.get("rcsb_assemblies") is None:
            items.append(f"RCSB assembly: {NO_RECORD}")
        for a in b.get("rcsb_assemblies") or []:
            items.append(f"RCSB assembly {a['id']}: {a['defined_by']}" + (f" ({a['method']})" if a.get("method") else "")
                         + f", polymer chains {', '.join(a['polymer_chains'])}, operators {', '.join(a['operators'])}"
                         + ("; same as the asymmetric unit" if a.get("equals_asymmetric_unit") else ""))
        p = b.get("pisa") or {}
        if p.get("multimer_sets") is not None:
            items += [f"PISA multimer set {x['set']} (one alternative; its assemblies coexist): "
                      + " + ".join(f"{a['composition']}: {a['verdict']}" for a in x["assemblies"])
                      for x in p["multimer_sets"]]
        for key in ("interfaces_error", "multimers_error"):
            if p.get(key):
                items.append(f"PISA {key.split('_')[0]}: not available ({p[key]})")
        doc.bullets(items)
        rows = []
        for c in b.get("contacts") or []:
            pi = c.get("pisa_interface")
            nums = ", ".join(str(r["uniprot"]) if r["uniprot"] is not None else f"{r['auth_seq_id']}(auth)"
                             for r in c["target_residues"])
            rows.append([f"{c['target_chain']}", f"{c['partner_chain']} ({c['partner_role']})",
                         Lines([c["symmetry"]["code"], c["symmetry"]["operator"] or NO_RECORD,
                                "asymmetric unit" if c["in_asymmetric_unit"] else "symmetry copy"]),
                         f"{c['n_target_residues']}: {nums}",
                         f"#{pi['id']}, css {_num(pi['css']) or NO_RECORD}, "
                         f"{_num(round(pi['area'])) if pi.get('area') is not None else NO_RECORD} Å², {pi['target_residues_in_pisa']}"
                         if pi else "none at this chain pair and operator" if c.get("pisa_checked") else "PISA not available"])
        if b.get("contacts_error"):
            doc.p(f"Contacts not computed: {b['contacts_error']}", small=True)
        doc.table(["Target chain", "Partner chain", "Symmetry", "Target residues (UniProt)", "PISA interface"],
                  rows, (0.09, 0.14, 0.16, 0.43, 0.18))


def _preparation(doc, m, entries):
    doc.h(2, "6. Structure preparation (Phase 2)")
    if m.get("unreadable"):
        doc.p(f"prep/manifest.json could not be read ({m['unreadable']}).")
        return
    src = m.get("source_dossier") or {}
    doc.p(f"Manifest prep/manifest.json created {m.get('created_at')} by prepare_structures from {src.get('file')} "
          f"(retrieved {src.get('retrieved_at')}). Variant filter: {m.get('mutation_filter') or 'none'}. "
          f"{m.get('selection') or ''}", small=True)
    now = next((d for d in entries.get("get_structures") or [] if d["file"] == src.get("file")), None)
    if not now:
        doc.p(f"Note: {src.get('file')} is not in the dossier any more; this section shows the manifest only.")
    elif now.get("retrieved_at") != src.get("retrieved_at"):
        doc.p(f"Note: {src.get('file')} in the dossier was retrieved {now.get('retrieved_at')}; the manifest was built "
              f"from the copy retrieved {src.get('retrieved_at')}. Sections 3 and 6 may show different selections.")

    rows, pocket_rows, n_instances, n_pockets, n_prepared = [], [], 0, 0, 0
    for s in m.get("structures") or []:
        prep = s.get("prepared") or {}
        instances = s.get("nonpolymer_instances") or []
        # only category "ligand" gets a pocket (clean.py); nucleotides/cofactors stay in the receptor
        if prep.get("ok"):
            n_prepared += 1
            n_instances += sum(1 for i in instances if i.get("category") == "ligand")
            receptor = prep["receptor"]["path"]
            chains = Lines(_chain_line(c) for c in prep["metrics"]["chains"])
            r_free = _num(prep["metrics"].get("r_free")) or "not reported"
        else:
            receptor = "not prepared: " + str(prep.get("error") or (s.get("file") or {}).get("error") or NO_RECORD)
            chains = r_free = "not checked"
        res = s.get("resolution_angstrom")
        rows.append([Link(s.get("pdb_id"), (s.get("links") or {}).get("rcsb")),
                     f"{_num(res)} Å" if res is not None else NO_RECORD, r_free, receptor, chains,
                     _mutation_cell(s.get("construct") or {})])
        by_label = {i.get("label_asym_id"): i for i in instances}
        for pk in prep.get("pockets") or []:
            pocket_rows.append(_pocket_row(s["pdb_id"], pk, by_label.get(pk.get("label_asym_id")) or {}))
            n_pockets += not pk.get("error")
    doc.h(3, "Structures in the manifest")
    doc.table(["PDB", "Resolution", "R-free", "Receptor file", "Target chains (UniProt numbering)", "Mutations"],
              rows, (0.07, 0.08, 0.07, 0.17, 0.41, 0.20))
    doc.h(3, "Ligand pockets")
    doc.p(f"{n_pockets} pockets computed of {n_instances} ligand instances in the {n_prepared} of "
          f"{len(m.get('structures') or [])} entries that were prepared (nucleotides, cofactors and additives are "
          "not counted as ligands; unprepared entries have no pockets). Pocket = target residues within 5 Å of the "
          "ligand instance; box = ligand + 5 Å per side, mmCIF frame.", small=True)
    doc.table(["PDB", "Ligand instance", "Pocket residues", "Box (size, centre)", "Flags", "AlphaFold mean pLDDT"],
              pocket_rows, (0.07, 0.22, 0.14, 0.25, 0.18, 0.14))
    without = (m.get("summary") or {}).get("ligands_without_pocket") or []
    if without:
        doc.p("Ligand instances without a pocket:", small=True)
        doc.bullets(without)
    af = m.get("alphafold")
    if af:
        f, vp = af.get("file") or {}, af.get("variant_position") or {}
        doc.p(Link(af.get("entry_id"), f.get("url")),
              f": file {f.get('path')}" if f.get("ok") else f": not downloaded ({f.get('error')})",
              f"; residue at UniProt {vp['uniprot_position']}: {vp.get('residue')} (pLDDT {_num(vp.get('plddt'))})"
              if vp.get("uniprot_position") else "", small=True)

    _interfaces(doc, m)

    doc.h(3, "Structure chosen by the user")
    uc, history = m.get("user_choice"), m.get("user_choice_history") or []
    if not uc:
        doc.p("No structure chosen (select_structure not called" + (", or the last choice was dropped" if history else "") + ").")
    else:
        pk = uc.get("pocket") or {}
        items = [f"Choice: {choice.label(uc)}" + (f" ({', '.join(pk['synonyms'])})" if pk.get("synonyms") else ""),
                 f"Recorded {uc.get('chosen_at')} by {uc.get('recorded_by')}"
                 + (f"; revalidated against the manifest of {uc['revalidated_at']}" if uc.get("revalidated_at") else ""),
                 f'Reason recorded by select_structure: "{uc.get("reason")}"']
        if uc.get("kind") == "experimental":
            mut = _mutation_cell(uc.get("construct") or {})
            items += ([f"Receptor: {uc.get('receptor')}"] + [f"Chain {_chain_line(c)}" for c in uc.get("chains") or []]
                      + ["Construct mutations: " + ("; ".join(mut) if isinstance(mut, Lines) else mut)])
        else:
            items.append(f"AlphaFold file: {uc.get('file')}")
        if pk:
            items += [f"Pocket: {pk.get('n_residues')} residues; box {_box(pk.get('box') or {})}",
                      f"Ligand pose (SDF): {pk.get('sdf') or 'not available'}"]
        if uc.get("changed_fields"):
            items.append("Changed when the manifest was rebuilt: " + ", ".join(uc["changed_fields"]))
        doc.bullets(items + [f"Warning: {w}" for w in uc.get("warnings") or []])
    if history:
        doc.p(f"{len(history)} earlier choice(s) in user_choice_history: " + "; ".join(
            choice.label(h) + (f" (dropped: {h['dropped_reason']})" if h.get("dropped_reason") else " (replaced)")
            for h in history), small=True)
    problems = list(m.get("failed_downloads") or []) + list(m.get("warnings") or [])
    if problems:
        doc.p("Failed downloads and warnings recorded in the manifest:", small=True)
        doc.bullets(problems)


def _methods(doc, entries, manifest=None):
    doc.h(2, "7. Methods and limitations")
    doc.p("This document is generated deterministically from the dossier files listed at the top. It contains "
          "no interpretation: values are quoted from the source records, every link is built from an identifier "
          f"returned by the source, and fields without data are shown as '{NO_RECORD}'. Selections are partial "
          "by design; each section states how many records were shown out of how many were found.")
    notes = []
    for tool, _ in SECTIONS:
        for d in entries.get(tool) or []:
            r = d["result"]
            for note in list(r.get("notes") or []) + list((r.get("experimental") or {}).get("notes") or []):
                if note not in notes:
                    notes.append(note)
    if notes:
        doc.p("Selection rules and caveats stated by the tools that produced the data:")
        doc.bullets(notes)
    if manifest and not manifest.get("unreadable"):
        doc.p("Structure preparation rules (prep/manifest.json, decisions):")
        doc.bullets(list(manifest.get("decisions") or []))
        doc.p("Not done yet (pending): " + ("; ".join(manifest.get("pending") or []) or "nothing listed") + ".")


# ------------------------------------------------------------ rendering
def render_markdown(doc: Doc) -> str:
    out = []
    for b in doc.blocks:
        if b[0] == "h":
            out.append("#" * b[1] + " " + b[2])
        elif b[0] == "p":
            text = _md(b[1])
            out.append(f"<small>{text}</small>" if b[2] else text)
        elif b[0] == "ul":
            out.append("\n".join("- " + _md(i) for i in b[1]))
        elif b[0] == "table":
            _, headers, rows, _w = b
            lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
            lines += ["| " + " | ".join(_md(c) for c in row) + " |" for row in rows]
            out.append("\n".join(lines))
    return "\n\n".join(out) + "\n"


_FONT = None


def _font():
    """Arial from Windows (has α, ≥, →, ′...). The built-in PDF fonts lack them."""
    global _FONT, _CMAP
    if _FONT:
        return _FONT
    folder = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    files = {"Report": "arial.ttf", "Report-Bold": "arialbd.ttf",
             "Report-Italic": "ariali.ttf", "Report-BoldItalic": "arialbi.ttf"}
    if all((folder / f).exists() for f in files.values()):
        for name, f in files.items():
            pdfmetrics.registerFont(TTFont(name, str(folder / f)))
        pdfmetrics.registerFontFamily("Report", normal="Report", bold="Report-Bold",
                                      italic="Report-Italic", boldItalic="Report-BoldItalic")
        _CMAP = set(pdfmetrics.getFont("Report").face.charToGlyph)
        _FONT = ("Report", "Report-Bold", None)
    else:
        log.warning("Arial not found in %s: using Helvetica; Greek letters and symbols may not render.", folder)
        _FONT = ("Helvetica", "Helvetica-Bold", "Arial not found: symbols such as α or ≥ may not render.")
    return _FONT


def render_pdf(doc: Doc, path: Path, title: str, footer: str):
    regular, bold, _ = _font()
    st = {
        "body": ParagraphStyle("body", fontName=regular, fontSize=9, leading=12, spaceAfter=4),
        "small": ParagraphStyle("small", fontName=regular, fontSize=7.5, leading=9.5, spaceAfter=4,
                                textColor=colors.HexColor("#444444")),
        "cell": ParagraphStyle("cell", fontName=regular, fontSize=7.2, leading=8.8),
        "head": ParagraphStyle("head", fontName=bold, fontSize=7.4, leading=9),
        1: ParagraphStyle("h1", fontName=bold, fontSize=16, leading=20, spaceAfter=6),
        2: ParagraphStyle("h2", fontName=bold, fontSize=12.5, leading=16, spaceBefore=10, spaceAfter=4,
                          textColor=colors.HexColor("#123b6d")),
        3: ParagraphStyle("h3", fontName=bold, fontSize=10, leading=13, spaceBefore=6, spaceAfter=3),
    }
    width = A4[0] - 32 * mm
    flow = []
    for b in doc.blocks:
        if b[0] == "h":
            flow.append(Paragraph(escape(_printable(b[2])), st[b[1]]))
        elif b[0] == "p":
            flow.append(Paragraph(_pdf(b[1]), st["small" if b[2] else "body"]))
        elif b[0] == "ul":
            flow.append(ListFlowable([ListItem(Paragraph(_pdf(i), st["body"]), leftIndent=10) for i in b[1]],
                                     bulletType="bullet", start="•", leftIndent=10, bulletFontSize=7))
        elif b[0] == "table":
            _, headers, rows, widths = b
            data = [[Paragraph(escape(_printable(h)), st["head"]) for h in headers]]
            data += [[Paragraph(_pdf(c), st["cell"]) for c in row] for row in rows]
            t = Table(data, colWidths=[w * width for w in widths], repeatRows=1)
            t.setStyle(TableStyle([
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#b0b7c3")),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8edf4")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ]))
            flow += [t, Spacer(1, 5)]

    def on_page(canvas, d):
        canvas.saveState()
        canvas.setFont(regular, 7)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(16 * mm, 9 * mm, footer)
        canvas.drawRightString(A4[0] - 16 * mm, 9 * mm, f"Page {d.page}")
        canvas.restoreState()

    template = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=16 * mm, rightMargin=16 * mm,
                                 topMargin=15 * mm, bottomMargin=15 * mm, title=title,
                                 author="protein-research MCP")
    template.build(flow, onFirstPage=on_page, onLaterPages=on_page)
    return template.page


# ------------------------------------------------------------ public API
def build_document(accession: str, entries: dict, generated: str, manifest=None) -> Doc:
    prof = _first(entries, "get_protein_profile")
    lit = _first(entries, "search_literature")
    up = prof["result"]["uniprot"] if prof else (lit["result"].get("protein") if lit else {}) or {}
    gene, name = up.get("gene"), up.get("protein_name")

    doc = Doc()
    doc.h(1, f"Target dossier: {gene or accession}" + (f", {name}" if name else ""))
    doc.p(Link(f"UniProtKB {accession}", up.get("source_url")),
          f" · organism: {up.get('organism') or NO_RECORD} · generated {generated}", small=True)
    doc.p("Data report built only from source records saved in the dossier. Interpretation is not part of this "
          "document.", small=True)
    _coverage(doc, accession, entries, manifest)
    if prof:
        _uniprot(doc, prof)
        _open_targets(doc, prof)
    if entries.get("get_structures"):
        _structures(doc, entries["get_structures"])
    if entries.get("get_ligands"):
        _ligands(doc, entries["get_ligands"])
    if entries.get("search_literature"):
        _literature(doc, entries["search_literature"])
    if manifest:
        _preparation(doc, manifest, entries)
    _methods(doc, entries, manifest)
    return doc


def build_report(uniprot_id: str) -> dict:
    accession, err = parse_accession(uniprot_id)
    if err:
        return err
    entries = load_dossier(accession)
    manifest = load_manifest(accession)
    if not entries:
        return {"ok": False, "uniprot_id": accession,
                "error": f"No dossier for {accession}. Run get_protein_profile, get_structures, get_ligands "
                         "and/or search_literature first: each one saves its result to the dossier."}
    generated = utc_now()
    doc = build_document(accession, entries, generated, manifest)
    folder = dossier.DOSSIER_DIR / accession
    md_path, pdf_path = folder / f"report_{accession}.md", folder / f"report_{accession}.pdf"
    gene = next((d["result"]["uniprot"].get("gene") for d in entries.get("get_protein_profile") or []), None)
    try:
        md_path.write_text(render_markdown(doc), encoding="utf-8")
        pages = render_pdf(doc, pdf_path, title=f"Target dossier {gene or ''} {accession}".strip(),
                           footer=f"{gene or accession} ({accession}) · generated {generated} from source records only")
    except OSError as exc:  # typically the PDF is open in a viewer on Windows
        log.exception("Could not write the report")
        return {"ok": False, "uniprot_id": accession,
                "error": f"Could not write the report ({exc}). If the PDF is open, close it and try again."}
    missing = [tool for tool, _ in SECTIONS if not entries.get(tool)]
    font_note = _font()[2]
    return {
        "ok": True,
        "uniprot_id": accession,
        "pdf_path": str(pdf_path),
        "markdown_path": str(md_path),
        "pages": pages,
        "generated": generated,
        "included": {tool: [d["file"] for d in entries[tool]] for tool, _ in SECTIONS if entries.get(tool)},
        "missing": missing,
        "preparation": ("prep/manifest.json" if manifest and not manifest.get("unreadable")
                        else "prep/manifest.json unreadable" if manifest else None),
        **({"font_warning": font_note} if font_note else {}),
    }
