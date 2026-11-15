"""Self-contained HTML report (figures embedded as base64, plotly.js inlined)."""

import base64
import html
import re
from datetime import date
from pathlib import Path

import pandas as pd

CATALOG = "https://www.pgscatalog.org"

SHORTLIST_COLS = {
    "pgs_id": "PGS", "name": "Name", "evidence": "Level", "n_variants": "Variants",
    "dev_publication": "Publication", "gwas_ancestry": "GWAS ancestry",
    "n_evals": "Evals", "n_external_pubs": "External pubs", "eval_ancestry": "Eval. ancestries",
    "max_eval_cases": "Max cases", "max_eval_n": "Max N", "auc_pgs_only": "AUC PGS only",
    "or_median": "Median OR", "hr_median": "Median HR", "r2_pgs_only": "R² PGS only",
    "notes": "Caveats", "other_phenotypes": "Other phenotypes evaluated",
    "source_traits": "Found via",
}

# Light theme only, whatever the OS / browser setting
CSS = """
:root{color-scheme:light;--bg:#fcfcfb;--card:#ffffff;--ink:#0b0b0b;--ink2:#52514e;
--line:#e6e5e1;--accent:#2a78d6;--A:#0d3a75;--B:#2a78d6;--C:#86b6ec;--D:#bcd4f2}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1280px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:21px;margin:48px 0 8px;padding-top:16px;
border-top:1px solid var(--line)}h3{font-size:16px;margin:28px 0 8px}
.muted{color:var(--ink2)}a{color:var(--accent)}
.kpis{display:flex;flex-wrap:wrap;gap:12px;margin:16px 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 16px;min-width:150px}
.kpi b{display:block;font-size:24px}.kpi span{color:var(--ink2);font-size:12px}
.tbl{overflow-x:auto;border:1px solid var(--line);border-radius:8px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{position:sticky;top:0;background:var(--card);color:var(--ink2);font-weight:600;white-space:nowrap}
td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.lvl{display:inline-block;width:22px;text-align:center;border-radius:4px;font-weight:700;color:#fff}
.lvl-A{background:var(--A)}.lvl-B{background:var(--B)}.lvl-C{background:var(--C);color:#0b0b0b}
.lvl-D{background:var(--D);color:#0b0b0b}
img{max-width:100%;border:1px solid var(--line);border-radius:8px;margin:8px 0}
.plot{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px}
details{margin:12px 0}summary{cursor:pointer;font-weight:600}
.box{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);
border-radius:8px;padding:12px 16px;margin:16px 0}
"""

GUIDE = """
<div class="box"><b>How to read this report</b><ul>
<li><b>Evidence level</b>: A = validated by an independent publication on ≥1,000 cases
(or N ≥ 10,000 for a quantitative trait); B = either a large sample or an external validation;
C = only small evaluations by the authors; D = never evaluated on the trait itself.</li>
<li><b>Metrics are only comparable within the same study</b> (same sample, phenotype and
covariates): see the head-to-head comparisons. An AUC "with covariates" includes age, sex…
and can be very high even if the PGS adds little.</li>
<li>OR / HR: check the unit in the publication (per standard deviation, top vs. bottom decile…).</li>
<li>An "external" evaluation comes from another publication than the one that developed the
score (it may still involve the same consortia). Evaluations on other phenotypes
(e.g. a body-weight PGS tested on hearing) are excluded from the metrics.</li>
<li>Before use, check the ancestry of your target population, the genome build and variant
coverage in your data, and the license.</li></ul></div>
"""


def _img(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{data}" alt="{html.escape(path.stem)}">'


def _cell(col, v):
    if col == "pgs_id":
        return f'<td><a href="{CATALOG}/score/{v}" target="_blank">{v}</a></td>'
    if col == "evidence":
        return f'<td><span class="lvl lvl-{v}">{v}</span></td>'
    if isinstance(v, bool):
        return f"<td>{'yes' if v else 'no'}</td>"
    if isinstance(v, (int, float)):
        if pd.isna(v):
            return '<td class="num muted">–</td>'
        txt = f"{v:,.0f}" if float(v).is_integer() or abs(v) >= 100 else f"{v:.3g}"
        return f'<td class="num">{txt}</td>'
    return f"<td>{html.escape(str(v)) if pd.notna(v) else ''}</td>"


def table(df: pd.DataFrame, cols: dict | None = None) -> str:
    cols = cols or {c: c for c in df.columns}
    head = "".join(f"<th>{html.escape(lbl)}</th>" for lbl in cols.values())
    body = "".join(
        "<tr>" + "".join(_cell(c, r[c]) for c in cols) + "</tr>" for _, r in df.iterrows()
    )
    return f'<div class="tbl"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _head_to_head_html(h2h: pd.DataFrame) -> str:
    if h2h.empty:
        return '<p class="muted">No study evaluates several of these scores side by side.</p>'
    parts = []
    for (author, pheno, anc, cohorts, metric), g in h2h.groupby(
        ["eval_author", "phenotype", "ancestry", "cohorts", "metric"], sort=False
    ):
        n = int(g["n"].max())
        cases = int(g["n_cases"].max())
        extra = f", {cases:,} cases" if cases else ""
        where = f"{cohorts}, " if cohorts else ""
        parts.append(
            f"<details><summary>{html.escape(author)} — {html.escape(pheno)} — "
            f"{html.escape(where + anc)} (N={n:,}{extra}) — {html.escape(metric)}</summary>"
            + table(g, {"rank": "Rank", "pgs_id": "PGS", "name": "Name", "estimate": "Estimate",
                        "ci_lower": "CI lower", "ci_upper": "CI upper"})
            + "</details>"
        )
    return "".join(parts)


OVERVIEW_JS = f"""<script>
document.getElementById("all-scores").on("plotly_click", e => {{
  const id = e.points[0].customdata[0];
  window.open("{CATALOG}/score/" + id, "_blank");
}});
</script>"""


def _anchor(key: str) -> str:
    return re.sub(r"\W+", "_", key)


def _members_html(trait, groups) -> str:
    groups = groups or {}
    members = groups.get(trait, {}).get("members") or [(trait, trait)]
    links = []
    for tid, label in members:
        page = "score" if tid.startswith("PGS") else "trait"
        name = tid if label in (tid, groups.get(trait, {}).get("label")) else f"{label} {tid}"
        links.append(f'<a href="{CATALOG}/{page}/{tid}" target="_blank">{html.escape(name)}</a>')
    return '<span class="muted" style="font-size:14px;font-weight:400"> ' + " · ".join(links) + "</span>"


def _filters(cfg) -> str:
    lo, hi, n_snp = cfg.get("date_min"), cfg.get("date_max"), cfg.get("n_snp")
    parts = []
    if lo is not None or hi is not None:
        parts.append(f"published {lo or '…'}–{hi or '…'}")
    if n_snp is not None:
        parts.append(f"at least {n_snp:,} variants")
    return ", ".join(parts) or "none"


def build(summary, h2h, figures, cfg, out_file: Path, overview_html: str = "", groups=None):
    """`figures`: {trait group: [(title, png_path), ...]}; `overview_html`: interactive figure
    of all scores shown first; `groups`: {trait group: {"label", "members": [(id, label)]}}."""
    traits = summary["query_trait"].unique()
    loci_cols = {c: c.replace("_", " ") for c in summary.columns
                 if c.endswith(("_snps", "_region_variants"))}
    target = cfg.target_ancestry
    sections = []
    for trait in traits:
        s = summary[summary["query_trait"] == trait]
        label = s["trait_label"].iloc[0]
        kpis = [
            (len(s), "scores"),
            (int(s["n_evals"].sum()), "evaluations of the trait"),
            (int((s["n_external_pubs"] > 0).sum()), "externally validated scores"),
            (int(s["target_evaluated"].sum()), f"scores evaluated in {target}"),
            (int((s["evidence"] == "A").sum()), "level A scores"),
        ]
        shortlist = s[s["evidence"].isin(["A", "B"])]
        cols = {**SHORTLIST_COLS, **loci_cols}
        cols = {k: v for k, v in cols.items() if k in s.columns and s[k].notna().any()}
        figs = "".join(f"<h3>{html.escape(t)}</h3>{_img(p)}" for t, p in figures.get(trait, []))
        sections.append(
            f'<h2 id="{_anchor(trait)}">{html.escape(label)} {_members_html(trait, groups)}</h2>'
            + '<div class="kpis">' + "".join(f'<div class="kpi"><b>{v}</b><span>{k}</span></div>'
                                             for v, k in kpis) + "</div>"
            + f"<h3>Shortlist (levels A and B) — {len(shortlist)} scores</h3>"
            + (table(shortlist, cols) if len(shortlist) else '<p class="muted">None.</p>')
            + "<h3>Head-to-head comparisons (same study, same sample)</h3>"
            + _head_to_head_html(h2h[h2h["query_trait"] == trait])
            + figs
            + f"<details><summary>All scores ({len(s)})</summary>{table(s, cols)}</details>"
        )
    nav = " · ".join(
        f'<a href="#{_anchor(t)}">'
        f'{html.escape(summary.loc[summary.query_trait == t, "trait_label"].iloc[0])}</a>'
        for t in traits
    )
    overview = (
        f"<h2 style='border:0;margin-top:24px'>All scores ({len(summary)})</h2>"
        "<p class='muted'>Hover over a point to see the score, click to open its PGS Catalog "
        "page. Evidence levels: see the box below.</p>"
        f"<div class='plot'>{overview_html}</div>{OVERVIEW_JS}"
        if overview_html else ""
    )
    empty = [info["label"] for g, info in (groups or {}).items() if g not in set(traits)]
    emptied = (" — no score left for: " + html.escape(", ".join(empty))) if empty else ""
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light"><title>PGScaler report</title>
<style>{CSS}</style></head><body><main>
<h1>PGScaler report</h1>
<p class="muted">Polygenic scores from the PGS Catalog — generated on {date.today():%Y-%m-%d} —
target population: {html.escape(target)} — traits: {nav}</p>
<p class="muted">Filters: {html.escape(_filters(cfg))}{emptied}</p>
{overview}
{GUIDE}{''.join(sections)}
</main></body></html>"""
    out_file.write_text(page, encoding="utf-8")
