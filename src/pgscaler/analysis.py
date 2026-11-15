"""Turn raw PGS Catalog tables into a per-score decision table."""

import re

import pandas as pd

ANCESTRY_CODES = {
    "European": "EUR",
    "East Asian": "EAS",
    "South Asian": "SAS",
    "African unspecified": "AFR",
    "African American or Afro-Caribbean": "AFR",
    "Sub-Saharan African": "AFR",
    "Hispanic or Latin American": "AMR",
    "Greater Middle Eastern (Middle Eastern, North African or Persian)": "GME",
    "Asian unspecified": "ASN",
    "Not reported": "NR",
    "NR": "NR",
}

# Metrics reported for the PGS alone vs. a model also including covariates
DISCRIMINATION = ["AUROC", "C-index", "PGS AUROC (no covariates)"]
VARIANCE = ["R²", "PGS R2 (no covariates)"]
EFFECTS = ["OR", "HR"]
HEAD_TO_HEAD_METRICS = DISCRIMINATION + VARIANCE + EFFECTS + ["β"]

# Words too generic to say that an evaluated phenotype is the queried trait
STOPWORDS = {
    "disease", "diseases", "disorder", "trait", "measurement", "level", "levels", "self",
    "reported", "risk", "onset", "time", "event", "body", "with", "and", "from", "type",
}


def ancestry_codes(text: str) -> set[str]:
    out = set()
    for part in re.split(r"\s*\|\s*|,\s*", str(text)):
        if part:
            out.add(ANCESTRY_CODES.get(part, "OTH"))
    return out


def _tokens(text: str) -> set[str]:
    """Word stems (5 first letters) so that e.g. depressive ~ depression, autism ~ autistic."""
    words = re.findall(r"[a-z]+", str(text).lower().replace("'s", ""))
    return {w[:5] for w in words if len(w) >= 4 and w not in STOPWORDS}


def _has_acronym(label: str, phenotype: str) -> bool:
    """e.g. 'major depressive disorder' -> MDD found in 'Incident MDD'."""
    acronym = "".join(w[0] for w in re.findall(r"[a-z]+", label.lower()))
    return len(acronym) >= 3 and re.search(rf"\b{acronym}\b", phenotype, re.IGNORECASE) is not None


def annotate_performance(perf, scores, aliases, patterns) -> pd.DataFrame:
    """Add query_trait, `primary` (evaluated phenotype = queried trait) and
    `external` (evaluated by another publication than the one developing the score).

    `aliases`: {trait group: [labels of the group, its traits and their children]}.
    `patterns`: {trait group: regex} overriding automatic matching."""
    s = scores[["id", "query_trait", "trait_reported", "publication.id"]].rename(
        columns={"id": "pgs_id", "publication.id": "dev_pub_id"}
    )
    p = perf.merge(s, on="pgs_id")
    p["external"] = p["eval_pub_id"] != p["dev_pub_id"]

    def is_primary(row):
        pattern = patterns.get(row.query_trait)
        if pattern:
            return re.search(pattern, row.phenotype, flags=re.IGNORECASE) is not None
        names = aliases.get(row.query_trait, [])
        ref = set().union(*map(_tokens, names)) | _tokens(row.trait_reported)
        return bool(ref & _tokens(row.phenotype)) or any(
            _has_acronym(n, row.phenotype) for n in names)

    flags = p.drop_duplicates(["ppm_id", "query_trait"]).copy()
    flags["primary"] = flags.apply(is_primary, axis=1)
    p = p.merge(flags[["ppm_id", "query_trait", "primary"]], on=["ppm_id", "query_trait"])
    # An AUROC/R² without covariates is a PGS-only metric
    no_cov = p["covariates"].isna() | (p["covariates"].astype(str).str.strip() == "")
    p.loc[no_cov & (p["metric"] == "AUROC"), "metric"] = "PGS AUROC (no covariates)"
    p.loc[no_cov & (p["metric"] == "R²"), "metric"] = "PGS R2 (no covariates)"
    return p


def _fmt_dist(dist) -> str:
    if not isinstance(dist, dict) or not dist:
        return ""
    return ", ".join(f"{k} {v:g}%" for k, v in sorted(dist.items(), key=lambda kv: -kv[1]))


def evidence_level(ext_pubs: int, max_cases: float, max_n: float) -> str:
    """A: externally validated in a large sample; B: large sample OR external;
    C: only small evaluations by the developers; D: no evaluation of the trait."""
    if pd.isna(max_n):
        return "D"
    large = max_cases >= 1000 if max_cases > 0 else max_n >= 10_000
    if ext_pubs and large:
        return "A"
    if ext_pubs or large:
        return "B"
    return "C"


def summarize(scores, perf, ancestry, labels, target_ancestry, locus_checks):
    """One row per (queried trait, PGS) with every criterion useful to pick a score."""
    prim = perf[perf["primary"]]
    rows = []
    for _, s in scores.iterrows():
        pid, trait = s["id"], s["query_trait"]
        ev = prim[(prim["pgs_id"] == pid) & (prim["query_trait"] == trait)]
        other = perf[(perf["pgs_id"] == pid) & (perf["query_trait"] == trait) & ~perf["primary"]]
        ext = ev[ev["external"]]
        per_eval = ev.drop_duplicates("ppm_id")
        ancs = set().union(*map(ancestry_codes, per_eval["ancestry"])) if len(ev) else set()
        dist = ancestry[ancestry["score_id"] == pid].set_index("stage")["dist"]

        def best(metrics, agg="max", ev=ev):
            v = ev.loc[ev["metric"].isin(metrics), "estimate"]
            return getattr(v, agg)() if len(v) else float("nan")

        max_cases = per_eval["n_cases"].max() if len(ev) else 0
        max_n = per_eval["n"].max() if len(ev) else float("nan")
        n_ext_pubs = ext["eval_pub_id"].nunique()
        year = s["publication.date_publication"].year
        row = {
            "query_trait": trait,
            "trait_label": labels.get(trait, trait),
            "pgs_id": pid,
            "name": s["name"],
            "trait_reported": s["trait_reported"],
            "source_traits": s["source_trait"],
            "method": s["method_name"],
            "n_variants": s["variants_number"],
            "dev_publication": f"{s['publication.firstauthor']} ({year})",
            "dev_date": s["publication.date_publication"],
            "dev_journal": s["publication.journal"],
            "dev_doi": s["publication.doi"],
            "gwas_ancestry": _fmt_dist(dist.get("gwas")),
            "dev_ancestry": _fmt_dist(dist.get("dev")),
            "n_evals": per_eval["ppm_id"].nunique(),
            "n_evals_other_traits": other["ppm_id"].nunique(),
            "other_phenotypes": " ; ".join(other["phenotype"].drop_duplicates().head(4)),
            "n_external_pubs": n_ext_pubs,
            "external_by": ", ".join(sorted(set(ext["eval_author"]))),
            "eval_ancestry": ",".join(sorted(ancs)),
            "target_evaluated": target_ancestry in ancs,
            "max_eval_n": max_n,
            "max_eval_cases": max_cases,
            "auc_pgs_only": best(["PGS AUROC (no covariates)"]),
            "auc_with_covariates": best(["AUROC", "C-index"]),
            "or_median": best(["OR"], "median"),
            "hr_median": best(["HR"], "median"),
            "r2_pgs_only": best(["PGS R2 (no covariates)"]),
            "r2_with_covariates": best(["R²"]),
            "harmonized_grch38": isinstance(
                s.get("ftp_harmonized_scoring_files.GRCh38.positions"), str
            ),
            "scoring_file": s["ftp_scoring_file"],
            "evidence": evidence_level(n_ext_pubs, max_cases, max_n),
        }
        for name, df in locus_checks.items():
            hit = df[df["pgs_id"] == pid]
            row[f"{name}_snps"] = hit["snps_found"].iloc[0] if len(hit) else ""
            row[f"{name}_region_variants"] = hit["n_region"].iloc[0] if len(hit) else None
        row["notes"] = "; ".join(_notes(row, target_ancestry))
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["query_trait", "evidence", "target_evaluated", "n_external_pubs", "max_eval_cases",
         "max_eval_n"],
        ascending=[True, True, False, False, False, False],
    ).reset_index(drop=True)


def _notes(r, target):
    notes = []
    if r["n_evals"] == 0:
        notes.append("never evaluated on the trait")
    elif r["n_external_pubs"] == 0:
        notes.append("only evaluated by its authors")
    if r["n_evals"] and not r["target_evaluated"]:
        notes.append(f"never evaluated in {target}")
    if pd.isna(r["auc_pgs_only"]) and pd.notna(r["auc_with_covariates"]):
        notes.append("AUC only with covariates")
    if 0 < r["max_eval_cases"] < 200:
        notes.append(f"at most {int(r['max_eval_cases'])} cases evaluated")
    if r["n_evals_other_traits"]:
        notes.append(f"{r['n_evals_other_traits']} evals on other phenotypes")
    if not r["harmonized_grch38"]:
        notes.append("no harmonized GRCh38 file")
    return notes


def head_to_head(perf, scores) -> pd.DataFrame:
    """Scores compared within the same publication, phenotype, sample and metric:
    the only setting where their metrics are directly comparable."""
    p = perf[perf["primary"] & perf["metric"].isin(HEAD_TO_HEAD_METRICS)].copy()
    keys = ["query_trait", "eval_pub_id", "eval_author", "phenotype", "ancestry", "cohorts",
            "metric"]
    p["cohorts"] = p["cohorts"].fillna("")
    p["n_scores"] = p.groupby(keys)["pgs_id"].transform("nunique")
    p = p[p["n_scores"] >= 2]
    p = p.merge(
        scores[["id", "query_trait", "name"]].rename(columns={"id": "pgs_id"}),
        on=["pgs_id", "query_trait"],
    )
    p["rank"] = p.groupby(keys)["estimate"].rank(ascending=False, method="min").astype(int)
    cols = keys + ["n", "n_cases", "pgs_id", "name", "estimate", "ci_lower", "ci_upper", "rank"]
    return p[cols].sort_values(keys + ["rank"]).reset_index(drop=True)
