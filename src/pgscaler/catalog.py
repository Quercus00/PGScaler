"""Fetch data from the PGS Catalog: scores, performance metrics, scoring-file content."""

from pathlib import Path

import pandas as pd
from pandaspgs import get_scores
from pandaspgs.client import get_performance, get_trait

PERF_URL = "https://www.pgscatalog.org/rest/performance/search?pgs_id={}"
TRAIT_URL = "https://www.pgscatalog.org/rest/trait/{}?include_children=1"


def fetch_trait(trait_id: str) -> dict:
    """Label and child traits (ontology descendants that have scores) of a trait.
    A PGS id is returned as a pseudo-trait without children."""
    if trait_id.startswith("PGS"):
        return {"id": trait_id, "label": trait_id, "children": []}
    res = get_trait(TRAIT_URL.format(trait_id), cached=True)
    res = (res[0] if res else {}) if isinstance(res, list) else (res or {})
    children = [c["id"] if isinstance(c, dict) else c for c in res.get("child_traits") or []]
    return {"id": trait_id, "label": res.get("label", trait_id), "children": children}


def group_traits(infos: dict[str, dict], explicit: dict[str, list[str]], auto: bool = True):
    """Group the queried traits.

    - explicit groups ({name: [trait ids]}) take precedence;
    - otherwise, if `auto`, a trait that is a descendant of another queried trait joins it
      (e.g. bipolar I disorder -> bipolar disorder); the top-most ancestor wins.

    Returns {trait_id: group_key} and {group_key: {"label", "members"}}.
    """
    group_of = {t: name for name, members in explicit.items() for t in members if t in infos}
    in_explicit = set(group_of)
    for t in infos:
        if t in in_explicit:
            continue
        parents = [p for p, info in infos.items()
                   if auto and p != t and t in info["children"] and p not in in_explicit]
        root = max(parents, key=lambda p: len(infos[p]["children"])) if parents else t
        group_of[t] = root
    groups = {}
    for t, g in group_of.items():
        label = g if g in explicit else infos[g]["label"]
        groups.setdefault(g, {"label": label, "members": []})["members"].append(t)
    return group_of, groups


def fetch_scores(infos: dict[str, dict], group_of: dict[str, str], include_children: bool):
    """Scores (one row per PGS and trait group) and their ancestry distribution
    (gwas/dev/eval stages). `source_trait` lists the trait ids the score was found under."""
    scores, ancestry = [], []
    for trait, info in infos.items():
        ids = [trait] + (info["children"] if include_children else [])
        for q in ids:
            res = get_scores(pgs_id=q) if q.startswith("PGS") else get_scores(trait_id=q)
            if res.scores is None or res.scores.empty:
                continue
            s = res.scores.copy()
            s["query_trait"] = group_of[trait]
            s["source_trait"] = q
            scores.append(s)
            ancestry.append(res.ancestry_distribution)
    scores = pd.concat(scores, ignore_index=True)
    sources = scores.groupby(["query_trait", "id"])["source_trait"].agg(
        lambda x: ", ".join(dict.fromkeys(x)))
    scores = scores.drop_duplicates(["query_trait", "id"]).drop(columns="source_trait")
    scores = scores.merge(sources.reset_index(), on=["query_trait", "id"])
    scores["variants_number"] = pd.to_numeric(scores["variants_number"], errors="coerce")
    scores["publication.date_publication"] = pd.to_datetime(
        scores["publication.date_publication"], errors="coerce"
    )
    ancestry = pd.concat(ancestry, ignore_index=True).drop_duplicates(["score_id", "stage"])
    return scores, ancestry


def filter_scores(scores: pd.DataFrame, date_min=None, date_max=None, n_snp=None):
    """Keep scores published between date_min and date_max (years, inclusive) and with at
    least n_snp variants. A None bound is not applied. Scores with an unknown date or number
    of variants are dropped when the corresponding filter is set."""
    year = scores["publication.date_publication"].dt.year
    keep = pd.Series(True, index=scores.index)
    if date_min is not None:
        keep &= year >= date_min
    if date_max is not None:
        keep &= year <= date_max
    if n_snp is not None:
        keep &= scores["variants_number"] >= n_snp
    return scores[keep].reset_index(drop=True)


def fetch_performance(pgs_ids) -> pd.DataFrame:
    """Long table: one row per (evaluation PPM, metric).

    Parses the raw REST JSON (through the pandaspgs client and its cache) because
    pandaspgs.get_performances() crashes when none of the samples list a cohort.
    """
    rows = []
    for pgs_id in pgs_ids:
        for rec in get_performance(PERF_URL.format(pgs_id), cached=True) or []:
            samples = rec["sampleset"]["samples"]
            base = {
                "ppm_id": rec["id"],
                "pgs_id": rec["associated_pgs_id"],
                "phenotype": rec["phenotyping_reported"],
                "covariates": rec["covariates"],
                "eval_pub_id": rec["publication"]["id"],
                "eval_author": rec["publication"]["firstauthor"],
                "eval_date": rec["publication"]["date_publication"],
                "n": sum(s["sample_number"] or 0 for s in samples),
                "n_cases": sum(s["sample_cases"] or 0 for s in samples),
                "ancestry": " | ".join(sorted({s["ancestry_broad"] or "NR" for s in samples})),
                "cohorts": ", ".join(
                    sorted({c["name_short"] for s in samples for c in s["cohorts"]})
                ),
            }
            for kind, metrics in rec["performance_metrics"].items():
                for m in metrics:
                    rows.append(
                        {
                            **base,
                            "kind": kind,
                            "metric": " ".join(m["name_short"].split()),
                            "metric_long": m["name_long"],
                            "estimate": m["estimate"],
                            "ci_lower": m.get("ci_lower"),
                            "ci_upper": m.get("ci_upper"),
                            "se": m.get("se"),
                        }
                    )
    columns = ["ppm_id", "pgs_id", "phenotype", "covariates", "eval_pub_id", "eval_author",
               "eval_date", "n", "n_cases", "ancestry", "cohorts", "kind", "metric",
               "metric_long", "estimate", "ci_lower", "ci_upper", "se"]
    return pd.DataFrame(rows, columns=columns)  # columns kept even when no evaluation


def check_locus(pgs_ids, scores: pd.DataFrame, locus: dict, cache_file: Path) -> pd.DataFrame:
    """Does each score contain a key locus (e.g. APOE for Alzheimer)?

    `locus` = {name, chr, snps: {rsid: GRCh38 position}, window}. Reports which of
    the key SNPs are in the score and how many score variants fall within
    +/- window of them. Streams the harmonized GRCh38 scoring files; results are
    cached in `cache_file` (genome-wide scores weigh tens of MB each).
    """
    cols = ["pgs_id", "locus", "snps_found", "n_region"]
    cache = pd.read_csv(cache_file) if cache_file.exists() else pd.DataFrame(columns=cols)
    done = set(cache.loc[cache["locus"] == locus["name"], "pgs_id"])
    urls = scores.drop_duplicates("id").set_index("id")[
        "ftp_harmonized_scoring_files.GRCh38.positions"
    ]
    chrom, snps = str(locus["chr"]), dict(locus["snps"])
    lo, hi = min(snps.values()) - locus["window"], max(snps.values()) + locus["window"]
    new = []
    for pgs_id in pgs_ids:
        url = urls.get(pgs_id)
        if pgs_id in done or not isinstance(url, str):
            continue
        found, n_region = set(), 0
        for chunk in pd.read_csv(
            url, sep="\t", comment="#", compression="gzip", chunksize=200_000, dtype=str,
            usecols=lambda c: c in {"rsID", "hm_rsID", "hm_chr", "hm_pos"},
        ):
            chunk = chunk[chunk["hm_chr"] == chrom]
            pos = pd.to_numeric(chunk["hm_pos"], errors="coerce")
            n_region += int(pos.between(lo, hi).sum())
            ids = set(chunk.get("hm_rsID", [])) | set(chunk.get("rsID", []))
            found |= {snp for snp, p in snps.items() if snp in ids or (pos == p).any()}
        new.append({"pgs_id": pgs_id, "locus": locus["name"],
                    "snps_found": ",".join(sorted(found)) or "none", "n_region": n_region})
        print(f"  {pgs_id}: {locus['name']} SNPs={new[-1]['snps_found']}, region={n_region}")
    if new:
        cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(cache_file, index=False)
    return cache[cache["locus"] == locus["name"]]
