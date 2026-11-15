import re
from pathlib import Path

import analysis
import catalog
import hydra
import pandas as pd
import plots
import report
from omegaconf import DictConfig, OmegaConf
from utils import print_found

REPO_ROOT = Path(__file__).resolve().parents[2]  # src/pgscaler/main.py -> repo root

FOREST_PLOTS = [
    # (metrics, title, x label, log scale, reference line)
    (["OR"], "Odds ratio per evaluation", "OR (log)", True, 1.0),
    (["HR"], "Hazard ratio per evaluation", "HR (log)", True, 1.0),
    (["PGS AUROC (no covariates)"], "AUC of the PGS alone", "AUC", False, 0.5),
    (["AUROC", "C-index"], "AUC / C-index with covariates", "AUC", False, 0.5),
    (["PGS R2 (no covariates)"], "R² of the PGS alone", "R²", False, None),
]


@hydra.main(version_base=None, config_path="../config", config_name="config")
def main(cfg: DictConfig):
    # A relative outputpath is relative to the repository root, whatever the working directory
    out = Path(cfg.outputpath)
    if not out.is_absolute():
        out = REPO_ROOT / out
    (out / "figures").mkdir(parents=True, exist_ok=True)
    print(OmegaConf.to_yaml(cfg))
    infos = {t: catalog.fetch_trait(t) for t in cfg.trait}
    group_of, groups = catalog.group_traits(
        infos, OmegaConf.to_container(cfg.get("trait_groups") or {}),
        auto=cfg.get("group_child_traits", True),
    )
    traits = list(groups)
    labels = {g: info["label"] for g, info in groups.items()}
    aliases = {g: [info["label"]] + [infos[t]["label"] for t in info["members"]]
               for g, info in groups.items()}

    scores, ancestry = catalog.fetch_scores(
        infos, group_of, cfg.get("include_child_traits", True))
    found = scores.groupby("query_trait").size()
    scores = catalog.filter_scores(
        scores, cfg.get("date_min"), cfg.get("date_max"), cfg.get("n_snp"))
    kept = scores.groupby("query_trait").size()
    print(f"Filters: publication years {cfg.get('date_min')}–{cfg.get('date_max')}, "
          f"at least {cfg.get('n_snp')} variants")
    for g, info in groups.items():
        members = ", ".join(f"{t} ({infos[t]['label']})" for t in info["members"])
        print_found(f"{info['label']} [{members}]: {found.get(g, 0)} scores found, "
                    f"{kept.get(g, 0)} kept")
    if scores.empty:
        print("No score left after filtering: nothing to report.")
        return

    # Per-trait options may be keyed by group name or by any trait id of the group
    def by_group(option):
        out = {}
        for key, value in (OmegaConf.to_container(cfg.get(option) or {})).items():
            g = key if key in groups else group_of.get(key)
            if g is not None:
                out.setdefault(g, []).append(value)
        return out

    patterns = {g: "|".join(f"(?:{p})" for p in ps) for g, ps in by_group("eval_patterns").items()}
    perf = catalog.fetch_performance(scores["id"].unique())
    perf = analysis.annotate_performance(perf, scores, aliases, patterns)
    print_found(f"{perf.ppm_id.nunique()} evaluations, {perf[perf.primary].ppm_id.nunique()} "
                "on the queried traits")

    locus_checks = {}
    for trait, loci_lists in by_group("key_loci").items():
        loci = [locus for lst in loci_lists for locus in lst]
        ids = scores.loc[scores.query_trait == trait, "id"]
        for locus in loci:
            print(f"Checking {locus['name']} in {len(ids)} scoring files (cached)...")
            res = catalog.check_locus(ids, scores, locus, out / "cache" / "key_loci.csv")
            locus_checks[locus["name"]] = pd.concat(
                [locus_checks.get(locus["name"]), res], ignore_index=True
            )

    summary = analysis.summarize(scores, perf, ancestry, labels, cfg.target_ancestry,
                                 locus_checks)
    h2h = analysis.head_to_head(perf, scores)

    summary.to_csv(out / "scores_summary.csv", index=False)
    perf.to_csv(out / "performance_metrics.csv", index=False)
    h2h.to_csv(out / "head_to_head.csv", index=False)

    figures = {}
    for trait in traits:
        s = summary[summary.query_trait == trait]
        p = perf[perf.query_trait == trait]
        if s.empty:
            continue
        stem = re.sub(r"\W+", "_", trait)
        figs = []
        path = out / "figures" / f"{stem}_overview.png"
        figs.append((plots.overview_and_years(s, path), path))
        for metrics, title, xlabel, log, ref in FOREST_PLOTS:
            path = out / "figures" / f"{stem}_{'_'.join(metrics)}.png".replace(" ", "")
            if plots.forest(p, s, metrics, xlabel, path, log=log, ref=ref):
                figs.append((title, path))
        figures[trait] = figs

    groups_html = {g: {"label": info["label"],
                       "members": [(t, infos[t]["label"]) for t in info["members"]]}
                   for g, info in groups.items()}
    report.build(summary, h2h, figures, cfg, out / "report.html",
                 overview_html=plots.all_scores_interactive(summary), groups=groups_html)
    print_found(f"Report written to {(out / 'report.html').resolve()}")


if __name__ == "__main__":
    main()
