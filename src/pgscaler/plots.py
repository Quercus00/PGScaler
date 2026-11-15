"""Static figures for the report (one set per queried trait)."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from analysis import ancestry_codes
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullFormatter

PLAIN = FuncFormatter(lambda v, _: f"{v:g}")

SURFACE, INK, INK_2, GRID = "#ffffff", "#0b0b0b", "#52514e", "#e6e5e1"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # categorical slots 1-3 (all-pairs safe)
OTHER = "#9a9893"
# Sequential (one hue, dark = strongest evidence)
EVIDENCE = {"A": "#0d3a75", "B": "#2a78d6", "C": "#86b6ec", "D": "#bcd4f2"}
ANCESTRY_COLORS = {"EUR": SERIES[0], "EAS": SERIES[1], "AFR": SERIES[2]}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
})


def _label(row) -> str:
    name = str(row["name"])
    return f"{row['pgs_id']}  {name[:28] + '…' if len(name) > 29 else name}"


def _anc_group(text: str) -> str:
    codes = ancestry_codes(text)
    if len(codes) == 1 and (c := next(iter(codes))) in ANCESTRY_COLORS:
        return c
    return "Other / mixed"


def _anc_legend(groups):
    order = [g for g in ["EUR", "EAS", "AFR", "Other / mixed"] if g in groups]
    handles = [Line2D([], [], marker="o", ls="", color=ANCESTRY_COLORS.get(g, OTHER), label=g)
               for g in order]
    handles += [
        Line2D([], [], marker="o", ls="", color=INK_2, label="external validation"),
        Line2D([], [], marker="o", ls="", mfc="none", color=INK_2, label="by the authors"),
    ]
    return handles


def all_scores_interactive(summary: pd.DataFrame, div_id: str = "all-scores") -> str:
    """Interactive plotly figure of every score: publication date x trait, colored by
    evidence level. Returns an HTML <div> embedding plotly.js (works offline)."""
    import plotly.express as px

    s = summary.copy()
    traits = list(dict.fromkeys(s["trait_label"]))
    # Spread scores of the same trait vertically so that points published together don't overlap
    s = s.sort_values("dev_date")
    rank = s.groupby("trait_label").cumcount()
    lanes = s.groupby("trait_label")["pgs_id"].transform("size").clip(upper=7)  # 1 lane if alone
    offset = ((rank % lanes) - (lanes - 1) / 2) * 0.1
    s["y"] = s["trait_label"].map({t: i for i, t in enumerate(traits)}) + offset
    s["Level"] = s["evidence"]
    s["notes"] = s["notes"].fillna("").str.replace("; ", "<br>")
    s["n_variants_txt"] = s["n_variants"].map(lambda v: f"{v:,.0f}" if pd.notna(v) else "?")
    # Marker diameter grows with log10(number of variants) (scores range from a few SNPs to
    # millions), rescaled between the smallest and largest score shown: 7 to 22 px
    log_n = np.log10(s["n_variants"].clip(lower=1))
    lo, hi = log_n.min(), log_n.max()
    rel = (log_n - lo) / (hi - lo) if hi > lo else log_n * 0 + 0.3
    px_size = dict(zip(s["pgs_id"], (7 + 15 * rel).fillna(7)))
    fig = px.scatter(
        s, x="dev_date", y="y", color="Level",
        category_orders={"Level": list(EVIDENCE)}, color_discrete_map=EVIDENCE,
        custom_data=["pgs_id", "name", "trait_reported", "dev_publication", "n_variants_txt",
                     "evidence", "n_evals", "n_external_pubs", "eval_ancestry", "notes"],
    )
    for trace in fig.data:  # one trace per evidence level
        trace.marker.size = [px_size[c[0]] for c in trace.customdata]
    fig.update_traces(
        marker={"opacity": 0.85, "line": {"width": 1, "color": INK_2}},
        hovertemplate=(
            "<b>%{customdata[0]}</b> — %{customdata[1]}<br>"
            "%{customdata[2]}<br>%{customdata[3]} · %{customdata[4]} variants<br>"
            "Level <b>%{customdata[5]}</b> · %{customdata[6]} evals · "
            "%{customdata[7]} external pubs<br>Ancestries: %{customdata[8]}<br>"
            "<i>%{customdata[9]}</i><extra></extra>"
        ),
    )
    fig.update_layout(
        template="plotly_white",
        height=160 + 70 * len(traits), margin={"l": 10, "r": 10, "t": 10, "b": 40},
        paper_bgcolor="#ffffff", plot_bgcolor="#ffffff",
        font={"family": "system-ui, -apple-system, Segoe UI, sans-serif", "color": INK_2},
        legend={"title": f"Evidence level — marker size: number of variants, log scale "
                         f"({10 ** lo:,.0f} to {10 ** hi:,.0f})",
                "orientation": "h", "y": 1.02, "yanchor": "bottom", "x": 0,
                "itemsizing": "constant"},
        hoverlabel={"align": "left", "bgcolor": "#ffffff", "font": {"color": INK}},
        xaxis={"title": "Publication date", "gridcolor": GRID},
        yaxis={"title": None, "tickvals": list(range(len(traits))), "ticktext": traits,
               "range": [len(traits) - 0.4, -0.6], "gridcolor": GRID, "zeroline": False},
    )
    return fig.to_html(full_html=False, include_plotlyjs=True, div_id=div_id,
                       config={"displaylogo": False, "responsive": True})


def overview_and_years(summary: pd.DataFrame, path) -> str:
    """Overview (score size vs. evaluations, left) and publications per year (right), with
    one shared evidence-level legend. The overview is left out when no score has been
    evaluated on the trait (it would be empty). Returns the figure title."""
    has_overview = bool((summary["n_evals"] > 0).any())
    if has_overview:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.5), width_ratios=[3, 2])
        _overview(summary, ax1)
    else:
        fig, ax2 = plt.subplots(figsize=(8, 4.5))
    _by_year(summary, ax2)
    levels = [lvl for lvl in EVIDENCE if (summary["evidence"] == lvl).any()]
    handles = [Line2D([], [], marker="s", ls="", ms=10, color=EVIDENCE[lvl], label=f"Level {lvl}")
               for lvl in levels]
    title = "Evidence level" + (" (overview: marker size ∝ √(max cases, or N/10))"
                                if has_overview else "")
    fig.legend(handles=handles, title=title, loc="lower center", ncol=len(levels),
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=3)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return "Overview and publications per year" if has_overview else "Publications per year"


def _overview(summary: pd.DataFrame, ax):
    """Each score: size of the score vs. how much it has been evaluated."""
    s = summary.copy()
    size = np.where(s["max_eval_cases"] > 0, s["max_eval_cases"], s["max_eval_n"] / 10)
    s["size"] = pd.Series(size, index=s.index).fillna(0)
    scale = s["size"].max() if s["size"].max() > 0 else 1
    for lvl, color in EVIDENCE.items():
        d = s[(s["evidence"] == lvl) & s["n_variants"].notna()]
        if len(d):
            ax.scatter(d["n_variants"], d["n_evals"], s=20 + 180 * np.sqrt(d["size"] / scale),
                       color=color, edgecolor=INK_2, linewidth=0.6, zorder=3)
    for _, r in s[s["evidence"] == "A"].iterrows():
        ax.annotate(r["pgs_id"], (r["n_variants"], r["n_evals"]), xytext=(5, 4),
                    textcoords="offset points", fontsize=8, color=INK_2)
    ax.set_xscale("log")
    ax.set_xlabel("Number of variants in the score (log)")
    ax.set_ylabel("Evaluations on the trait")
    ax.set_title("Score size vs. evaluations", loc="left")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_axisbelow(True)


def forest(perf: pd.DataFrame, summary: pd.DataFrame, metrics, xlabel, path, log=True, ref=1.0):
    """One point per evaluation (with CI), colored by ancestry, filled when external."""
    p = perf[perf["primary"] & perf["metric"].isin(metrics)]
    p = p.merge(summary[["pgs_id", "name"]], on="pgs_id")
    if p.empty:
        return False
    order = summary[summary["pgs_id"].isin(p["pgs_id"])]
    labels = {r["pgs_id"]: _label(r) for _, r in order.iterrows()}
    ypos = {pid: i for i, pid in enumerate(reversed(order["pgs_id"].tolist()))}
    fig, ax = plt.subplots(figsize=(9, 1.2 + 0.32 * len(ypos)))
    p = p.assign(grp=p["ancestry"].map(_anc_group))
    jitter = p.groupby("pgs_id").cumcount()
    count = p.groupby("pgs_id")["pgs_id"].transform("size")
    y = p["pgs_id"].map(ypos) + (jitter - (count - 1) / 2) * min(0.6 / count.max(), 0.12)
    for (_, r), yy in zip(p.iterrows(), y):
        c = ANCESTRY_COLORS.get(r["grp"], OTHER)
        if pd.notna(r["ci_lower"]) and pd.notna(r["ci_upper"]):
            ax.plot([r["ci_lower"], r["ci_upper"]], [yy, yy], color=c, lw=1.2, alpha=0.7)
        ax.plot(r["estimate"], yy, "o", ms=6, color=c, mfc=c if r["external"] else SURFACE,
                mew=1.5, zorder=3)
    if ref is not None:
        ax.axvline(ref, color=INK_2, lw=1, ls="--")
    if log:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(PLAIN)
        ax.xaxis.set_minor_formatter(PLAIN)
        if ax.get_xlim()[1] / ax.get_xlim()[0] > 20:
            ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(list(ypos.values()), [labels[k] for k in ypos])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(xlabel)
    ax.set_title(f"{' / '.join(metrics)} per evaluation", loc="left")
    ax.legend(handles=_anc_legend(set(p["grp"])), loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def _by_year(summary: pd.DataFrame, ax):
    """Scores published per year, stacked by evidence level."""
    s = summary.assign(year=summary["dev_publication"].str.extract(r"\((\d{4})\)")[0].astype(int))
    tab = s.pivot_table(index="year", columns="evidence", values="pgs_id", aggfunc="count",
                        fill_value=0)
    tab = tab.reindex(range(tab.index.min(), tab.index.max() + 1), fill_value=0)
    bottom = np.zeros(len(tab))
    for lvl in [lv for lv in EVIDENCE if lv in tab.columns]:
        ax.bar(tab.index, tab[lvl], bottom=bottom, color=EVIDENCE[lvl], width=0.7,
               edgecolor=SURFACE, linewidth=2)
        bottom += tab[lvl].to_numpy()
    ax.set_ylabel("Scores published")
    ax.grid(axis="x", visible=False)
    ax.set_title("Scores per publication year", loc="left")
    years = tab.index
    step = max(1, len(years) // 12)
    ax.set_xticks(years[::step])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
    mid, half = (years.min() + years.max()) / 2, max((years.max() - years.min()) / 2 + 0.6, 2.5)
    ax.set_xlim(mid - half, mid + half)  # a single year doesn't fill the whole width
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_axisbelow(True)
