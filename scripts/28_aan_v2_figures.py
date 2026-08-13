"""Generate the final non-writing V2 scientific figure package.

Every plotted value is read from a tracked V2 result table. Release-value
assertions intentionally stop the script if a primary source table changes.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import os
import tempfile
import textwrap
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pd_envtox_matplotlib")
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v2" / "aan_figures"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE_PATHS = [
    "results/v2/human_validation/frozen21_human_PD_summary.tsv",
    "results/v2/human_validation/frozen21_matched_null_summary.tsv",
    "results/v2/gse46798/GSE46798_module_tests.tsv",
    "results/v2/gse17542_sn_vta/GSE17542_SN_VTA_module_tests.tsv",
    "results/v2/continuous_network_test/continuous_network_predictive_summary.tsv",
    "results/v2/reactome_convergence/Reactome_targeted_PD_family.tsv",
    "results/v2/reactome_convergence/Reactome_strict_three_model_recurrence.tsv",
    "results/v2/gse116280_rotenone/GSE116280_proteasome_all_contrasts.tsv",
    "results/v2/gse150005_multitoxicant/GSE150005_proteasome_CAMERA_summary.tsv",
    "results/v2/gse196190_heldout/GSE196190_primary_mechanism_result.tsv",
    "results/v2/heldout_external_sensitivity/external_sensitivity_primary_union_results.tsv",
    "results/v2/heldout_external_sensitivity/GSE4773_primary_union_sensitivity.tsv",
    "results/v2/human_sn_proteasome/donor_pathway_scores.tsv",
    "results/v2/human_sn_proteasome/pathway_model_results.tsv",
    "results/v2/human_sn_proteasome/validation_summary.tsv",
]

NAVY = "#17324D"
BLUE = "#2878B5"
TEAL = "#11857A"
GOLD = "#D89B2B"
PURPLE = "#6750A4"
SLATE = "#64748B"
MID_GREY = "#A8B2BD"
LIGHT_GREY = "#E8EDF2"
PALE_BLUE = "#EAF3F8"
PALE_TEAL = "#DDF3EF"
PALE_GOLD = "#FAEBC8"
PALE_PURPLE = "#EDE8F8"
PALE_GREY = "#EEF1F4"
TEXT = "#17212B"
WHITE = "#FFFFFF"

FIXED_DATE = dt.datetime(2026, 8, 13, 12, 0, 0, tzinfo=dt.timezone.utc)
PDF_METADATA = {
    "Author": "pd-envtox",
    "Creator": "scripts/28_aan_v2_figures.py",
    "CreationDate": FIXED_DATE,
    "ModDate": FIXED_DATE,
}
PNG_METADATA = {"Software": "scripts/28_aan_v2_figures.py"}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.labelcolor": TEXT,
            "axes.edgecolor": SLATE,
            "axes.linewidth": 0.8,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "text.color": TEXT,
            "figure.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tsv(relative_path: str) -> pd.DataFrame:
    path = ROOT / relative_path
    if not path.exists():
        raise FileNotFoundError(f"Missing figure source: {relative_path}")
    return pd.read_csv(path, sep="\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    clean = frame.copy()
    object_columns = clean.select_dtypes(include="object").columns
    clean[object_columns] = clean[object_columns].fillna("NA").replace("", "NA")
    clean.to_csv(path, sep="\t", index=False, na_rep="NA")


def assert_close(observed: float, expected: float, label: str, tol: float = 1e-10) -> None:
    if not np.isclose(observed, expected, atol=tol, rtol=tol):
        raise AssertionError(f"{label}: expected {expected}, observed {observed}")


def format_p(value: float) -> str:
    if value <= 0:
        return "<0.001"
    if value < 0.001:
        return f"{value:.1e}"
    if value < 0.1:
        return f"{value:.4f}"
    return f"{value:.3f}"


def clean_axes(ax: mpl.axes.Axes, keep_left: bool = True, keep_bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    ax.spines["bottom"].set_visible(keep_bottom)


def save_figure(fig: mpl.figure.Figure, stem: str, title: str) -> None:
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.08,
        metadata=PNG_METADATA,
    )
    pdf_meta = dict(PDF_METADATA)
    pdf_meta["Title"] = title
    fig.savefig(
        pdf,
        bbox_inches="tight",
        pad_inches=0.08,
        metadata=pdf_meta,
    )
    plt.close(fig)


def release_values() -> dict[str, object]:
    human_panel = read_tsv(SOURCE_PATHS[0]).iloc[0]
    human_null = read_tsv(SOURCE_PATHS[1]).iloc[0]
    continuous = read_tsv(SOURCE_PATHS[4])
    reactome = read_tsv(SOURCE_PATHS[5])
    strict = read_tsv(SOURCE_PATHS[6])
    rotenone = read_tsv(SOURCE_PATHS[7])
    multitoxicant = read_tsv(SOURCE_PATHS[8])
    heldout = read_tsv(SOURCE_PATHS[9]).iloc[0]
    external = read_tsv(SOURCE_PATHS[10])
    gse4773 = read_tsv(SOURCE_PATHS[11])
    donors = read_tsv(SOURCE_PATHS[12])
    human_models = read_tsv(SOURCE_PATHS[13])
    human_summary = read_tsv(SOURCE_PATHS[14]).iloc[0]

    proteasome = reactome.loc[reactome["pathway"] == "Proteasome assembly"].iloc[0]
    strict_proteasome = strict.loc[strict["pathway"] == "Proteasome assembly"].iloc[0]
    primary_rotenone = rotenone.loc[rotenone["contrast"] == "D15_50nM_24h"].iloc[0]
    multitox_rot = multitoxicant.loc[multitoxicant["Contrast"] == "Rotenone"].iloc[0]
    human_primary = human_models.loc[
        human_models["analysis"] == "GSE178265_primary_adjusted"
    ].iloc[0]
    human_replication = human_models.loc[
        human_models["analysis"] == "GSE243639_replication_adjusted"
    ].iloc[0]

    assert int(human_panel["meta_FDR_sig_n"]) == 0
    assert_close(
        float(human_null["empirical_p_direction_agreement_upper"]),
        0.19738026197380262,
        "human panel matched-null p",
    )
    assert_close(
        float(proteasome["meta_FDR"]),
        0.002029077321127518,
        "proteasome cross-model meta-FDR",
    )
    assert_close(
        float(primary_rotenone["nominal_p"]),
        0.09578544061302682,
        "GSE116280 frozen primary p",
    )
    assert_close(
        float(multitox_rot["PValue_one_sided_down"]),
        0.240082166890404,
        "GSE150005 rotenone one-sided p",
    )
    assert str(heldout["frozen_outcome"]) == "NOT_CONFIRMED"
    assert_close(
        float(heldout["CAMERA_p_one_sided_down"]),
        0.122137297655125,
        "GSE196190 frozen primary p",
    )
    assert str(human_summary["overall_outcome"]) == "PRIMARY_ONLY"
    assert_close(
        float(human_primary["p_one_sided_down"]),
        0.03522407025103581,
        "GSE178265 one-sided p",
    )
    assert_close(
        float(human_replication["p_one_sided_down"]),
        0.5095340240511492,
        "GSE243639 one-sided p",
    )

    primary_donors = donors.loc[donors["analysis"] == "GSE178265_primary_adjusted"]
    replication_donors = donors.loc[
        donors["analysis"] == "GSE243639_replication_adjusted"
    ]
    assert (int((primary_donors["pd"] == 1).sum()), int((primary_donors["pd"] == 0).sum())) == (6, 8)
    assert (
        int((replication_donors["pd"] == 1).sum()),
        int((replication_donors["pd"] == 0).sum()),
    ) == (7, 8)

    return {
        "human_panel": human_panel,
        "human_null": human_null,
        "continuous": continuous,
        "proteasome": proteasome,
        "strict_proteasome": strict_proteasome,
        "rotenone": rotenone,
        "primary_rotenone": primary_rotenone,
        "multitoxicant": multitoxicant,
        "multitox_rot": multitox_rot,
        "heldout": heldout,
        "external": external,
        "gse4773": gse4773,
        "donors": donors,
        "primary_donors": primary_donors,
        "replication_donors": replication_donors,
        "human_models": human_models,
        "human_primary": human_primary,
        "human_replication": human_replication,
        "human_summary": human_summary,
    }


def draw_node(
    ax: mpl.axes.Axes,
    x: float,
    y: float,
    text: str,
    fill: str,
    edge: str,
    linestyle: str = "-",
    width: float = 16.0,
    height: float = 10.0,
) -> None:
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.45,rounding_size=1.3",
        linewidth=1.4,
        edgecolor=edge,
        facecolor=fill,
        linestyle=linestyle,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=8.7, linespacing=1.25, zorder=3)


def draw_arrow(ax: mpl.axes.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.1,
        color=MID_GREY,
        shrinkA=3,
        shrinkB=3,
        zorder=1,
    )
    ax.add_patch(arrow)


def figure_1_evidence_design(v: dict[str, object]) -> pd.DataFrame:
    proteasome = v["proteasome"]
    primary_rotenone = v["primary_rotenone"]
    heldout = v["heldout"]
    human_primary = v["human_primary"]
    human_replication = v["human_replication"]
    human_panel = v["human_panel"]
    human_null = v["human_null"]
    continuous = v["continuous"]
    max_partial = continuous.iloc[continuous["partial_spearman_r"].abs().argmax()]

    rows = [
        ("Fixed panel", "Original hypothesis", "21-gene panel", "Original", ""),
        (
            "Fixed panel",
            "External reanalysis",
            "GSE46798 + GSE17542: not confirmed",
            "Not confirmed",
            "0 FDR-significant panel genes",
        ),
        (
            "Fixed panel",
            "Human disease evidence",
            "GSE20141 + GSE7621: not confirmed",
            "Not confirmed",
            f"meta-FDR genes = {int(human_panel['meta_FDR_sig_n'])}; matched-null p = {human_null['empirical_p_direction_agreement_upper']:.3f}",
        ),
        ("Network score", "Original hypothesis", "Continuous proximity score", "Original", ""),
        (
            "Network score",
            "External reanalysis",
            "Negligible, non-generalizing association",
            "Mixed/weak",
            f"maximum |partial r| = {abs(max_partial['partial_spearman_r']):.3f}; 3 of 4 outcomes null",
        ),
        (
            "Proteasome",
            "Retrospective convergence",
            "Three-model Reactome convergence",
            "Retrospective signal",
            f"meta-FDR = {proteasome['meta_FDR']:.5f}",
        ),
        (
            "Proteasome",
            "Frozen external test",
            "GSE116280 mature LUHMES",
            "Not confirmed",
            f"one-sided-down p = {primary_rotenone['nominal_p']:.4f}",
        ),
        (
            "Proteasome",
            "Human primary",
            "GSE178265 DA nuclei",
            "Confirmed",
            f"one-sided-down p = {human_primary['p_one_sided_down']:.4f}",
        ),
        (
            "Proteasome",
            "Independent replication",
            "GSE243639 DA neurons",
            "Not confirmed",
            f"one-sided-down p = {human_replication['p_one_sided_down']:.4f}",
        ),
        (
            "Cilium/Hedgehog/trafficking",
            "Retrospective convergence",
            "GSE150005 nomination",
            "Nominated",
            "retrospective only",
        ),
        (
            "Cilium/Hedgehog/trafficking",
            "Frozen external test",
            "GSE196190 MPP+",
            "Not confirmed",
            f"one-sided-down p = {heldout['CAMERA_p_one_sided_down']:.4f}",
        ),
        (
            "Cilium/Hedgehog/trafficking",
            "Independent replication",
            "External LUHMES sensitivities",
            "Not confirmed",
            "GSE229460 p = 0.401; GSE287941 p = 0.392",
        ),
    ]
    flow = pd.DataFrame(
        rows, columns=["mechanism", "stage", "label", "status", "statistic"]
    )
    write_tsv(flow, OUT / "fig1_evidence_flow.tsv")

    fig, ax = plt.subplots(figsize=(12.2, 7.1))
    ax.set_xlim(-5, 105)
    ax.set_ylim(0, 100)
    ax.axis("off")

    columns = {
        "Original": 8,
        "Retrospective": 29,
        "Frozen external": 51,
        "Human primary": 73,
        "Independent": 95,
    }
    headers = {
        "Original": "Original\nhypothesis",
        "Retrospective": "Retrospective\nreanalysis",
        "Frozen external": "Frozen\nexternal test",
        "Human primary": "Human disease\nprimary",
        "Independent": "Independent\nreplication / sensitivity",
    }
    for key, x in columns.items():
        ax.text(x, 95, headers[key], ha="center", va="center", fontsize=9.5, fontweight="bold", color=NAVY)

    lane_y = {"Fixed panel": 78, "Network score": 59, "Proteasome": 38, "Cilium": 16}
    for y in lane_y.values():
        ax.plot([0, 103], [y - 8.3, y - 8.3], color=LIGHT_GREY, lw=0.8, zorder=0)

    draw_node(ax, columns["Original"], lane_y["Fixed panel"], "21-gene panel\nnetwork-prioritized", WHITE, NAVY)
    draw_node(ax, columns["Retrospective"], lane_y["Fixed panel"], "GSE46798 + GSE17542\npanel not confirmed", PALE_GREY, SLATE)
    draw_node(ax, columns["Human primary"], lane_y["Fixed panel"], "GSE20141 + GSE7621\n0 meta-FDR genes", PALE_GREY, SLATE)
    draw_arrow(ax, (16, lane_y["Fixed panel"]), (21, lane_y["Fixed panel"]))
    draw_arrow(ax, (37, lane_y["Fixed panel"]), (65, lane_y["Fixed panel"]))

    draw_node(ax, columns["Original"], lane_y["Network score"], "Continuous network\nproximity score", WHITE, NAVY)
    draw_node(
        ax,
        columns["Retrospective"],
        lane_y["Network score"],
        f"Tiny maximum effect\n|partial r| = {abs(max_partial['partial_spearman_r']):.3f}\n3 of 4 outcomes null",
        PALE_GOLD,
        GOLD,
        linestyle="--",
    )
    draw_arrow(ax, (16, lane_y["Network score"]), (21, lane_y["Network score"]))

    draw_node(
        ax,
        columns["Retrospective"],
        lane_y["Proteasome"],
        f"Proteasome assembly\n3-model convergence\nmeta-FDR = {proteasome['meta_FDR']:.4f}",
        PALE_PURPLE,
        PURPLE,
        linestyle="--",
    )
    draw_node(
        ax,
        columns["Frozen external"],
        lane_y["Proteasome"],
        f"GSE116280\nnot confirmed\np = {primary_rotenone['nominal_p']:.4f}",
        PALE_GREY,
        SLATE,
    )
    draw_node(
        ax,
        columns["Human primary"],
        lane_y["Proteasome"],
        f"GSE178265\nconfirmed\np = {human_primary['p_one_sided_down']:.4f}",
        PALE_TEAL,
        TEAL,
    )
    draw_node(
        ax,
        columns["Independent"],
        lane_y["Proteasome"],
        f"GSE243639\nnot confirmed\np = {human_replication['p_one_sided_down']:.4f}",
        PALE_GREY,
        SLATE,
    )
    draw_arrow(ax, (37, lane_y["Proteasome"]), (43, lane_y["Proteasome"]))
    draw_arrow(ax, (59, lane_y["Proteasome"]), (65, lane_y["Proteasome"]))
    draw_arrow(ax, (81, lane_y["Proteasome"]), (87, lane_y["Proteasome"]))

    draw_node(
        ax,
        columns["Retrospective"],
        lane_y["Cilium"],
        "Cilium / Hedgehog /\ntrafficking nominated\nin GSE150005",
        PALE_PURPLE,
        PURPLE,
        linestyle="--",
    )
    draw_node(
        ax,
        columns["Frozen external"],
        lane_y["Cilium"],
        f"GSE196190\nnot confirmed\np = {heldout['CAMERA_p_one_sided_down']:.4f}",
        PALE_GREY,
        SLATE,
    )
    draw_node(
        ax,
        columns["Independent"],
        lane_y["Cilium"],
        "External LUHMES\nsensitivities\nnot confirmed",
        PALE_GREY,
        SLATE,
        linestyle=":",
    )
    draw_arrow(ax, (37, lane_y["Cilium"]), (43, lane_y["Cilium"]))
    draw_arrow(ax, (59, lane_y["Cilium"]), (87, lane_y["Cilium"]))

    legend = [
        Patch(facecolor=PALE_PURPLE, edgecolor=PURPLE, label="Retrospective or nominated"),
        Patch(facecolor=PALE_TEAL, edgecolor=TEAL, label="Confirmed frozen primary"),
        Patch(facecolor=PALE_GOLD, edgecolor=GOLD, label="Mixed or weak"),
        Patch(facecolor=PALE_GREY, edgecolor=SLATE, label="Not confirmed"),
    ]
    ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.07), ncol=4, frameon=False, fontsize=8.7)
    fig.suptitle("Evidence progression after the original network-panel hypothesis", x=0.5, y=1.01, fontsize=16, fontweight="bold")
    ax.text(0, 99.7, "Frozen tests retain their original outcomes; sensitivity analyses do not relabel them.", fontsize=10, color=SLATE)
    save_figure(fig, "fig1_evidence_design", "Evidence progression after the original network-panel hypothesis")
    return flow


def figure_2_cross_model(v: dict[str, object]) -> pd.DataFrame:
    row = v["proteasome"]
    strict = v["strict_proteasome"]
    values = pd.DataFrame(
        [
            {
                "model": "Human corrected DA neurons\nparaquat/maneb",
                "NES": row["NES_GSE46798_exposure_corrected_limma"],
                "nominal_p": row["nominal_p_GSE46798_exposure_corrected_limma"],
            },
            {
                "model": "Mouse SN dopamine neurons\n10-day MPTP",
                "NES": row["NES_GSE17542_SN_MPTP10_limma"],
                "nominal_p": row["nominal_p_GSE17542_SN_MPTP10_limma"],
            },
            {
                "model": "Human postmortem SN\nPD meta-analysis",
                "NES": row["NES_human_PD_postmortem_meta_limma"],
                "nominal_p": row["nominal_p_human_PD_postmortem_meta_limma"],
            },
        ]
    )
    values["cross_model_meta_z"] = float(row["meta_z"])
    values["cross_model_meta_p"] = float(row["meta_p"])
    values["cross_model_meta_FDR"] = float(row["meta_FDR"])
    values["shared_leading_genes"] = strict["shared_leading_genes_all_3"]
    write_tsv(values, OUT / "fig2_cross_model_values.tsv")

    fig, ax = plt.subplots(figsize=(9.3, 5.5))
    y = np.arange(len(values))[::-1]
    nes = values["NES"].to_numpy(dtype=float)
    ax.axvspan(-3.6, 0, color=PALE_BLUE, alpha=0.55, zorder=0)
    ax.axvline(0, color=SLATE, lw=1.1, zorder=1)
    for yi, effect, pval in zip(y, nes, values["nominal_p"]):
        ax.hlines(yi, 0, effect, color=BLUE, lw=3.2, zorder=2)
        ax.scatter(effect, yi, s=105, color=NAVY, edgecolor=WHITE, linewidth=1.2, zorder=3)
        ax.text(0.08, yi, f"NES {effect:.2f}   p {format_p(float(pval))}", va="center", ha="left", fontsize=9.2)

    ax.set_yticks(y)
    ax.set_yticklabels(values["model"])
    ax.set_xlim(-3.65, 1.05)
    ax.set_ylim(-0.8, 2.8)
    ax.set_xlabel("Normalized enrichment score (negative indicates lower pathway expression)")
    ax.grid(axis="x", color=LIGHT_GREY, linewidth=0.8)
    clean_axes(ax)
    ax.set_title("Proteasome assembly converges retrospectively across three model classes", loc="left", pad=14)
    ax.text(
        0.985,
        0.965,
        f"Signed meta-z = {row['meta_z']:.2f}\nmeta-FDR = {row['meta_FDR']:.4f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        color=NAVY,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": PALE_PURPLE, "edgecolor": PURPLE},
    )
    genes = str(strict["shared_leading_genes_all_3"]).replace(";", ", ")
    fig.text(0.11, 0.015, f"Shared leading edge: {genes}", fontsize=8.8, color=SLATE)
    save_figure(fig, "fig2_proteasome_cross_model", "Proteasome assembly cross-model convergence")
    return values


def figure_3_rotenone(v: dict[str, object]) -> pd.DataFrame:
    values = v["rotenone"].copy()
    label_map = {
        "D15_50nM_24h": "Day 15 | 50 nM | 24 h (frozen primary)",
        "D15_50nM_12h": "Day 15 | 50 nM | 12 h",
        "D15_100nM_24h": "Day 15 | 100 nM | 24 h",
        "D8_50nM_24h": "Day 8 | 50 nM | 24 h",
        "D8_50nM_12h": "Day 8 | 50 nM | 12 h",
        "D8_100nM_24h": "Day 8 | 100 nM | 24 h",
    }
    values["display_label"] = values["contrast"].map(label_map)
    values["analysis_role"] = np.where(
        values["contrast"] == "D15_50nM_24h", "Frozen primary", "Declared sensitivity"
    )
    values["nominal_p_lt_0_05"] = values["nominal_p"] < 0.05
    write_tsv(values, OUT / "fig3_rotenone_values.tsv")

    fig, ax = plt.subplots(figsize=(10.2, 5.9))
    y = np.arange(len(values))[::-1]
    colors = []
    edges = []
    for _, row in values.iterrows():
        if row["analysis_role"] == "Frozen primary":
            colors.append(GOLD)
            edges.append("#8A5B00")
        elif row["nominal_p_lt_0_05"]:
            colors.append(BLUE)
            edges.append(NAVY)
        else:
            colors.append(MID_GREY)
            edges.append(SLATE)

    ax.axvspan(-1.85, 0, color=PALE_BLUE, alpha=0.5, zorder=0)
    ax.axvline(0, color=SLATE, lw=1.0)
    bars = ax.barh(y, values["NES"], height=0.58, color=colors, edgecolor=edges, linewidth=1.0)
    for yi, (_, row), bar in zip(y, values.iterrows(), bars):
        ax.text(1.23, yi, f"p = {format_p(float(row['nominal_p']))}", ha="right", va="center", fontsize=9)
        if row["analysis_role"] == "Frozen primary":
            ax.text(float(row["NES"]) - 0.05, yi, "Primary", ha="right", va="center", fontsize=8.4, fontweight="bold", color="#6B4600")

    ax.set_yticks(y)
    ax.set_yticklabels(values["display_label"])
    ax.set_xlim(-1.85, 1.32)
    ax.set_xlabel("Proteasome assembly NES")
    ax.grid(axis="x", color=LIGHT_GREY, linewidth=0.8)
    clean_axes(ax)
    ax.set_title("Frozen mature-LUHMES rotenone replication did not meet its threshold", loc="left", pad=14)
    legend = [
        Patch(facecolor=GOLD, edgecolor="#8A5B00", label="Frozen primary"),
        Patch(facecolor=BLUE, edgecolor=NAVY, label="Sensitivity with nominal p < 0.05"),
        Patch(facecolor=MID_GREY, edgecolor=SLATE, label="Other sensitivity"),
    ]
    ax.legend(handles=legend, loc="lower left", bbox_to_anchor=(0, -0.25), ncol=3, frameon=False, fontsize=8.7)
    fig.text(0.11, 0.015, "Nominal p values are shown. Sensitivity contrasts do not relabel the frozen primary outcome.", fontsize=8.8, color=SLATE)
    save_figure(fig, "fig3_rotenone_replication", "Frozen GSE116280 rotenone replication")
    return values


def donor_panel(
    ax: mpl.axes.Axes,
    data: pd.DataFrame,
    title: str,
    genes: int,
) -> None:
    for group, x, color, label in [(0, 0, BLUE, "Control"), (1, 1, GOLD, "PD")]:
        subset = data.loc[data["pd"] == group, "pathway_score"].sort_values().to_numpy()
        jitter = np.linspace(-0.13, 0.13, len(subset)) if len(subset) > 1 else np.array([0.0])
        ax.scatter(
            np.full(len(subset), x) + jitter,
            subset,
            s=48,
            color=color,
            edgecolor=WHITE,
            linewidth=0.8,
            alpha=0.95,
            zorder=3,
        )
        mean = float(np.mean(subset))
        ax.plot([x - 0.22, x + 0.22], [mean, mean], color=TEXT, lw=2.1, zorder=4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Control\nn={int((data['pd'] == 0).sum())}", f"PD\nn={int((data['pd'] == 1).sum())}"])
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("Donor pathway score")
    ax.set_title(title, loc="left", pad=10)
    ax.text(0.98, 0.03, f"{genes} measured genes", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=SLATE)
    ax.grid(axis="y", color=LIGHT_GREY, linewidth=0.8)
    clean_axes(ax)


def figure_4_human(v: dict[str, object]) -> pd.DataFrame:
    primary = v["human_primary"]
    replication = v["human_replication"]
    estimates = pd.DataFrame(
        [
            {
                "cohort": "GSE178265 primary",
                "beta_PD": primary["beta_PD"],
                "ci95_low": primary["ci95_low"],
                "ci95_high": primary["ci95_high"],
                "p_one_sided_down": primary["p_one_sided_down"],
                "n_PD": int(primary["n_PD"]),
                "n_control": int(primary["n_control"]),
                "measured_genes_n": int(primary["measured_genes_n"]),
                "outcome": primary["outcome"],
            },
            {
                "cohort": "GSE243639 replication",
                "beta_PD": replication["beta_PD"],
                "ci95_low": replication["ci95_low"],
                "ci95_high": replication["ci95_high"],
                "p_one_sided_down": replication["p_one_sided_down"],
                "n_PD": int(replication["n_PD"]),
                "n_control": int(replication["n_control"]),
                "measured_genes_n": int(replication["measured_genes_n"]),
                "outcome": replication["outcome"],
            },
        ]
    )
    write_tsv(estimates, OUT / "fig4_human_model_values.tsv")

    fig = plt.figure(figsize=(11.7, 6.5))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.28], wspace=0.42)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[0, 2])

    donor_panel(ax1, v["primary_donors"], "A. GSE178265 primary", int(primary["measured_genes_n"]))
    donor_panel(ax2, v["replication_donors"], "B. GSE243639 replication", int(replication["measured_genes_n"]))

    y = np.array([1, 0])
    effects = estimates["beta_PD"].to_numpy(dtype=float)
    lows = estimates["ci95_low"].to_numpy(dtype=float)
    highs = estimates["ci95_high"].to_numpy(dtype=float)
    colors = [TEAL, SLATE]
    ax3.axvline(0, color=SLATE, lw=1.0, zorder=0)
    for yi, effect, low, high, color, (_, row) in zip(y, effects, lows, highs, colors, estimates.iterrows()):
        ax3.errorbar(
            effect,
            yi,
            xerr=[[effect - low], [high - effect]],
            fmt="o",
            ms=8,
            color=color,
            ecolor=color,
            elinewidth=2.0,
            capsize=4,
            zorder=3,
        )
        ax3.text(
            0.58,
            yi,
            f"β = {effect:+.3f}\none-sided p = {format_p(float(row['p_one_sided_down']))}",
            ha="right",
            va="center",
            fontsize=8.8,
        )
    ax3.set_yticks(y)
    ax3.set_yticklabels(["GSE178265\nprimary", "GSE243639\nreplication"])
    ax3.set_xlim(-1.02, 0.62)
    ax3.set_ylim(-0.65, 1.65)
    ax3.set_xlabel("Adjusted PD coefficient (95% CI)")
    ax3.set_title("C. Covariate-adjusted effects", loc="left", pad=10)
    ax3.grid(axis="x", color=LIGHT_GREY, linewidth=0.8)
    clean_axes(ax3)

    fig.suptitle("Human substantia nigra proteasome validation was primary-only", x=0.06, ha="left", y=1.01, fontsize=16, fontweight="bold")
    fig.text(
        0.06,
        0.015,
        "Dots are donors; horizontal bars are raw group means. Model estimates adjust for age, sex, and postmortem interval.",
        fontsize=8.8,
        color=SLATE,
    )
    save_figure(fig, "fig4_human_sn_validation", "Human substantia nigra proteasome validation")
    return estimates


def figure_5_matrix(v: dict[str, object]) -> pd.DataFrame:
    human_panel = v["human_panel"]
    continuous = v["continuous"]
    proteasome = v["proteasome"]
    primary_rotenone = v["primary_rotenone"]
    multitox_rot = v["multitox_rot"]
    heldout = v["heldout"]
    human_primary = v["human_primary"]
    human_replication = v["human_replication"]
    external = v["external"]

    matrix_rows = [
        ("21-gene network panel", "Retrospective nomination", "Nominated", "21 fixed genes", "results/original20_table1.csv"),
        ("21-gene network panel", "Additional model evidence", "Not confirmed", "0 FDR genes in GSE46798/GSE17542", "results/v2/limma_confirmation/README.md"),
        ("21-gene network panel", "Frozen primary test", "Not tested", "", ""),
        ("21-gene network panel", "Human disease evidence", "Not confirmed", f"{int(human_panel['meta_FDR_sig_n'])} meta-FDR genes", SOURCE_PATHS[0]),
        ("21-gene network panel", "Independent replication", "Not tested", "", ""),
        ("Continuous network score", "Retrospective nomination", "Nominated", "Minimum-distance proximity", "results/robust_sporadic_candidates_ranked.csv"),
        ("Continuous network score", "Additional model evidence", "Mixed/weak", f"max |partial r| = {continuous['partial_spearman_r'].abs().max():.3f}", SOURCE_PATHS[4]),
        ("Continuous network score", "Frozen primary test", "Not tested", "", ""),
        ("Continuous network score", "Human disease evidence", "Not confirmed", f"partial r = {continuous.loc[continuous['outcome']=='human_PD_postmortem_meta','partial_spearman_r'].iloc[0]:.3f}", SOURCE_PATHS[4]),
        ("Continuous network score", "Independent replication", "Not tested", "", ""),
        ("Proteasome assembly", "Retrospective nomination", "Retrospective signal", f"meta-FDR = {proteasome['meta_FDR']:.4f}", SOURCE_PATHS[5]),
        ("Proteasome assembly", "Additional model evidence", "Mixed/weak", f"GSE150005 rotenone p = {multitox_rot['PValue_one_sided_down']:.3f}; day-8 nominal support", SOURCE_PATHS[8]),
        ("Proteasome assembly", "Frozen primary test", "Not confirmed", f"GSE116280 p = {primary_rotenone['nominal_p']:.4f}", SOURCE_PATHS[7]),
        ("Proteasome assembly", "Human disease evidence", "Confirmed", f"GSE178265 p = {human_primary['p_one_sided_down']:.4f}", SOURCE_PATHS[13]),
        ("Proteasome assembly", "Independent replication", "Not confirmed", f"GSE243639 p = {human_replication['p_one_sided_down']:.4f}", SOURCE_PATHS[13]),
        ("Cilium / Hedgehog / trafficking", "Retrospective nomination", "Nominated", "GSE150005 pathway family", "results/v2/gse150005_multitoxicant/GSE150005_exploratory_mechanism_summary.tsv"),
        ("Cilium / Hedgehog / trafficking", "Additional model evidence", "Not tested", "", ""),
        ("Cilium / Hedgehog / trafficking", "Frozen primary test", "Not confirmed", f"GSE196190 p = {heldout['CAMERA_p_one_sided_down']:.4f}", SOURCE_PATHS[9]),
        ("Cilium / Hedgehog / trafficking", "Human disease evidence", "Not tested", "", ""),
        ("Cilium / Hedgehog / trafficking", "Independent replication", "Not confirmed", f"LUHMES p = {external['CAMERA_p_one_sided_down'].min():.3f} to {external['CAMERA_p_one_sided_down'].max():.3f}", SOURCE_PATHS[10]),
    ]
    values = pd.DataFrame(matrix_rows, columns=["claim", "evidence_stage", "status", "statistic", "source"])
    write_tsv(values, OUT / "fig5_validation_matrix.tsv")

    claims = [
        "21-gene\nnetwork panel",
        "Continuous\nnetwork score",
        "Proteasome\nassembly",
        "Cilium / Hedgehog /\ntrafficking",
    ]
    claim_keys = [
        "21-gene network panel",
        "Continuous network score",
        "Proteasome assembly",
        "Cilium / Hedgehog / trafficking",
    ]
    stages = [
        "Retrospective\nnomination",
        "Additional model\nevidence",
        "Frozen primary\ntest",
        "Human disease\nevidence",
        "Independent\nreplication",
    ]
    stage_keys = [
        "Retrospective nomination",
        "Additional model evidence",
        "Frozen primary test",
        "Human disease evidence",
        "Independent replication",
    ]
    status_style = {
        "Nominated": (PALE_PURPLE, PURPLE),
        "Retrospective signal": (PALE_PURPLE, PURPLE),
        "Confirmed": (PALE_TEAL, TEAL),
        "Mixed/weak": (PALE_GOLD, GOLD),
        "Not confirmed": (PALE_GREY, SLATE),
        "Not tested": (WHITE, LIGHT_GREY),
    }

    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    ax.set_xlim(0, len(stages))
    ax.set_ylim(0, len(claims))
    ax.invert_yaxis()
    ax.axis("off")

    for col, stage in enumerate(stage_keys):
        ax.text(col + 0.5, -0.18, stages[col], ha="center", va="bottom", fontsize=9.4, fontweight="bold", color=NAVY)
    for row, claim in enumerate(claim_keys):
        ax.text(-0.08, row + 0.5, claims[row], ha="right", va="center", fontsize=9.5, fontweight="bold", color=TEXT)
        for col, stage in enumerate(stage_keys):
            item = values.loc[(values["claim"] == claim) & (values["evidence_stage"] == stage)].iloc[0]
            fill, edge = status_style[item["status"]]
            rect = Rectangle((col + 0.03, row + 0.04), 0.94, 0.92, facecolor=fill, edgecolor=edge, linewidth=1.1)
            ax.add_patch(rect)
            status = str(item["status"])
            statistic = str(item["statistic"])
            if status == "Retrospective signal":
                status = "Retrospective\nsignal"
            elif status == "Not confirmed":
                status = "Not\nconfirmed"
            elif status == "Mixed/weak":
                status = "Mixed / weak"
            elif status == "Not tested":
                status = "Not tested"
            cell_text = status
            if statistic and statistic != "nan":
                compact = "\n".join(
                    textwrap.fill(
                        part,
                        width=26,
                        break_long_words=False,
                        break_on_hyphens=False,
                    )
                    for part in statistic.split("; ")
                )
                cell_text += f"\n{compact}"
            ax.text(col + 0.5, row + 0.5, cell_text, ha="center", va="center", fontsize=7.9, linespacing=1.2, color=TEXT)

    ax.set_title("Validation outcomes preserve the difference between nomination and confirmation", loc="left", pad=46)
    legend = [
        Patch(facecolor=PALE_PURPLE, edgecolor=PURPLE, label="Nominated / retrospective signal"),
        Patch(facecolor=PALE_TEAL, edgecolor=TEAL, label="Confirmed"),
        Patch(facecolor=PALE_GOLD, edgecolor=GOLD, label="Mixed / weak"),
        Patch(facecolor=PALE_GREY, edgecolor=SLATE, label="Not confirmed"),
        Patch(facecolor=WHITE, edgecolor=LIGHT_GREY, label="Not tested"),
    ]
    ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.18), ncol=5, frameon=False, fontsize=8.4)
    fig.text(0.105, 0.015, "Outcomes are assigned from the frozen plans and tracked result tables; null tests are retained.", fontsize=8.8, color=SLATE)
    save_figure(fig, "fig5_validation_outcome_matrix", "V2 validation outcome matrix")
    return values


def write_manifests() -> None:
    source_rows = []
    for relative in SOURCE_PATHS:
        path = ROOT / relative
        source_rows.append(
            {"source": relative, "sha256": sha256(path), "bytes": path.stat().st_size}
        )
    write_tsv(pd.DataFrame(source_rows), OUT / "source_manifest.tsv")

    output_rows = []
    for path in sorted(OUT.iterdir()):
        if path.name in {"figure_manifest.tsv", "source_manifest.tsv"} or path.suffix not in {".png", ".pdf", ".tsv"}:
            continue
        row = {
            "file": path.name,
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "width_px": "",
            "height_px": "",
        }
        if path.suffix == ".png":
            with Image.open(path) as image:
                row["width_px"], row["height_px"] = image.size
        output_rows.append(row)
    write_tsv(pd.DataFrame(output_rows), OUT / "figure_manifest.tsv")


def main() -> None:
    configure_style()
    values = release_values()
    figure_1_evidence_design(values)
    figure_2_cross_model(values)
    figure_3_rotenone(values)
    figure_4_human(values)
    figure_5_matrix(values)
    write_manifests()
    print(f"Saved AAN V2 figure package to {OUT}")


if __name__ == "__main__":
    main()
