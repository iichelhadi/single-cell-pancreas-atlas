#!/usr/bin/env python3

import itertools
import os

import loompy
import numpy as np
import pandas as pd


FILES = {
    "50": "revision_out/scenic_consensus50.loom",
    "70": "revision_out/scenic_consensus70.loom",
    "80": "revision_out/scenic_consensus80.loom",
}

OUT_DIR = "revision_out"


def load_auc(path):
    with loompy.connect(path, mode="r", validate=False) as ds:
        auc = ds.ca.RegulonsAUC
        names = list(auc.dtype.names)

        matrix = np.column_stack([auc[name] for name in names])

        if "CellID" in ds.ca:
            cells = np.asarray(ds.ca.CellID)
        else:
            cells = np.asarray([f"Cell_{i+1}" for i in range(ds.shape[1])])

    return pd.DataFrame(matrix, index=cells, columns=names)


auc = {threshold: load_auc(path) for threshold, path in FILES.items()}

summary = pd.DataFrame({
    "threshold_pct": [50, 70, 80],
    "n_regulons": [
        auc["50"].shape[1],
        auc["70"].shape[1],
        auc["80"].shape[1],
    ],
})

summary.to_csv(
    os.path.join(OUT_DIR, "post_ctx_regulon_counts.tsv"),
    sep="\t",
    index=False,
)

pairwise_rows = []
correlation_rows = []

for a, b in itertools.combinations(["50", "70", "80"], 2):
    regs_a = set(auc[a].columns)
    regs_b = set(auc[b].columns)

    shared = sorted(regs_a & regs_b)
    union = regs_a | regs_b

    pairwise_rows.append({
        "threshold_A": a,
        "threshold_B": b,
        "n_regulons_A": len(regs_a),
        "n_regulons_B": len(regs_b),
        "n_shared_regulons": len(shared),
        "regulon_name_jaccard": len(shared) / len(union),
    })

    for regulon in shared:
        x = auc[a][regulon]
        y = auc[b][regulon]

        correlation_rows.append({
            "threshold_A": a,
            "threshold_B": b,
            "regulon": regulon,
            "pearson_r": x.corr(y, method="pearson"),
            "spearman_r": x.corr(y, method="spearman"),
        })

pairwise = pd.DataFrame(pairwise_rows)
correlations = pd.DataFrame(correlation_rows)

pairwise.to_csv(
    os.path.join(OUT_DIR, "post_ctx_regulon_overlap.tsv"),
    sep="\t",
    index=False,
)

correlations.to_csv(
    os.path.join(OUT_DIR, "post_ctx_shared_regulon_auc_correlations.tsv"),
    sep="\t",
    index=False,
)

correlation_summary = (
    correlations
    .groupby(["threshold_A", "threshold_B"])
    .agg(
        n_shared_regulons=("regulon", "size"),
        median_pearson_r=("pearson_r", "median"),
        median_spearman_r=("spearman_r", "median"),
        min_spearman_r=("spearman_r", "min"),
    )
    .reset_index()
)

correlation_summary.to_csv(
    os.path.join(OUT_DIR, "post_ctx_auc_correlation_summary.tsv"),
    sep="\t",
    index=False,
)

print("\nRegulon counts:")
print(summary.to_string(index=False))

print("\nRegulon-name overlap:")
print(pairwise.to_string(index=False))

print("\nAUCell correlation summary:")
print(correlation_summary.to_string(index=False))
