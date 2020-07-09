import pandas as pd
import numpy as np


def max_phase_assignment(df_preprocessed, df_cc_markers):
    # df_preprocessed: n_genes * n_cells (index=["gene1", "gene2", ...], columns=["cell1", "cell2", ...])
    # df_cc_markers: n_cc_genes * 1 (index=["cc_gene1", "cc_gene2", ...], columns=["peaktime"])

    # df_preprocessed_cc: n_cc_genes * n_cells (index=["cc_gene1", "cc_gene2", ...], columns=["cell1", "cell2", ...])
    df_preprocessed_cc = pd.merge(df_cc_markers, df_preprocessed, left_index=True, right_index=True)

    # df_preprocessed_cc_filtered: n_filtered_cc_genes * n_cells (index=["filtered_cc_gene1", "filtered_cc_gene2", ...], columns=["cell1", "cell2", ...])
    df_preprocessed_cc_filtered = cc_markers_correlation_filter(df_preprocessed_cc)

    # df_preprocessed_cc_filtered_mean: n_phases * n_cells (index=["S", "G1", "G2/M", ..., columns=["cell1", "cell2", ...]])
    #phases = list(set(df_preprocessed_cc_filtered["Peaktime"]))
    #phases = ["G1", "G1/S", "S", "G2", "G2/M", "M"]
    #phases = ['S', 'G2', 'G2/M', 'G1', 'M', 'G1/S']
    #print("phases after filtering: {}".format(phases))
    df_preprocessed_cc_filtered_mean = df_preprocessed_cc_filtered.groupby("Peaktime").mean()#.reindex(phases)

    # max_phase_assignment
    max_phase = []
    for i in df_preprocessed_cc_filtered_mean.iteritems():
        expr = i[1][1:]
        max_phase.append(expr.index[np.argmax(i[1][1:])])
    df_max_phase = pd.DataFrame(max_phase).T
    df_max_phase.columns = df_preprocessed_cc_filtered_mean.columns
    df_max_phase.index = ["phase"]
    # df_max_phase: 1 * n_cells (index=["phase"], columns=["cell1", "cell2", ...])
    return df_max_phase

def cc_markers_correlation_filter(df_preprocessed_cc):
    print("Filtering out cell cycle marker genes with low correlation...")
    less_corr_markers = []
    corr_markers = []
    #phases = list(set(df_preprocessed_cc["Peaktime"]))
    #phases = ["G1", "G1/S", "S", "G2", "G2/M", "M"]
    #phases = ['S', 'G2', 'G2/M', 'G1', 'M', 'G1/S']
    #print("phases before filtering: {}".format(phases))
    df_preprocessed_cc_mean = df_preprocessed_cc.groupby("Peaktime").mean()#.reindex(phases)
    for i in df_preprocessed_cc.iterrows():
        peaktime = i[1][0]
        expr = i[1][1:].astype(float)
        mean_expr = df_preprocessed_cc_mean.loc[peaktime, :]
        corr = expr.corr(mean_expr)
        if corr <= 0.2:
            less_corr_markers.append(i[0])
        else:
            corr_markers.append(i[0])
    df_preprocessed_cc_filtered = df_preprocessed_cc[df_preprocessed_cc.index.isin(corr_markers)]
    print("Filtered df shape is {}".format(df_preprocessed_cc_filtered.shape))
    return df_preprocessed_cc_filtered
