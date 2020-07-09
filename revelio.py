import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from modules import cell_cycle_marker
from modules import preprocessing
from modules import vg
from modules import phase_assignment
from modules import rotation
from modules import visualize

parser = argparse.ArgumentParser(description="REVELIOの前半。渡されたdfを前処理し、PCAをおこない各PCのcell cycle scoreを算出")
parser.add_argument('filename', help="遺伝子発現行列: 行=遺伝子名, 列=細胞名")
parser.add_argument('organism')
parser.add_argument('--n_genes_filter', help="その遺伝子が発現している細胞数が一定数以上の遺伝子のみを残す" , type=int, default=0)
parser.add_argument('--n_cells_filter', help="その細胞で発現している遺伝子数が一定数以上の細胞のみを残す", type=int, default=0)
parser.add_argument('--log_transformation', help='対数変換をおこなうか', action='store_true')
parser.add_argument('--variable_genes', help="発現変動遺伝子の抽出をおこなうか", action='store_true')
parser.add_argument('--cell_cycle_assignment_file', help="既知の細胞周期割り当てファイル")
parser.add_argument('--n_components', type=int, default=50)

args = parser.parse_args()

os.makedirs("output", exist_ok=True)

# 汎用性高く書き直す
# inputは、gene * cellsの発現データ
# output1は、PC cell cycle score
# output2は、after rotationのPC cell cycle scoreとDC1, DC2からの色付き描画
# log変換の有無
# 発現変動遺伝子の出し方（qPCRでは同じ基準は使えない）
# phaseの割当方法の違い。マーカーか既にあるラベリングか
# 複数回のRotation。
# あとはグラフの描画
# まあとりあえず一から見直してわかりやすく、再利用性を高めて書く。


# class REVELIO:


#     def __init__(self, file):
#         self.df = 
    # input: DataFrame (genes * cells)

# def cc_markers_filter(df, df_cc_markers):
#     # input
#     # df: genes * cells (n_genes * n_cells)
#     # df_cc: genes * peaktime (n_genes * 1)
#     # indexでの結合。df, df_cc共に、indexにGene nameもしくは、Gene IDを用意する必要がある。

#     # output
#     # df: cell cycle genes * cells (n_cc_genes * n_cells)
#     return pd.merge(df_cc_markers, df, left_index=True, right_index=True)


# def phase_assign(df, df_cc_markers, df_filtered_log):
#     # df:
#     # df_cc_markers: n_gene * 1(peaktime)
#     df_cc_mean_corr = ave_cc_marker_expr_for_each_cell(df, df_cc_markers, df_filtered_log)
#     max_phase = []
#     for i in df_cc_mean_corr.iteritems():
#         expr = i[1][1:]
#         max_phase.append(expr.index[np.argmax(i[1][1:])])
#     df_max_phase = pd.DataFrame(max_phase).T
#     df_max_phase.columns = df_cc_mean_corr.columns
#     df_max_phase.index = ["phase"]
#     return df_max_phase # output: phase * cells (1 * n_cells)

# def ave_cc_marker_expr_for_each_cell(df, df_cc_markers, df_filtered_log, threshold=0.2):
#     phases = list(set(df["Peaktime"]))
#     df_cc_mean = df.groupby("Peaktime").mean().reindex(phases)
#     cc_markers_corr_list = corr_cc_markers(df, df_cc_mean, threshold)
#     df_cc_markers_corr = df_cc_markers[df_cc_markers.index.isin(cc_markers_corr_list)]
#     df_cc_mean_corr = pd.merge(df_cc_markers_corr, df_filtered_log, left_index=True, right_index=True).groupby("Peaktime").mean().reindex(phases)
#     return df_cc_mean_corr # output: phases * cells (n_phases * n_cells)
    
# def corr_cc_markers(df, df_cc_mean, threshold=0.2):
#     # less correlation ccMarkersを除く。
#     less_corr_markers = []
#     corr_markers = []
#     for i in df.iterrows():
#         peaktime = i[1][0]
#         expr = i[1][1:].astype(float)
#         mean_expr = df_cc_mean[df_cc_mean.index==peaktime].loc[peaktime]
#         corr = expr.corr(mean_expr)
#         if corr <= threshold:
#             less_corr_markers.append(i[0])
#         else:
#             corr_markers.append(i[0])
#     return corr_markers


# def cb_six_phases():
#     phases = ["G1", "G1/S", "S", "G2", "G2/M", "M"]
#     colors = ["red", "orange", "green", "lightblue", "blue", "purple"]
#     return phases, colors

# def three_phases():
#     return ["G1", "S", "G2M"]


# pc1-pc2のscatter plotをphaseラベル付きで描画
def plot_scatter_with_cc_phase(df_preprocessed_pca_phase, pc1, pc2):
    phases, colors = cb_six_phases()
    for i in df_preprocessed_pca_phase.loc[[pc1, pc2, "phase"], :].iteritems():
        plt.scatter(i[1][0], i[1][1], c=colors[phase.index(i[1]["phase"])], s=15)

    plt.title("scatter plot of PC{} and PC{}".format(pc1, pc2))
    plt.xlabel("PC{}".format(pc1))
    plt.ylabel("PC{}".format(pc2))
    plt.savefig("hoge")

    




# =======================================
# input: gene name + phase * cells ((n_genes + 1) * n_cells)

# n_genes, n_cells
def main(filename, organism, n_genes_filter, n_cells_filter, is_log_transformation, is_variable_genes, cell_cycle_assignment_file, n_components):
    # n_genes * n_cellsでindex=Gene name, columns=Cell nameのtsvファイルを要求
    df = pd.read_csv(filename, sep="\t", index_col=0)
    
    # Filtering
    df_cell_filtered = preprocessing.n_genes_filter(df, n_genes_filter)
    df_cell_gene_filtered = preprocessing.n_cells_filter(df_cell_filtered, n_cells_filter)
    df_preprocessed = df_cell_gene_filtered
    # Log Transformation
    if is_log_transformation: 
        df_filtered_log = preprocessing.log_transformation(df_cell_gene_filtered)
        df_preprocessed = df_filtered_log

    
    # Cell Cycle Assignment
    if cell_cycle_assignment_file:
        # df_phase: 1 * n_cells (index=["phase"], columns=["cell1", "cell2", ...])
        df_phase = pd.read_csv(cell_cycle_assignment_file, sep="\t", index_col=0)
        df_phase.index = ["phase"]
    else:
        # get cc markers
        if organism == "human":
            print("Getting {} cell cycle markers ...".format(organism))
            cb_human = cell_cycle_marker.cb_human_cc_markers()
            df_cc_markers = cb_human[["Matched name", "Peaktime"]].set_index("Matched name")
        elif organism == "mouse":
            print("Getting {} cell cycle markers ...".format(organism))
            df_cc_markers = cell_cycle_marker.cb_mouse_cc_markers()
        else:
            print("No Cyclebase data for the {}".format(organism))
            return 0
        print("The number of cell cycle marker is ... {}".format(df_cc_markers.shape[0]))
        
        # phase assign by cc markers
        print("Cell cycle phase assigning...")
        df_phase = phase_assignment.max_phase_assignment(df_preprocessed, df_cc_markers)


    # Variable Genes (ここでの発現変動遺伝子同定方法はlog変換が前提)
    if is_log_transformation and is_variable_genes:
        variable_genes = vg.variable_genes(df_cell_gene_filtered, df_filtered_log, 20)
        df_filtered_log_vg = df_filtered_log[df_filtered_log.index.isin(variable_genes)]
        df_preprocessed = df_filtered_log_vg

    # PCA
    print("Doing PCA...")
    pca = PCA(n_components=n_components)
    pca.fit(df_preprocessed.T)
    # df_preprocessed_pca: n_components * n_cells
    df_preprocessed_pca = pd.DataFrame(pca.transform(df_preprocessed.T).T)
    df_preprocessed_pca.columns = df_preprocessed.columns


    # PC cell cycle score
    # df_preprocessed_pca_phase_mean: n_phase * n_components
    df_preprocessed_pca_phase = pd.concat([df_preprocessed_pca, df_phase])
    df_preprocessed_pca_phase_mean = df_preprocessed_pca_phase.T.set_index("phase").astype(float).groupby("phase").mean()
    df_preprocessed_pca_phase_var = df_preprocessed_pca_phase_mean.var()
    fig = plt.figure()
    plt.scatter(np.arange(n_components), df_preprocessed_pca_phase_var)
    plt.title("PC cell cycle score")
    plt.savefig("output/PC_cc_score_before_rotation.png".format(filename))
    print("Figure: PC cell cycle score before rotation is saved!")

    # Scatter plot of PCs
    cc_scores = np.array(df_preprocessed_pca_phase_var).argsort()[::-1]
    pc0 = cc_scores[0]
    pc1 = cc_scores[1]
    pc2 = cc_scores[2]
    visualize.scatter(df_preprocessed_pca_phase, pc0, pc1, "output/scatter_pc{}_pc{}_before_rotation.png".format(pc0, pc1))

    # 上の結果をみて、pcの値を決める
    # Rotation
    df_preprocessed_pca_rotated = rotation.rotation(df_preprocessed_pca, df_phase, pc0, pc1, pc2, n_components)

    # PC cell cycle score
    # df_preprocessed_pca_rotated_phase_mean: n_phase * n_components
    df_preprocessed_pca_rotated_phase = pd.concat([df_preprocessed_pca_rotated, df_phase])
    df_preprocessed_pca_rotated_phase_mean = df_preprocessed_pca_rotated_phase.T.set_index("phase").astype(float).groupby("phase").mean()
    df_preprocessed_pca_rotated_phase_var = df_preprocessed_pca_rotated_phase_mean.var()
    fig = plt.figure()
    plt.scatter(np.arange(n_components), df_preprocessed_pca_rotated_phase_var)
    plt.title("PC cell cycle score")
    plt.savefig("output/PC_cc_score_after_one_rotation.png".format(filename))
    print("Figure: PC cell cycle score after one rotation is saved!")

    # Scatter plot of PCs
    visualize.scatter(df_preprocessed_pca_rotated_phase, pc0, pc1, "output/scatter_pc{}_pc{}_after_rotation.png".format(pc0, pc1))
    # More rotation


if __name__ == "__main__":
    #args = sys.argv
    filename = args.filename
    organism = args.organism
    n_genes_filter = args.n_genes_filter
    n_cells_filter = args.n_cells_filter
    is_log_transformation = args.log_transformation
    is_variable_genes = args.variable_genes
    cell_cycle_assignment_file = args.cell_cycle_assignment_file
    n_components = args.n_components
    main(filename, organism, n_genes_filter, n_cells_filter, is_log_transformation, is_variable_genes, cell_cycle_assignment_file, n_components)


# 問題PCAの再現性
# phaseの違いによるPC cell cycleの違い
