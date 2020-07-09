import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# cb = Cyclebase




# 今はGene name, peaktimeを返す
def cb_human_cc_markers():

    # Cyclebase human_peridoci with peaktime 手動データ
    cb_human = pd.read_csv("cc_markers_data/cyclebase_human.tsv", sep="\t")
    cb_human = cb_human[~cb_human["Peaktime"].str.contains("non-periodic")]
    #cb_human = cb_human[["Matched name", "Peaktime"]].set_index("Matched name")
    #cb_human = cb_human.rename(columns={"Matched name" : "Gene name"})
    return cb_human

def human_to_mouse():

    # humanとmouseのGeneID対応 from Ensembl Biomart
    # G1:59, G1/S:37 ,S:36, G2:111, G2/M:39, M:96, total: 378
    human_to_mouse = pd.read_csv("cc_markers_data/mart_export.txt", sep="\t")
    human_to_mouse = human_to_mouse[["Protein stable ID", "Mouse gene stable ID", "Mouse gene name"]] # 必要があれば変える
    return human_to_mouse

def cb_mouse_cc_markers():
    # mouse cell cycle marker with peaktime
    # G1:59, G1/S:39, S:30, G2:100, G2/M:40, M:89, total: 357
    cb_human = cb_human_cc_markers()
    human2mouse = human_to_mouse()
    cb_mouse = pd.merge(cb_human, human2mouse, left_on="Identifier", right_on="Protein stable ID")
    cb_mouse = cb_mouse[["Mouse gene stable ID", "Mouse gene name", "Peaktime", "rank", "Phenotypes"]].drop_duplicates("Mouse gene stable ID")
    cb_mouse = cb_mouse[["Mouse gene name", "Peaktime"]].set_index("Mouse gene name")
    return cb_mouse



