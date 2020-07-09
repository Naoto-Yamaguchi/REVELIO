import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def scatter(df, pc1, pc2, filename):
    # phases = list(set(df_preprocessed_pca_phase["phase"]))
    # colors 

    # for i in df_preprocessed_pca_phase.loc[[pc1, pc2, "phase"], :].iteritems():
    #     plt.scatter(i[1][pc1], i[1][pc2], c=colors[phases.index(i[1]["phase"])], s=15)
    # plt.savefig(filename)
    df = df.loc[[pc1, pc2, "phase"], :].T
    phases = df["phase"]

    fig = plt.figure()
    for phase in phases.unique():
        plt.scatter(df.loc[df["phase"] == phase, pc1], df.loc[df["phase"] == phase, pc2], label=phase)
    plt.title("PC1 - PC2")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    #plt.show()
    plt.savefig(filename)
