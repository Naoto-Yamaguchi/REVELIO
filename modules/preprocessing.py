import pandas as pd
import numpy as np

def n_genes_filter(df, n_genes):
    print("Filtering out cells which have less than {} expressing genes ...".format(n_genes))
    print("Input df shape is {}".format(df.shape))
    cell_few_genes = []
    for column in df.columns:
        if len(df[column].value_counts()) < n_genes:
            cell_few_genes.append(column)
    df = df.drop(cell_few_genes, axis=1)
    print("Filtered df shape is {}".format(df.shape))
    return df

def n_cells_filter(df, n_cells):
    print("Filtering out genes which are expressed less than {} cells ...".format(n_cells))
    print("Input df shape is {}".format(df.shape))
    gene_few_cells = []
    for gene in df.T.columns:
        if len(df.T[gene].value_counts()) < n_cells:
            gene_few_cells.append(gene)
    df = df.drop(gene_few_cells, axis=0)
    print("Filtered df shape is {}".format(df.shape))
    return df

def log_transformation(df):
    print("Log transforming ...")
    log_transformation = lambda x: np.log(x+1)
    return df.apply(log_transformation)

def scaling(df):
    return df
