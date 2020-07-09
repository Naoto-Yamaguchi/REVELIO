import math
import pandas as pd
import numpy as np

# Variable_genes detection in REVELIO

# input: df (n_genes * n_cells)
# output: { gene1: ave_expr1, gene2: ave_expr2, ... }
def ave_expr(df_log_transformed):
    ave_exprs = {}
    for gene in df_log_transformed.iterrows():
        expr_sum = 0
        for count in gene[1]:
            #print(cell)
            expr_sum += np.exp(count)
        ave_expr = np.log(expr_sum / len(gene[1]))
        ave_exprs[gene[0]] = ave_expr

    return ave_exprs

def var_expr(df_count):
    aves = {}
    for gene in df_count.iterrows():
        expr_sum = 0
        for count in gene[1]:
            expr_sum += count
        ave_expr = expr_sum / len(gene[1])
        aves[gene[0]] = ave_expr

    var_exprs = {}
    for gene in df_count.iterrows():
        mse_sum = 0
        for count in gene[1]:
            mse_sum += (count - aves[gene[0]]) ** 2
        mse_ave = mse_sum / len(gene[1])
        var_exprs[gene[0]] = np.log(mse_ave / aves[gene[0]])

    return var_exprs


def variable_genes(df_count, df_log_transformed, n_bucket=20):
    print("Finding variable genes ...")
    ave_exprs = ave_expr(df_log_transformed)
    var_exprs = var_expr(df_count)

    # Normalizing vars
    min_expr = min(ave_exprs.values())
    max_expr = max(ave_exprs.values())
    step_size = (max_expr - min_expr) / n_bucket

    bucket2gene = {}
    gene2bucket = {}
    for i in range(1,n_bucket+1):
        bucket2gene[i] = []
    for gene in ave_exprs:
        for i in range(1,n_bucket+1): 
            if min_expr + (i-1) * step_size <= ave_exprs[gene] and ave_exprs[gene] <= min_expr + i * step_size:
                bucket2gene[i].append(gene)
                gene2bucket[gene] = i

    normalized_ds = {}
    for gene in var_exprs:
        d = var_exprs[gene]
        bucket_num = gene2bucket[gene]
        genes_in_the_bucket = bucket2gene[bucket_num]
        
        sm = 0
        for gene_in_the_bucket in genes_in_the_bucket:
            sm += var_exprs[gene_in_the_bucket]
        sm = sm / len(genes_in_the_bucket)
        
        mse_sm  = 0
        for gene_in_the_bucket in genes_in_the_bucket:
            mse_sm += (var_exprs[gene_in_the_bucket] - sm) ** 2
        mse_sm = mse_sm / len(genes_in_the_bucket)
        
        normalized_d = (d-sm) / math.sqrt(mse_sm)
        normalized_ds[gene] = normalized_d

    variable_genes = []
    for gene in ave_exprs:
        if 0.2 < ave_exprs[gene] and ave_exprs[gene] < 4 and 0.5 < normalized_ds[gene] and normalized_ds[gene] < 10:
            variable_genes.append(gene)

    print("{} variable genes are found".format(len(variable_genes)))
    return variable_genes
