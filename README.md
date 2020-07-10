# REVELIO
論文: D.Schwabe et al., The transcriptome dynamics of single cells during the cell cycle, bioRxiv, 2019
url: https://www.biorxiv.org/content/10.1101/2019.12.23.887570v1.full#disp-formula-17
まとめスライド: https://docs.google.com/presentation/d/1F0MEIoLcgEodNyPm6SZAXRDj5GVAP2YVCRqw2RMkOr0/edit#slide=id.g88d8237a0e_2_59


# フォルダ構成
メイン（一部省略）

├── cc_markers_data
│   ├── cyclebase_human.tsv
│   ├── human_periodic.tsv
│   ├── mart_export.txt
│   └── not_used
│       ├── 
├── expr_data
│   ├── h9.tsv
│   ├── h9_phase.tsv
│   ├── mESC_phase.tsv
│   └── mESC_tpm.tsv
├── modules
│   ├── cell_cycle_marker.py
│   ├── phase_assignment.py
│   ├── preprocessing.py
│   ├── rotation.py
│   ├── vg.py
│   ├── visualize.py
│   ├── __pycache__
│   │   ├──
├── output
├── revelio.py


# 実行
`python3 revelio.py <遺伝子発現行列> <生物種> --n_genes_filter <フィルタリングする遺伝子数> --n_cells_filter <フィルタリングする細胞数> --log_transformation (対数変換する場合) --variable_genes (発現変動遺伝子抽出する場合) --cell_cycle_assignment <細胞周期割り当てファイル> (細胞周期割り当てが予めわかっている場合)`

ex1)
マウスのmESCのtpmデータ。フィルタリングあり、対数変換あり、発現変動遺伝子抽出なし、細胞周期割り当ては、Cyclebaseに要録されている細胞周期マーカーの遺伝子発現量を利用。
`python3 revelio.py expr_data/mESC_tpm.tsv mouse --n_genes_filter 500 --n_cells_filter 5 --log_transformation`

ex2)
マウスのmESCのtpmデータ。フィルタリングあり、対数変換あり、発現変動遺伝子抽出あり、細胞周期割り当ては、Hoechst染色で得られた既知の情報を利用。
`python3 revelio.py expr_data/mESC_tpm.tsv mouse --n_genes_filter 500 --n_cells_filter 5 --log_transformation --variable_genes --cell_cycle_assignment expr_data/mESC_phase.tsv`

ex3)
ヒトの培養細胞H9のqPCEデータ。フィルタリングなし、対数変換なし、発現変動遺伝子抽出なし、細胞周期割り当ては、実験で得られた既知の情報を利用。
`python3 revelio.py expr_data/hg.tsv human --cell_cycle_assignment expr_data/h9_phase.tsv`

# 備考
複数回Rotationする場合は、rotation.rotation.pyを呼び出せば簡単な呼び出し方法は未実装。
