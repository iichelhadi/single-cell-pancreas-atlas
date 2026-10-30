Order to run:

Step			Run
Inputs			Figure5H_prepare_pySCENIC_inputs.R
GRNBoost2, 50 seeds	run_grnboost2_multi_v3.py --n-runs 50
Consensus networks	build_consensus_pyscenic_v2.R
ctx + aucell		run_pyscenic_ctx_aucell_v2.pbs
Threshold comparison	compare_thresholds_post_ctx.py
Lineage groups		Figure5H_define_C3_lineage_groups.R
Statistics		Figure5H_SCENIC_consensus_analysis.R
Figures			Figure5H_SCENIC_figures.R
