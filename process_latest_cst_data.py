
import proc_cst

# # transfer from original data format into simplfied tsv
# ########################################################################
# proc_cst.save_lung_data_to_tsv('phospho')
# proc_cst.save_lung_data_to_tsv('AcK')
# proc_cst.save_lung_data_to_tsv('Kme1')
# proc_cst.save_lung_data_to_tsv('Rme1')

# # calculate ratios of cell lines to normal tissue and save to new files
# ########################################################################
# proc_cst.calc_ratios('phospho')
# proc_cst.calc_ratios('AcK')
# proc_cst.calc_ratios('Kme1')
# proc_cst.calc_ratios('Rme1')

# combine ratios from all plexes together
proc_cst.combine_ratios('phospho')
proc_cst.combine_ratios('Ack')
proc_cst.combine_ratios('Kme1')
proc_cst.combine_ratios('Rme1')


