library(tidyverse)

cc_df = read_csv("../shared_data/Combined_Eds1_Rgt1_Lys14_CCTargets.csv")

# per https://doi.org/10.1016/S0079-6603(03)01008-0
# pdr1 is the profile similarity gene from the GI map
pdr_targets = c(
  'RPN4',
  'RTA1',
  'YOR1',
  'HXT9',
  'HXT11',
  'PDR5',
  'PDR5',
  'PDR10',
  'PDR10',
  'SNQ2',
  'DPM1',
  'LCB2',
  'SUR2',
  'HXT12',
  'TPO1',
  'TPO4'
)

pdr_targets_id = c(
  'YDL020C',
  'YGR213C',
  'YGR281W',
  'YJL219W',
  'YOL156W',
  'YOR153W',
  'YOR153W',
  'YOR328W',
  'YOR328W',
  'YDR011W',
  'YPR183W',
  'YDR062W',
  'YDR297W',
  'YIL170W',
  'YLL028W',
  'YOR273C'
)

# there is no evidence for binding of these genes

View(filter(cc_df, id %in% pdr_targets_id ))
