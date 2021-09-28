import pandas as pd
import numpy as np

base_path = 'code/data/Fig2/'

# Load the RM absence-prescence from the table data
data_RM = pd.read_table(base_path+'NCBI_data/pairs_out.txt', index_col=0).drop(columns=['count_pairs'])

# Store the motif list
pd.DataFrame(data_RM.columns, columns=['R-M motif']).to_csv(base_path + 'motifs.csv', index=False)

# Load the taxonomy table
data_taxonomy = pd.read_table(base_path+'NCBI_data/taxa_out.txt', index_col=0)


# Filter out entries with genus in square brackets
data_taxonomy = data_taxonomy[~data_taxonomy.genus.str.contains('\[')]
data_taxonomy = data_taxonomy[~data_taxonomy.genus.str.contains('\]')]

# Filter out entries with the genus "Candidatus"
data_taxonomy = data_taxonomy[~data_taxonomy.genus.str.contains('candidatus', case=False)]

# Filter out the entries 'sp._SM30' and  'sp._JZB09'
data_taxonomy = data_taxonomy[~data_taxonomy.genus.str.contains('sp._SM30', case=False)]
data_taxonomy = data_taxonomy[~data_taxonomy.genus.str.contains('sp._JZB09', case=False)]

# Filter out the genomes whose genus are not clear
data_taxonomy = data_taxonomy.drop(index=['GCF_003288255.1_ASM328825v1',
                                          'GCF_000012225.1_ASM1222v1',
                                          'GCF_000023725.1_ASM2372v1',
                                          'GCF_000198515.1_ASM19851v1',
                                          'GCF_000242935.2_ASM24293v3',
                                          'GCF_000287335.1_ASM28733v1',
                                          'GCF_000287355.1_ASM28735v1',
                                          'GCF_000299095.1_ASM29909v1',
                                          'GCF_000299115.1_ASM29911v1',
                                          'GCF_000304455.1_CCh_cEper1',
                                          'GCF_000330845.1_ASM33084v1',
                                          'GCF_000342265.1_ASM34226v1',
                                          'GCF_000473305.1_ASM47330v1',
                                          'GCF_000597885.1_ASM59788v1',
                                          'GCF_000730245.1_ASM73024v1',
                                          'GCF_000746585.1_ASM74658v2',
                                          'GCF_000801295.1_ASM80129v1',
                                          'GCF_000814825.1_ASM81482v1',
                                          'GCF_000828975.1_ASM82897v1',
                                          'GCF_000829235.1_ASM82923v1',
                                          'GCF_000972765.1_ASM97276v1',
                                          'GCF_001023575.1_ASM102357v1',
                                          'GCF_001518995.2_ASM151899v2',
                                          'GCF_001547755.1_ASM154775v1',
                                          'GCF_001584725.1_ASM158472v1',
                                          'GCF_001679665.1_ASM167966v1',
                                          'GCF_001688705.2_ASM168870v2',
                                          'GCF_001688905.2_ASM168890v2',
                                          'GCF_001688965.2_ASM168896v2',
                                          'GCF_001712875.1_ASM171287v1',
                                          'GCF_001717545.1_ASM171754v1',
                                          'GCF_001880285.1_ASM188028v1',
                                          'GCF_001953215.1_ASM195321v1',
                                          'GCF_002009425.1_ASM200942v1',
                                          'GCF_002119605.1_ASM211960v1',
                                          'GCF_002151545.1_ASM215154v1',
                                          'GCF_002284855.1_ASM228485v1',
                                          'GCF_002284875.1_ASM228487v1',
                                          'GCF_002284895.1_ASM228489v1',
                                          'GCF_002284915.1_ASM228491v1',
                                          'GCF_002356635.1_ASM235663v1',
                                          'GCF_002838185.1_ASM283818v1',
                                          'GCF_002838205.1_ASM283820v1',
                                          'GCF_002850435.1_ASM285043v1',
                                          'GCF_002863805.1_ASM286380v1',
                                          'GCF_002892535.1_ASM289253v1',
                                          'GCF_002903045.1_ASM290304v1',
                                          'GCF_002903065.1_ASM290306v1',
                                          'GCF_002903085.1_ASM290308v1',
                                          'GCF_002998355.1_ASM299835v1',
                                          'GCF_002998695.1_ASM299869v1',
                                          'GCF_002999035.1_ASM299903v1',
                                          'GCF_003008555.1_ASM300855v1',
                                          'GCF_003063705.1_ASM306370v1',
                                          'GCF_003072645.1_ASM307264v1',
                                          'GCF_003151175.1_ASM315117v1',
                                          'GCF_003285105.1_ASM328510v1',
                                          'GCF_003330885.1_ASM333088v1',
                                          'GCF_003351745.1_ASM335174v1',
                                          'GCF_003351905.1_ASM335190v1',
                                          'GCF_003443635.1_ASM344363v1',
                                          'GCF_003544895.1_ASM354489v1',
                                          'GCF_003573975.1_ASM357397v1',
                                          'GCF_003573995.1_ASM357399v1',
                                          'GCF_003574135.1_ASM357413v1',
                                          'GCF_003584665.1_ASM358466v1',
                                          'GCF_003584705.1_ASM358470v1',
                                          'GCF_003932735.1_ASM393273v1',
                                          'GCF_003966565.1_ASM396656v1',
                                          'GCF_004010935.1_ASM401093v1',
                                          'GCF_004214875.1_ASM421487v1',
                                          'GCF_004295685.1_ASM429568v1',
                                          'GCF_004296475.1_ASM429647v1',
                                          'GCF_004296515.1_ASM429651v1',
                                          'GCF_004296535.1_ASM429653v1',
                                          'GCF_004299845.1_ASM429984v1',
                                          'GCF_004323595.1_ASM432359v1',
                                          'GCF_900090215.1_TRABTM',
                                          'GCF_900161835.1_HBAv1',
                                          'GCF_900343015.1_ARADv1'])


# Clean up the indicies
# Some indicies have either suffix delimted by either ' ' or '_'.
idx = data_RM.index
data_RM.index = idx.str.split('(?<!GCF)( |_)').str[0]   # Look for '_' or ' ' which are not proceeded by 'GCF'

idx = data_taxonomy.index
data_taxonomy.index = idx.str.split('(?<!GCF)( |_)').str[0]

# Clean up genus names
data_taxonomy.genus = data_taxonomy.genus.str.replace("'", "")

# Filter by intersection
intersec = data_RM.index.intersection(data_taxonomy.index)
data_taxonomy = data_taxonomy.loc[intersec]
data_RM       = data_RM.loc[intersec]

# Use any connection as indicator
A_ij = data_RM.to_numpy() > 0
A_ij = A_ij.astype(int)

# Create unique ID for each stratification
data_taxonomy['group_id']          = data_taxonomy.groupby('group').grouper.group_info[0]
data_taxonomy['subgroup_id']       = data_taxonomy.groupby('subgroup').grouper.group_info[0]
data_taxonomy['genus_id']          = data_taxonomy.groupby('genus').grouper.group_info[0]
data_taxonomy['species_strain_id'] = data_taxonomy.groupby('species_strain').grouper.group_info[0]

# Save data
np.savetxt(base_path+'A_ij.csv', A_ij, delimiter=',', fmt='%d')
data_taxonomy.to_csv(base_path+'data_taxonomy.csv')