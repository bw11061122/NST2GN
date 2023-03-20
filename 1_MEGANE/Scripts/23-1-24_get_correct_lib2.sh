#!/usr/bin/env bash
# I noticed that the file from Pio had either 8 spaces or a tab
# In this script, I change every instance 8 spaces occur to a tab
## Note: this worked and fixed the script I was running!


# I am using this tutorial to swap 8 spaces for a tab
# https://www.geeksforgeeks.org/sed-command-in-linux-unix-with-examples/

# first I am making a copy to work on
cp Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph.rep Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_copy.rep 

# This replaces 8 spaces with a tab
sed -e 's/        /\t/g' Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_copy.rep > Lib3/lib2_renamed_nodups_fixedunkowns_plusgraph_copy_tab.rep

# In the fofn, change Pio's crsID to mine
sed -e 's/pas211/bw450/g' data/M.zeb_on_AstCal_67.fofn > data/M.zeb_on_AstCal_67_bw450.fofn

## Create a test file with 5 sequences to run on slurm
head -10 data/M.zeb_on_AstCal_67_bw450.fofn > data/test_M.zeb_on_AstCal_67_bw450.fofn