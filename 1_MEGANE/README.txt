##### README
##### 2023-03-16 
##### Barbara Walkowiak 
##### NST2GN Part II Project 

##### PACKAGE VERSIONS / SOFTWARE
megane_v1.0.1.beta.sif
https://github.com/shohei-kojima/MEGAnE 

In this folder, I included scripts which were used to run the MEGANE pipeline
In order to generate a vcd file of polymorphic TE insertions for 2764 Malawi cichlid samples 

#### Inputs 
The lists to the input files (cichlid genome samples) are located in Fofn_2764 folder
The TE annotation library used in the pipeline is located in NST2GN_FINAL/Data/TE annotation library folder and was provided by Pio Sierra
To run the pipeline, I used AstCal1.2 reference genome (genome folder provided by Pio Sierra)

#### Outputs
The output of the pipeline is located in data/vcf folder and includes two vcf files (insertions shared with the reference and insertions present in the reference are included in separate vcd folders)

#### Scripts
MEGANE step 0: megane_step0_a_calliptera.sh 
Modify library prep: 23-1-24_get_correct_lib2.sh
MEGANE step 1 prep: 23-1-23_megane_launcher_step1_prep_a_calliptera.sh
MEGANE step 1: 23-1-24_test_megane_launch.sh
MEGANE step 2: 23-1-26_megane_step2.sh
MEGANE step 3: 23-1-27_megane_step3.sh

#### Fofn_2674
These are the list of file names that were provided from Bettina Fischer 
Which I used to generate the vcf file
(On the cluster, they are present in the "data" folder)

#### Lists of directories
These files contain list of directories which I generated as part of the MEGANE pipeline 