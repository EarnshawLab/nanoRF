###Script to run analysis step by step
### 0. include the relevant libraries. if not installed, run:
install.packages(RandomForest)
###otherwise proceed with
library(randomForest)
###1.loading protein groups table. substitute the path to import the right table.
protein.groups.path<- "/Users/calrician07/Dropbox/PhD/Edinburgh/PhD/1st_Year-masters/second_rotation/data/shinya/20150320_shinyadata/20140702_All SMC KO proteinGroups2.txt"
protein.groups<-read.delim(protein.groups.path) 
###the previous line assumes that the delimiter is tab and that the first row is the column names (as in maxquant's protein groups result. in case the field delimiter is a comma, you must add sep=',' to the command.
###protein groups is assumed to be a data frame because it has different kinds of info.
#
###2. selecting the columns that contain the relevant ratios for each experiment.
###to get them, it is suggested to do colnames(protein.groups). then one can add them into the c() below.
###for this particular experiment, we are interested in using the raw ratios first. these are (chosen from colnames):
###[69] "Ratio.H.L.CAPD3_4"      
###[73] "Ratio.H.L.CAPD3_5" 
###[77] "Ratio.H.L.CAPH_1"
###[81] "Ratio.H.L.CAPH_2"  
###[85] "Ratio.H.L.SCC1_2"   
###[89] "Ratio.H.L.SCC1_3"
###[93] "Ratio.H.L.SCC1_4"   
### [97] "Ratio.H.L.SMC2_1" 
###[101] "Ratio.H.L.SMC2_4"   
###[105] "Ratio.H.L.SMC5_1"   
###[109] "Ratio.H.L.SMC5_3"  
 
experiment.ratio.column.indices<-c(69,73,77,81,85,89,93,97,101,105,109)
#
###3. importing protein cohort training sets from a file and creating training factors.
####### in the protein groups table, we want to identify the proteins we will work with:
####### that is, our initial guesses of protein cohorts and proteins that aren't so.
#######However, doing it by hand would be a pain in the ass every time, specially if we are unsure.
###to automate this, we will:
###a) create string vectors that contain all the possible IDs of every kind of our proteins of interest
###b) store the instructions to make those vectors in a script so we can 'renovate' again and again later.
#######this has been done in a file like "load_training_set_ids.R") where the IDs of condensin, cohesin, smc5/6, cpc and cytosolic proteins can be found as vectors.
###c) actually run the script to get the vectors in our workspace
#######in our particular case, the vector script is loaded by running these files (modify the path accordingly):
source("/Users/calrician07/Dropbox/PhD/Edinburgh/PhD/1st_Year-masters/second_rotation/data/shinya/20150320_shinyadata/load_training_damy_set_ids20140717forLuis.R")
source("/Users/calrician07/Dropbox/PhD/Edinburgh/PhD/1st_Year-masters/second_rotation/data/shinya/20150320_shinyadata/load_training_set_ids20141007.R")
####### Running this script brought the following objects into the workspace:
# "cytosol.raw.ids"
# "condensin.trimmed.ids"       
#  "condensinI.trimmed.ids"            
#  "condensinII.trimmed.ids" 
#  "cpc.ids" 
#  "smc56.trimmed.ids" 
#  "scaffold.trimmed.ids"  
#kinetochore.trimmed.ids
#telomere.trimmed.ids
#histone.trimmed.ids
#CCAN.trimmed.ids
#Nup.Ran.trimmed.ids
# ribosomes.ids
#cytosol.raw.ids
#histoneH1.ids
   
######(having them like this will allow us to freely and painlessly modify our ID list to include the proteins we want, and then we just rerun the import line again.)

###d) search the protein groups table for the rows that contain the proteins of each ID list.
######to this end we will use the retrieve.proteins.identified function, which just searches the protein.groups table and retrieves the rows
######that contain the indices of the matching rows. we will work with these index lists later on.
######get our index lists through the following command. doesn't matter if the lists are redundant (i.e. contain repeated row numbers) for now.
cytosol.protIndices<- retrieve.proteins.identified(cytosol.raw.ids, protein.groups)
condensin.protIndices<- retrieve.proteins.identified(condensin.trimmed.ids, protein.groups)
condensinI.protIndices<- retrieve.proteins.identified(condensinI.trimmed.ids, protein.groups)
condensinII.protIndices<- retrieve.proteins.identified(condensinII.trimmed.ids, protein.groups)
cpc.protIndices<- retrieve.proteins.identified(cpc.ids, protein.groups)
smc56.protIndices<- retrieve.proteins.identified(smc56.trimmed.ids, protein.groups)
cohesin.protIndices<- retrieve.proteins.identified(cohesin.trimmed.ids, protein.groups)
kinetochore.protIndices<-retrieve.proteins.identified(kinetochore.trimmed.ids, protein.groups)
telomere.protIndices<-retrieve.proteins.identified(telomere.trimmed.ids, protein.groups)
histone.protIndices<-retrieve.proteins.identified(histone.trimmed.ids, protein.groups)
CCAN.protIndices<-retrieve.proteins.identified(CCAN.trimmed.ids, protein.groups)
Nup.Ran.protIndices<-retrieve.proteins.identified(Nup.Ran.trimmed.ids, protein.groups)
ribosomes.protIndices<-retrieve.proteins.identified(ribosomes.ids, protein.groups)
histoneH1.protIndices<- retrieve.proteins.identified(histoneH1.ids, protein.groups)

###
###e) create training factors indicating TRUE proteins of a class, FALSE proteins not belonging to a class and UNKNOWN proteins
######for this, we need to merge the indices of 2 protein groups from above. thus, 
###### one group will be the positive and the other will be the negative training set.
###### the rest of the proteins must be declared as unknown. Importantly, the function must assign
###### the negative training set first so the indices of our group of interest (positive training) are not lost/overwritten.
##### for example, if an important protein was misassigned as negative, because we have it in the positive set
##### it will remain in our training. The function returns a vector with the same number of elements as the rows in protein groups, 
###where each element is T (positive), F(negative) or ?(unassigned) according to our wishes. for each hypothesis, we need a training set.
###we thus get these training factors for each of the protein complexes, with cytosol proteins as a negative training set in every case.
cpc.tf<-merge.vectors(positive=cpc.protIndices, negative=cytosol.protIndices)
condensin.tf<-merge.vectors(positive=condensin.protIndices, negative=cytosol.protIndices)
condensinI.tf<-merge.vectors(positive=condensinI.protIndices, negative=cytosol.protIndices)
condensinII.tf<-merge.vectors(positive=condensinII.protIndices, negative=cytosol.protIndices)
cohesin.tf<- merge.vectors(positive=cohesin.protIndices, negative=cytosol.protIndices)
scaffold.tf<-merge.vectors(positive=scaffold.protIndices, negative=cytosol.protIndices)
smc56.tf<-merge.vectors(positive=smc56.protIndices, negative=cytosol.protIndices)
kinetochore.tf<-merge.vectors(positive=kinetochore.protIndices, negative=cytosol.protIndices)
telomere.tf<-merge.vectors(positive=telomere.protIndices, negative=cytosol.protIndices)
histone.tf<-merge.vectors(positive=histone.protIndices, negative=cytosol.protIndices)
CCAN.tf<-merge.vectors(positive=CCAN.protIndices, negative=cytosol.protIndices)
Nup.Ran.tf<-merge.vectors(positive=Nup.Ran.protIndices, negative=cytosol.protIndices)
ribosomes.tf<-merge.vectors(positive=ribosomes.protIndices, negative=cytosol.protIndices)
histoneH1.tf<-merge.vectors(positive=histoneH1.protIndices, negative=cytosol.protIndices)

### 4. running the random forest workflow for each training/hypothesis
###Now that we have each of the training factors we just need to run the actual analysis. In short, we need to:
#a) run the actual random forest algorithm which will give probabilities for each protein (row) to belong to the positive training set selected.
#in the best case, all the members of the complex in question will have the highest values, and all the negatives will be  below the lowest of the positives.
#this means there is a perfect split between the two classes. in the worst case, both classes will be copletely intertwined. 
#b) Calculate the quality of the category split for each possible "border" , to choose the best border/threshold.
#Because the RF will retrieve probability of belonging to T for each protein, we must decide on a border.
#How to choose that border? the best border for a RF will be the one that best divides our Ts from our Fs.
#For example, maybe for a RF all Ts are at .7 or above and all Fs are below . then .7 is the best border.
#To calculate this automatically, we need to evaluate all possible borders from 0 to 1 in an efficient way. 
#That is the purpose of the MCC. at some point within the probability range we will maximize the MCC. that probability (the max MCC) will automatically be chosen as our border.
#c)store the max MCC value to compare separation quality.
#if we were to order all possible cases  from the worst to the best case, their max MCCs would order from 0 (worst case) to 1 (best case).
#Thus, the max MCC value is informative of our success. 
#d)return all the information in an ordered way.
#The function rf.workflow does all this and its output contains all the important info. Keep in mind
#that its input is a simplified table that only contains the training factor and the columns of the ratios.

cpc.rf<-rf.workflow(cbind(cpc.tf, protein.groups[,experiment.ratio.column.indices]))
condensin.rf<-rf.workflow(cbind(condensin.tf, protein.groups[,experiment.ratio.column.indices]))
condensinI.rf<-rf.workflow(cbind(condensinI.tf, protein.groups[,experiment.ratio.column.indices]))
condensinII.rf<-rf.workflow(cbind(condensinII.tf, protein.groups[,experiment.ratio.column.indices]))
scaffold.rf<-rf.workflow(cbind(scaffold.tf, protein.groups[,experiment.ratio.column.indices]))
smc56.rf<-rf.workflow(cbind(smc56.tf, protein.groups[,experiment.ratio.column.indices]))
cohesin.rf<-rf.workflow(cbind(cohesin.tf, protein.groups[,experiment.ratio.column.indices]))
kinetochore.rf<-rf.workflow(cbind(kinetochore.tf, protein.groups[,experiment.ratio.column.indices]))
telomere.rf<-rf.workflow(cbind(telomere.tf, protein.groups[,experiment.ratio.column.indices]))
histone.rf<-rf.workflow(cbind(histone.tf, protein.groups[,experiment.ratio.column.indices]))
CCAN.rf<-rf.workflow(cbind(CCAN.tf, protein.groups[,experiment.ratio.column.indices]))
Nup.Ran.rf<-rf.workflow(cbind(Nup.Ran.tf, protein.groups[,experiment.ratio.column.indices]))
ribosomes.rf<-rf.workflow(cbind(ribosomes.tf, protein.groups[,experiment.ratio.column.indices]))
histoneH1.rf<-rf.workflow(cbind(histoneH1.tf, protein.groups[,experiment.ratio.column.indices]))


###5. export results into a table
#We want to produce a table that has all the relevant information. This includes, for each experiment, which proteins were used as training factors and the RF (RF) scores.



resultsTable<-cbind(
protein.groups[,c(8,9)],
cpc.rf$PREDICTORS.SCORES[,c("training.factor","T")],
condensin.rf$PREDICTORS.SCORES[,c("training.factor","T")],
condensinI.rf$PREDICTORS.SCORES[,c("training.factor","T")],
condensinII.rf$PREDICTORS.SCORES[,c("training.factor","T")],
smc56.rf$PREDICTORS.SCORES[,c("training.factor","T")],
cohesin.rf$PREDICTORS.SCORES[,c("training.factor","T")],
kinetochore.rf$PREDICTORS.SCORES[,c("training.factor","T")],
telomere.rf$PREDICTORS.SCORES[,c("training.factor","T")],
histone.rf$PREDICTORS.SCORES[,c("training.factor","T")],
CCAN.rf$PREDICTORS.SCORES[,c("training.factor","T")],
Nup.Ran.rf$PREDICTORS.SCORES[,c("training.factor","T")],
ribosomes.rf$PREDICTORS.SCORES[,c("training.factor","T")],
histoneH1.rf$PREDICTORS.SCORES[,c("training.factor","T")]


)

###Because this table has non-specific names, we create a vector with names we want it to have.
colnamevec<- c(
"protein.name",
"fasta.header",
"cpc.training.factor",
"cpc.RF.score",
"condensin.training.factor",
"condensin.RF.score", 
"condensinI.training.factor",
"condensinI.RF.score",
 "condensinII.training.factor",
 "condensinII.RF.score",
  "smc56.training.factor",
  "smc56.RF.score", 
  "cohesin.training.factor",
  "cohesin.RF.score",
    "kinetochore.training.factor",
  "kinetochore.RF.score",
   "telomere.training.factor", 
    "telomere.RF.score",  
    "histone.training.factor", 
    "histone.RF.score",  
    "CCAN.training.factor", 
    "CCAN.RF.score", 
    "NupRan.training.factor", 
    "NupRan.RF.score",   
    "ribosomes.training.factor", 
    "ribosomes.RF.score", 
    "histoneH1.training.factor", 
    "histoneH1.RF.score" 
  
)

###we change the column names of the table to match the custom names that we created above.
colnames(resultsTable)<- colnamevec

write.csv(resultsTable, file="20150321_RF_results_table.csv")
