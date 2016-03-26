###All proteomics nano Random Forest custom functions for R
###by Luis F. Montano-Gutierrez. University of Edinburgh 2015.
###All these functions are needed to run the random forest pipelines in the step_by_step scripts in this project
###Copy and paste this document into your R workspace.

### run this script to get all the functions needed for random forest.

###PATH:  /Users/calrician07/Documents/Edinburgh/PhD/1st_Year-masters/second_rotation/data/Maxquant_results/130122/txt_files/


####now we get a color vector out of the training factors

###we use the function true2color next described

true2color<-function(logical, color){
 if(logical==TRUE){
 
 output<- color
 }else{
output<-"black"
}
return(output)
 }
 
 ###################Roc curves
 roc.curve<- function(matrix,step=0.01, true.value='T', false.value='F', color="black", line.width=1, plot.symbol=20){

###first get a matrix with 
#a one column containing observation values ranked from highest to lowest
#bthe second column indicating whether the observation with that value is a T or a F.

###then divide the range of column one into tiny steps.
v.range<- range(as.numeric(matrix[,1]))
steps<- rev(seq(v.range[1], v.range[2], step))
#cat(steps)
###count the number of Ts and Fs
trues<-length(which(matrix[,2] == true.value))
falses<-length(which(matrix[,2]==false.value))
TPs<- rep(NA,length(steps))
FPs<- rep(NA,length(steps))


for (i in 1: length(steps)){
TPs[i]<- length(which(matrix[which(as.numeric(matrix[,1])>= steps[i]),2]==true.value))/trues
FPs[i]<- length(which(matrix[which(as.numeric(matrix[,1])>= steps[i]),2]==false.value))/falses

}
TPs<- c(0,TPs)
FPs<- c(0,FPs)
steps<-c(0,steps)
plot(TPs~FPs, xlim= c(0,1), ylim=c(0,1), type="s", col=color, pch= plot.symbol,lwd=line.width, ylab="Sensitivity", xlab="1-Specificity")
par(new=F)
abline(0,1, col="grey")

return(cbind(steps,TPs,FPs))
}

plot.rocmatrix<-function(matrix, color="black", plot.type="s", plot.ch=NA){
#######function assumes it will get a roc matrix outputted from roc.curve
TPs<-matrix[,2]
FPs<-matrix[,3]
plot(TPs~FPs, xlim= c(0,1), ylim=c(0,1), type=plot.type, col=color, pch=plot.ch)
par(new=F)
abline(0,1, col="grey")


}

##################Convert a class (T or F to color)

class2color.TF<-function(x){
###this function imports a ramp palette
##so it gets a value from 0 to 1 fro x
##and then retrieves that color from the ramp.
output<-NA
#output<- color(1000)[(round(x*100))]
if(x=='?'){
output<- "grey"
}else{
if(x=="F"){
output<-"green"
}
   else{if(x=="T"){
output<-"red"
}

}
}
return(output)
}


###########run random forest

get.trues<- function(table, class.col=1){

return(table[which(table[,class.col]=='T'),])


####selecting experiments at random and look for the ones that create the best classification
}


get.false<- function(table, class.col=1){

return(table[which(table[,class.col]=='F'),])



}


get.unknowns<- function(table, class.col=1){

return(table[which(table[,class.col]=='?'),])



}




get.trainings<- function(table, class.col=1, true.value='T', false.value='F'){

####getting all the lines whose training value is true or false.
 only.trainings<-table[sort(c(grep(true.value, table[,1]), grep(false.value, table[,1]))),]

#as we lost at least one class in the factor, we have to shrink the value''s levels just by redefining the factor.
only.trainings[,class.col]<-as.factor(as.vector(only.trainings[,class.col]))


return(only.trainings)
}
######February 14 2013. Creation of a wrapper function for random forest starting from a data frame.





###a function to get all the training values: that is, all values that are true or false, not '?'.


###the inputs of the function: 
###a data frame containing the training factor variable and all the other variables. The data frame must have column names and rownames
###an indication of which is the training factor

run.random.forest<- function(dataframe, training.factor=1, trees=3000){

###fill all nas with the median values for each column. make sure that the training factor doesn''t declare
###unknowns as NA, otherwise they will be given the name of the most common value.
cat("----filling NAs----\n")

dataframe.no.nas<-na.roughfix(dataframe)

###get all the training cases to create the forest and make a sub data frame with them
cat("----Getting training data frame----\n")
training.frame<-get.trainings(dataframe.no.nas)

###as well as another sub data frame with the unknowns.
cat("----Getting unknowns data frame----\n")
all.unknowns<-get.unknowns(dataframe.no.nas)
 
###create the forest with the training cases.
cat("----Performing forest----\n")

spare.name<-colnames(training.frame)[training.factor]
colnames(training.frame)[training.factor]<- "TF.factor"

the.forest<- randomForest(TF.factor~., data=training.frame, ntrees=trees)

###df.rf <- randomForest(chidito~., data=training.values.T.F.log2, importance=TRUE)
###predict the unknown cases, knowing the percentage of trees used to predict.
cat("----Predicting unknowns----\n")
df.pred<-predict(the.forest, all.unknowns[, (training.factor+1):dim(all.unknowns)[2] ], norm.votes=TRUE, type="prob")
### merge both the unknowns and the training cases.
### add the score of each case to the data frame order it with the trues at the top.
cat("----Returning value----\n")

probs<-rbind(the.forest$votes, df.pred)

probs<- probs[order(as.numeric((rownames(probs)))),]
#print(probs)
#print(dataframe[,1])



#out

out<-list(DATAFRAME=cbind(dataframe, probs), FOREST=the.forest)   
colnames(out$DATAFRAME)[training.factor]<-"training.factor"


return(out)



} 

#####################


####function to calculate the matthews correlation coefficient
###so in order to calculate the matthews correlation coefficient we need to get the confusion matrix:
### that is, the number of TPs, TNs, FPs and FNs.
###we will get the roc matri, that is, predictor values and real classification for each of the cases.
###then, given a threshold above which things will be called positives and below negatives, 
### we will get the confusion matrix measures and therefore we will calculate 
mcc<-function(roc.matrix, threshold, true.values='T', false.values= 'F'){
####we assume the cases are ordered from the highest to the lowest, the highest being the truest.

##true positives will be the cases above the threshold which are called T
tp<- length(which(roc.matrix[which(roc.matrix[,1]>threshold),2]==true.values))

#cat("---",tp, "---\n")
##truenegatives will be the cases below the threshold which are called F
tn<-length(which(roc.matrix[which(roc.matrix[,1]<threshold),2]==false.values))
#cat("---",tn, "---\n")
###false positives will be the cases above the threshold which are called F
fp<- length(which(roc.matrix[which(roc.matrix[,1]>threshold),2]==false.values))
#cat("---",fp, "---\n")
##falsenegatives will be the cases below the threshold which are called T
fn<- length(which(roc.matrix[which(roc.matrix[,1]<threshold),2]==true.values))
#cat("---",fn, "---\n")


numerator<- tn*tp-fp*fn
#pduct<-tp*tn*fp*fn
pduct<-tp
#cat("pduct=" ,pduct, "---\n")
if(pduct==0){
determinator<-1
}else{
determinator<- (log10(tp+fp)+log10(tp+fn)+log10(tn+fp)+log10(tn+fn))/2
determinator=10^determinator
}

return(numerator/determinator)

}


###now we make a wrapper that will include multiple objects as an output, inserted  a list.
###we assume that the first column of the data frame will be the training factor.
### the training factor will contain only values labeled as T, F, or ? depending of its nature. 
### the rest of the columns will be assumed to be predictors (variables)
###the output will be a list containing multiple objects after processing
rf.workflow<- function(dataframe, col="red", max.sensitivity=TRUE){
####run random forest, with output TF predictors properly formated and training factor column called as such.
	rf.result<-run.random.forest(dataframe)
	dataframe.with.predictors<- rf.result$DATAFRAME
 
    the.forest= rf.result$FOREST
	
	importance<-importance(the.forest)
	
	roc.input.matrix<- dataframe.with.predictors[order(dataframe.with.predictors$T, decreasing=T), c("T", "training.factor") ]
	
	
	roc.output.matrix<-roc.curve(roc.input.matrix, color=col)
	par(new=T)
	classif.color<-plot.predictor.quality(dataframe.with.predictors[, "T"],dataframe.with.predictors[,"training.factor"])
	
	mcc.vector<- rep(NA,101)
	names(mcc.vector)<- seq(1,0,-.01)
	
	for( i in 1:101){
		
		mcc.vector[i]<-mcc(roc.input.matrix, threshold=names(mcc.vector)[i])
		
	}
	
	
	max.cutoffs.group<-which(mcc.vector==max(mcc.vector, na.rm=T))
	
	if( max.sensitivity==TRUE ) {
	cutoff<-as.numeric(names(mcc.vector)[max.cutoffs.group[length(max.cutoffs.group)]])
}
	else{
	cutoff<-as.numeric(names(mcc.vector)[max.cutoffs.group[length(max.cutoffs.group)[1]]])
	}
	max.mcc<-max(mcc.vector, na.rm=T)[1]
	
	auc=trapezium.auc(roc.output.matrix, height=.01)
	
	return(list(PREDICTORS.SCORES=dataframe.with.predictors, ROC.INPUT= roc.input.matrix, ROC.OUTPUT=roc.output.matrix, MCC=mcc.vector, MAX.MCC= max.mcc, CUTOFF=cutoff, AUC=auc, IMPORTANCE=importance))
	
}
###Random forest function and predictor quality functions ready.
###Now we want to make an automatic roc matrix creation
###the roc matrix then generates a roc curve, which is plotted.
### the mathews correlation coeffitient is calculated for that roc matrix.

###function to remove a random element from a vector.

remove.element<-function(sample.vector){

index2remove<- sample(1:length(sample.vector),1)
   #cat("index to remove is", index2remove, "\n")   
      if(index2remove==1){
      
      
      return( sample.vector[2:length(sample.vector)])
      }else{
      if(index2remove<length(sample.vector)){
      before<-index2remove-1
      after<- index2remove+1
      
      return(sample.vector[c(1:before, after:length(sample.vector))])
      
      }else{
      if(index2remove==length(sample.vector)){
      end<- length(sample.vector)-1
      return( sample.vector[1:end])
      
      }
      }
      }
      } 




############hill climbing of the mcc value to select the best set of classifiers
#####take a random set of classifiers
#####process the random forest. 

mutate<-function(sample.vector, big.vector){
###the Idea of the mutate function is to take elements of the vector and add, remove, or swap elements of it with elements of the bigger vector.


operation<-sample(c("add","add", "remove", "swap"), 1)



if(operation=="add"){
###this is the number of experiments we cannot exceed
if(length(sample.vector)==26){
return(sample(sample.vector, 12))
}
new.element<-sample(setdiff(sample.vector, big.vector))

return(c(sample.vector, new.element))
}else{
      if(operation=="remove"){
      if(length(sample.vector)<=2){
      
      return(sample(big.vector, 8))
      }else{
      return(remove.element(sample.vector))
         }   
      
      }else{
      
             if(operation=="swap"){
             new.element<-sample(setdiff(sample.vector, big.vector))
            return(c(remove.element(sample.vector), new.element))
             }
           }
     
     
 }    
}

#February 20 2013. To accomplish the goals set in the entry of the same day in the lab journal, I will do the following:


#1. a function that takes a matrix of IDs/names related to factors, and whose output is the IDs/names related to a stated category


get.category.ids<-function(cat.matrix, category){
###[,1] ids. [,2] category factor
return( as.vector(cat.matrix[grep(category, cat.matrix[,2]),1 ]))
}


#2. a function that takes the output IDs/names from above and looks for them in a table. its output is a T F vector that indicates that

###phase 1. just worry about the identifier slot
retrieve.proteins.identified<- function(id.list, prot.group.table){
output<- rep(NA, dim(prot.group.table)[1])
query<-paste(id.list, collapse='|')



output[grep( query, prot.group.table[,1])]<-'T'

return(output)
}

#3. a function to merge two lists of the same size, to call elements coming from one "T" and from another one "F". the NAs will be turned into '?'

merge.vectors<-function(positive, negative){

final.vector<-rep('?', length(positive)) #both are the same length

final.vector[which(negative=='T')]<-'F'
final.vector[which(positive=='T')]<-'T'

return(final.vector)
}

#4. a function to run mcc hill climbing on a protein groups table, which outputs the highest mcc value. and the experiments used to achieve this.

#get.category.ids
#retrieve.proteins.identified
#merge.vectors
#optimize.mcc


optimize.mcc<- function(hill.climbing.training.factor, prot.groups.table, max.increases=15, max.cycles=50, mutate=TRUE){
###the protein groups table must have the experiments that we want to sample.




###this will be our compilation of experiments wich we want to sample
all.experiments<- colnames(prot.groups.table)
####the first experiment sample we will get will be the cell paper experiments (c2-c5) plus the structural ones. just because.
###this table already has



if(length(table(hill.climbing.training.factor))<3){


mcc.hill.climbing<- NA
cat("WARNING:ONLY ONE CATEGORY IN TRAINING SET. RETURN VALUE IS", mcc.hill.climbing, "\n")
return(mcc.hill.climbing)
}else{
mcc.hill.climbing<- list()
sampled.experiments<-all.experiments
forest<-rf.workflow(cbind(hill.climbing.training.factor, log2(prot.groups.table[,sampled.experiments])))
###now we are storing the experiments used for this run


mcc.hill.climbing[[1]]<-list(mcc=max(forest$MCC, na.rm=T),  samples=sampled.experiments, scores=forest$PREDICTORS.SCORES[,dim(forest$PREDICTORS.SCORES)[2]],  training.factor=forest$PREDICTORS.SCORES[,1], roc.output=forest$ROC.OUTPUT)
mcc.hill.climbing[[1]]["cutoff"]<- names(forest$MCC)[which(forest$MCC==mcc.hill.climbing[[1]]$mcc)][1]


cat("Initial mcc value: --", max(forest$MCC, na.rm=T), "--\n")

i<-2
flag<-0
while(i< max.increases&&flag<max.cycles){
cat("Climb # --", i, "--\n")
###we need to refer to the last max mcc

x<- i-1

if(mutate==TRUE){
sampled.experiments<-mutate(mcc.hill.climbing[[x]]$samples, all.experiments)
}
forest<- rf.workflow(cbind(hill.climbing.training.factor, log2(protein.groups[,sampled.experiments])))
##print(forest$PREDICTORS.SCORES[,'T'])
####output from rf.workflow to check what we need.
####return(list(PREDICTORS.SCORES=dataframe.with.predictors, ROC.INPUT= roc.input.matrix, ROC.OUTPTUT=roc.output.matrix, MCC=mcc.vector, CLASSIF.COLORS=classif.color))
mcc.hill.climbing[[i]]<-list(mcc=max(forest$MCC, na.rm=T),  samples=sampled.experiments, scores=forest$PREDICTORS.SCORES[,dim(forest$PREDICTORS.SCORES)[2]],  training.factor=forest$PREDICTORS.SCORES[,1], roc.output=forest$ROC.OUTPUT)
mcc.hill.climbing[[i]]["cutoff"]<- names(forest$MCC)[which(forest$MCC==mcc.hill.climbing[[i]]$mcc)][1]
#cat("cutoff-------\n")
#print( mcc.hill.climbing[[i]]["cutoff"]) 
#cat( "\n")
if(mcc.hill.climbing[[i]]$mcc<=mcc.hill.climbing[[x]]$mcc){
cat("MCC value unsatisfactory. Retrying...\n")

}else{

cat("**Successfully climbed to a higher mcc**: ", mcc.hill.climbing[[i]]$mcc, "\n")
i<-i+1
}

flag<-flag+1
cat("Total number of attempts --", flag, "--\n")
}
#return(mcc.hill.climbing[1:length(mcc.hill.climbing)-1])
return(mcc.hill.climbing)
}
}

########this is a function to "format" the input table, which will have a column of ids related to a column of gene categories.
##what can happen is that a gene is referred to as 2 or more ids, the id can be empty, or the category can be empty.
## therefore what this function wil do is create single entries for each of the duplicate IDs, because: a) we don't know which will appear or will be detected.
#b) in case both point at the same identified protein, there won't be any problem. 

id.category.format.check<-function(matrix){
new.lines<-matrix(ncol=2)

weird.lines<- vector()
c<-1
for(i in 1:dim(matrix)[1]){
#cat("entered for\n")
d<-length(grep("or|;", matrix[i,1]))
if(d>0){
#cat("found an or\n")
new.lines<- rbind(new.lines,cbind(strsplit(matrix[i,1][1], split=" or |;")[[1]], matrix[i,2]))

weird.lines[c]<-i
c<-c+1

}
if(nchar(matrix[i,1])*nchar(matrix[i,2])==0){
#cat("found an empty string\n")
weird.lines[c]<-i
c<-c+1
}
}
#cat("dimensions of new lines:", dim(new.lines),"\n")

#####trim the initial NA NA line from new lines.

trim.length<- dim(new.lines)[1]-1
new.lines<-tail(new.lines, trim.length)
####remove weird lines from the matrix 1: removing any index that might have been annotated twice.
weird.lines<-unique(weird.lines)
#cat("weird lines reduced\n")
matrix<-matrix[setdiff(1:dim(matrix)[1], weird.lines), ]
#cat("set difference performed\n")
###adding the new lines to the formatted matrix
return(rbind(matrix, new.lines))



}




######################id list format check: slice id lines into every independent identifier.

convert.trf.elements<-function(training.factor, elements, convert.to="F"){

training.factor[elements]<-convert.to
return(training.factor)

}

get.trues.vector<- function(vector){

return(which(vector=='T'))


####selecting experiments at random and look for the ones that create the best classification
}

did.it.improve<-function(previous,current){
if( previous==current){

return(FALSE)}else{

if(previous<current){
return(TRUE)

}else{
if(previous>current){
return(FALSE)

}


}
}
}













#########RANDOM FOREST OPTIMIZATION give a training set
oprafor.gts<- function(id.cat.matrix, variable.table, positive, negative, number.of.cycles=100, variable.indices, set.mutate=F){
cat("Formatting table...\n")
formatted.table<- id.category.format.check(id.cat.matrix)
cat("getting IDs from category ", positive, " to make the positive training set...\n")
some.category.ids.positive<-get.category.ids(formatted.table, positive)
cat("Identifying positive proteins in the protein table...")

positive.prots<-retrieve.proteins.identified(id.list=some.category.ids.positive, prot.group.table=variable.table)

######now we have to do it for the negative.

##cell.paper.ids.categories<-cell.paper.protein.table[,c(1,3)]

#formatted.table.negative<- id.category.format.check(cell.paper.ids.categories)
cat("getting IDs from category ", negative, " to make the negative training set...\n")

some.category.ids.negative<-get.category.ids(formatted.table, negative)
cat("Identifying negative training proteins in the protein table...")
negative.prots<- retrieve.proteins.identified(id.list=some.category.ids.negative, prot.group.table=variable.table)
cat("Generating the training factor")





new.training.factor<- merge.vectors(positive.prots,negative.prots)
cat("Proceeding to perform the random forest optimization")


out<-optimize.mcc(new.training.factor, prot.groups.table=variable.table[,variable.indices], max.cycles=number.of.cycles , mutate=set.mutate)

if(is.na(out)==TRUE){
return(out)
}else{

return(tail(out,2)[-2])
}
}





LOO.experiment<-function(rf.object, predictor.table, optimization=NA, cycles=3){

if(is.na(optimization)==FALSE){

optimization<-optimize.mcc(rf.object$training.factor, predictor.table[,rf.object$samples], max.increases=15, max.cycles=15)


}else{
optimization<-rf.object
#cat("optimization is rf.object\n")
}
mcc.to.improve<-optimization[length(rf.object)][[1]]$mcc
#cat("flag1\n")
loo.summary<- matrix(ncol=2, nrow=length(rf.object[length(rf.object)][[1]]$samples)+1)

colnames(loo.summary)<- c("final.mcc","improved")

rownames(loo.summary)<- c(NA,rf.object[length(rf.object)][[1]]$samples)

#cat("flag2\n")
print(dim(loo.summary))

loo.summary[1,]<- c(mcc.to.improve,as.character(did.it.improve(mcc.to.improve,mcc.to.improve)))

output<- list()
#cat("flag22\n")
output$BEST.PREDICTOR.SCORES<- optimization[length(optimization)][[1]]$scores
#cat("flag222\n")
output$BEST.ROC.OUTPUT<- optimization[length(optimization)][[1]]$roc.output
#cat("flag2222\n")
output$BEST.TRAINING.FACTOR<- optimization[length(optimization)][[1]]$training.factor
#cat("flag3\n")



for(i in 1:(dim(loo.summary)[1]-1)){


#cat("flag4\n")
new.experiment.set<-rf.object[length(rf.object)][[1]]$samples[-i]
#cat("flag44\n")
optimization<- optimize.mcc(rf.object[length(rf.object)][[1]]$training.factor, predictor.table[,new.experiment.set], max.cycles=cycles, mutate=F)
#cat("flag444\n")
loo.summary[i+1,c(1,2)]<- c(optimization[length(optimization)][[1]]$mcc, did.it.improve(as.numeric(mcc.to.improve),as.numeric(optimization[length(optimization)][[1]]$mcc)))
#cat("flag5\n")
if(loo.summary[i+1, 2]==TRUE){
output$BEST.PREDICTOR.SCORES<- optimization[length(optimization)][[1]]$scores
output$BEST.ROC.OUTPUT<- optimization[length(optimization)][[1]]$roc.output
output$BEST.TRAINING.FACTOR<- optimization[length(optimization)][[1]]$training.factor
output$BEST.EXPERIMENT.SET<- new.experiment.set
output$LEFT.OUT<- rf.object$samples[-i]
#cat("flag6\n")

}
#cat("flag7\n")
}
#cat("flag8\n")
output$SUMMARY<-loo.summary
return(output)





}


#############
plot.predictor.quality<- function(predictor, response, height=0.05){

#predictor<- rf.test1[,29]
#response<- rf.test1[,1]


response<-as.vector(response)
pred.colors<- sapply(as.vector(response[order(predictor, decreasing=T)]), class2color.TF, USE.NAMES=F)
barplot(rep(height,length(predictor)), border=NA, col=pred.colors, yaxt="n", ylim=c(0,1))
names(pred.colors)<- order(predictor, decreasing=T)
return(pred.colors)
}



###############plot functions


class2color.TF<-function(x, nacolor="grey"){

output<-NA
#output<- color(1000)[(round(x*100))]
if(x=='?'){
output<- nacolor
}else{
if(x=="F"){
output<-"green"
}
   else{if(x=="T"){
output<-"red"
}

}
}
return(output)
}


class2color.TF2<-function(x, colors=c("light blue", "orange"), nacolor="grey"){

output<-NA
#output<- color(1000)[(round(x*100))]
if(x=='?'){
output<- nacolor
}else{
if(x=="F"){
output<-colors[1]
}
   else{if(x=="T"){
output<-colors[2]
}

}
}

class.score2color<-function(tf, score, cutoff,colors=c("light blue", "red", "orange"), nacolor=NA){
output=rep(nacolor, length(score))
for(i in 1:length(tf)){
if(tf[i]=='?'){
output[i]<- nacolor
}else{
if(tf[i]=="F"){
output[i]<-colors[1]
}
   else{if(tf[i]=="T"){
output[i]<-colors[2]
}

}
}

if(is.na(output[i])&& score[i]>cutoff){
output[i]= colors[3]

}

}





return(output)
}


class2color.2TFs<-function(tf1, tf2, refcolor="light blue", truecolors=c("blue", "red"), mixcolor="magenta", nacolor="grey"){

output<-rep(NA, length(tf1))
for( i in 1:length(tf1)){

#output<- color(1000)[(round(x*100))]

if(tf1[i]=='?'||tf2[i]=='?'){
output[i]<- nacolor
}else{
if(tf1[i]=='F'||tf2[i]=='F'){
output[i]<-refcolor
}
}   
if(tf1[i]=="T"&&tf2[i]=="T"){
output[i]<-mixcolor
}else{
if(tf1[i]=="T"){
output[i]<-truecolors[1]
}
if(tf2[i]=="T"){
output[i]<-truecolors[2]
}
}


}



return(output)
}
 class2color.2TFs(cpc.tf, cohesin.tf)




class2pch.2TFs<-function(tf1, tf2, refcolor=21, truecolors=c(22, 24), mixcolor="magenta", nacolor=21){

output<-rep(NA, length(tf1))
for( i in 1:length(tf1)){

#output<- color(1000)[(round(x*100))]

if(tf1[i]=='?'||tf2[i]=='?'){
output[i]<- nacolor
}else{
if(tf1[i]=='F'||tf2[i]=='F'){
output[i]<-refcolor
}
}   
if(tf1[i]=="T"&&tf2[i]=="T"){
output[i]<-mixcolor
}else{
if(tf1[i]=="T"){
output[i]<-truecolors[1]
}
if(tf2[i]=="T"){
output[i]<-truecolors[2]
}
}


}



return(output)
}










#####plot classified proteins: 1 dimension

plot.classified.proteins<-function( rf, name.vector=NA, distance.from.border=.2,plot.character=20, char.size=.5, box="o", offset=0.005, symbol.size=1, plot.trues=T, plot.training=F, na.color=NA){
classifier<-rf$PREDICTORS.SCORES[,"T"]
cutoff<-rf$CUTOFF
training.factor<-rf$PREDICTORS.SCORES[,"training.factor"]

classifier=as.numeric(classifier)
randoms<- sapply(rep(1, length(classifier)), jitter)
#plot(randoms~classifier, col=sapply(classifier, sign2color), pch=plot.character)

par(bty=box)
plot(randoms~classifier,xlim=c(0, 1), xlab="", ylab="", yaxt="n", col= class.score2color(tf=training.factor, score=classifier, cutoff=cutoff), pch=plot.character, cex=symbol.size)

abline(v=cutoff, col="blue", lty="dotted")

for ( i in 1:length(randoms)){


if(plot.training==T){
if(training.factor[i]=="T"){
text(y=randoms[i]-offset, x=classifier[i], name.vector[i], offset=0, cex=char.size)
}}else{
if(plot.trues==T){
if(classifier[i]>0||training.factor[i]=="T"){
text(y=randoms[i]-offset, x=classifier[i], name.vector[i], offset=0, cex=char.size)
}}
}



}
return(cbind(name.vector, classifier,randoms))
}

#### plot classified proteins: 2 dimensions


sign2color<- function(x, colors=c("red", "blue")){
if(x*(-1)<0){
return(colors[1])
}else{
return(colors[2])
}
}


true2color<-function(logical, color){
 if(logical){
 
 output<- color
 }else{
output<-"black"
}
return(output)
 }
 
 
 trainings2colors<- function(tf1,tf2, color1, color2){
 
 tf1<-class2color(tf1, )
  tf1<-class2color(tf1)
 
 
 }

plot.classified.proteins2<-function(rf1, rf2, color.vector=class2color.2TFs(rf1$PREDICTORS.SCORES[,"training.factor"], rf2$PREDICTORS.SCORES[,"training.factor"], nacolor=NA), name.vector=NA, distance.from.border=.2,plot.character=class2pch.2TFs(rf1$PREDICTORS.SCORES[,"training.factor"], rf2$PREDICTORS.SCORES[,"training.factor"], nacolor=NA), char.size=.5, text.size=char.size, box="o", offset=0.005, plot.names=T, plot.trues=T, plot.ids=F, id.vector=NA, border.center=F, significance.gradient=F){

if(border.center==TRUE){
classifier1<-as.numeric(rf1$PREDICTORS.SCORES[,"T"])-as.numeric(rf1$CUTOFF)
classifier2<- as.numeric(rf2$PREDICTORS.SCORES[,"T"])-as.numeric(rf2$CUTOFF)
cutoff1<-0
cutoff2<-0
}else{
classifier1<-as.numeric(rf1$PREDICTORS.SCORES[,"T"])
classifier2<- as.numeric(rf2$PREDICTORS.SCORES[,"T"])
cutoff1<-rf1$CUTOFF
cutoff2<-rf2$CUTOFF

}


par(bty=box)
plot(         classifier2~classifier1, 
               xlim=c(0,1),
               xlab="", ylab="", col=color.vector, pch=plot.character, cex=char.size)

abline(v=cutoff1, h=cutoff2, col="blue", lty="dotted")
if(plot.names==T){
for ( i in 1:length(classifier2)){

if(plot.trues==T){

if(rf1$PREDICTORS.SCORES[i,"training.factor"]=='T'){
##color 
##col= rgb(red=classifier2[i], green=0,blue=classifier1[i])

text(y=classifier2[i]-offset, x=classifier1[i], name.vector[i], offset=0, cex=text.size, col='blue' )
}else{
if(rf2$PREDICTORS.SCORES[i,"training.factor"]=='T'){

text(y=classifier2[i]-offset, x=classifier1[i], name.vector[i], offset=0, cex=text.size,  col= 'red')


}else{

if(significance.gradient==TRUE){

kol= rgb(red=classifier2[i], green=0,blue=classifier1[i])
}else{
kol="black"
}

if(classifier2[i]>cutoff2||classifier1[i]>cutoff1){
text(y=classifier2[i]-offset, x=classifier1[i], name.vector[i], offset=0, cex=text.size, col=kol)


}else{

if(plot.ids==T&&is.na(id.vector)==F){
if(classifier2[i]>cutoff2||classifier1[i]>cutoff1||rf1$PREDICTORS.SCORES[i,"training.factor"]=='T'||rf2$PREDICTORS.SCORES[i,"training.factor"]=='T'){
text(y=classifier2[i]-offset, x=classifier1[i], name.vector[i], offset=0, cex=text.size  )


}

}


}




}

}
}

if((classifier2[i]>cutoff2&&rf2$PREDICTORS.SCORES[i,"training.factor"]=='?')||(classifier1[i]>cutoff1&&rf1$PREDICTORS.SCORES[i,"training.factor"]=='?')){

points(y=classifier2[i]-offset, x=classifier1[i], cex=char.size , col= "grey"  )


}

}




}

if(is.na(id.vector)==F){
return(cbind(id.vector,name.vector,classifier1,classifier2))


}else{

return(cbind(name.vector,classifier1,classifier2))

}
}

plot.classified.proteins2(CCAN.rf, Nup.Ran.rf, name.vector=protein.groups[,8])


####area under the curve calculation from the roc.output

trapezium.area<-function(major.base, minor.base,height){
	
	return(((major.base+minor.base)*height)/2)
	
}

trapezium.auc<- function(roc.output.matrix, height=.001){
	
###we move to the first non zero point of the false positives to avoid problems.
	roc.output.matrix[,3]= round(roc.output.matrix[,3], digits=3)
	roc.output.matrix[,2]= round(roc.output.matrix[,2], digits=3)
    roc.output.matrix=complete.roc.output(roc.output.matrix)
#base1<-0
	
	sum.vector<-vector()
	
	i<-1
	
#for( k in seq(0,1,height)){
		
#base2= min(roc.output.matrix[which(roc.output.matrix[,3]>=k),2])[1]
		
#sum.vector[i]<- ((base1+base2)*height)/2
		
#cat("base 1 .- ", base1,"base 2 .-", base2, "result .- ", sum.vector[i], "\n")
		
#base1<-base2
		
#i<-i+1
		
#}
	
	for(i in 2:nrow(roc.output.matrix))
	
	sum.vector[i-1]<- ((roc.output.matrix[i,2]+roc.output.matrix[i-1,2])*height)/2
	
	
	return(sum(sum.vector))
	
}
#trapezium.auc(cpc.rf$ROC.OUTPUT)

complete.roc.output= function(roc.output.matrix, digits=3){
	
	roc.output.matrix[,3]= round(roc.output.matrix[,3], digits=3)
	roc.output.matrix[,2]= round(roc.output.matrix[,2], digits=3)
	

	step.vector<- seq(0,1, 1/10^digits)
	tps<-rep(NA, length(step.vector))
	first.nonzero.fp=min(which(roc.output.matrix[,3]!=0))
	k<-1

	i<-first.nonzero.fp
	
		last.tp<- roc.output.matrix[first.nonzero.fp,2]
	while(i<=nrow(roc.output.matrix)){
#cat('entered')
	
	if(roc.output.matrix[i,3]<step.vector[k]){
	i<-i+1
		next
	}
	
	while(roc.output.matrix[i,3]>step.vector[k]){
		tps[k]<-last.tp

	k<-k+1	
		
		
	}
	
	if(roc.output.matrix[i,3]==step.vector[k]){
#cat('entered')
		tps[k]<-roc.output.matrix[i,2]
		last.tp<-roc.output.matrix[i,2]
		
		k<-k+1
	}
	
	i<-i+1

	}
	
return(rbind(c(0,0),cbind(step.vector,tps)))
	
}
#x=complete.roc.output(cpc.rf$ROC.OUTPUT)
#plot(x[,2]~x[,1], type='s')






generateRandomClassifiers<- function(numObservations, numClassifiers){
out=runif(numObservations)
for( i in 1:numClassifiers){

cbind(out, runif(numObservations))
}
return(out)
}



% go protein by protein in target table. look for the column on source table. add it on final vetor

reindex.experiment<-function(source.table, target.table,target.table.name.index, source.table.name.index, source.table.value.index){

out<- rep(NA, nrow(target.table))
for(i in 1: nrow(target.table)){
 result<-grep (target.table[i,target.table.name.index], source.table[,source.table.name.index])
cat(result,"\n\n")
if(length(result)==0){
out[i]<-NA
}else{
out[i]<-result

}


}
return(out)

}

reindex.experiment(rotation.protein.groups, protein.groups, 1,1,212)





random.ISCoF.iterations<- function( dataset=1, rows=1, cols=1, querySize=10, refSize=425, numEvents=10000){
if (length(dataset==1)){
toyDataset<-matrix(nrow= rows, ncol=cols, data=runif(rows*cols)  )

rows<-nrow(toyDataset)
cols=ncol(toyDataset)

}else{

toyDataset<-dataset
rows<-nrow(dataset)
cols=ncol(dataset)
}



aucvec<-rep(NA, numEvents)
mccvec<-rep(NA, numEvents)
for(i in 1:numEvents){
centroids<- sample(1:rows, querySize+refSize)
query<- centroids[1:querySize]
reference<-centroids[querySize+1:refSize]
####out of these centroids the last 240 will be the reference. the first 60 we divide them into 5 groups of 12.
###we make the values of each sequence of 12 similar to the first one of the sequence, but slightly jittered.


queryIndices<- indicestoTF(query, toyDataset)
referenceIndices<- indicestoTF(reference, toyDataset)
fake1.tf<-merge.vectors(positive=queryIndices, negative=referenceIndices)
fake1.rf<-rf.workflow(cbind(fake1.tf, as.data.frame(toyDataset)))

mccvec<-c(mccvec, fake1.rf$MAX.MCC)
aucvec<-c(aucvec, fake1.rf$AUC)
}

return(list(MAX.MCC=mccvec, AUC=aucvec))
}


random.ISCoF<- function( dataset=1, rows=1, cols=1, querySize=10, refSize=425, numEvents=10000){
if (length(dataset)==1){
toyDataset<-matrix(nrow= rows, ncol=cols, data=runif(rows*cols)  )


}else{

toyDataset<-dataset

}
rowz<-nrow(toyDataset)
colz<-ncol(toyDataset)

#print(rows)
#print(cols)
aucvec<-rep(NA, numEvents)
mccvec<-rep(NA, numEvents)

centroids<- sample(1:rowz, querySize+refSize) ##These centroids will be both positive and negative training sets.
query<- centroids[1:querySize]
reference<-centroids[querySize+1:refSize]
####out of these centroids the last chunk will be the reference. the first 60 we divide them into 5 groups of 12.
###we make the values of each sequence of 12 similar to the first one of the sequence, but slightly jittered.


queryIndices<- indicestoTF(query, toyDataset)
referenceIndices<- indicestoTF(reference, toyDataset)
fake1.tf<-merge.vectors(positive=queryIndices, negative=referenceIndices)
fake1.rf<-rf.workflow(cbind(fake1.tf, as.data.frame(toyDataset)))




return(fake1.rf)
}



#random.ISCoF(dataset=protein.groups[,experiment.ratio.column.indices])

randomTrainingSets<- random.ISCoF(dataset=protein.groups[,experiment.ratio.column.indices], numEvents=5000)



ISCoF.iterations<-function(tf, dataset, numEvents=500){
out<-rep(NA, numEvents)
for(i in 1:numEvents){
rf<-rf.workflow(cbind(tf, as.data.frame(dataset)))
out[i]<- rf$MAX.MCC

}
return(out)

}


indicestoTF<- function(indices, proteingroups){

out<- rep(NA, nrow(proteingroups))
out[indices]<-'T'
return(out)
}


plotMCCdensities<-function(queryVector, randomVector, bandw=.001){
plot(density(queryVector, bw= bandw), xlim=c(0,1), col='red', main='')
par(new=T)
plot(density(randomVector, bw= bandw), xlim=c(0,1), main='',col='blue', xlab='', ylab='', xaxp=NULL, yaxp=NULL)

}




plotErrorBars<-function(matr){

errbar(x=1:ncol(matr), y=colMeans(matr), yplus=colMeans(matr)+apply(matr,1, sd), yminus=colMeans(matr)-apply(matr,1, sd))

}



