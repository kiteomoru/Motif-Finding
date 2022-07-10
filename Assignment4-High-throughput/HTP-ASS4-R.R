###HTP-Homework4
#Relative Path
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


#Required libraries
library(dplyr)
library(tidyverse)
library(reshape)
library(stringr)
require(data.table)
library(bioseq)


##Q1
#load data 
pwmfile<- read.table('argR-counts-matrix.txt', header = F, sep = "\t",  row.names =1)
pwmfile<- pwmfile[,2:19]
pwmfile<-as.matrix(pwmfile)
 view(pwmfile)
#Create PWN function  for working out the position weight matrix value
pwm<- function(matrix){
  
  #Frequency with pseudocounts
p<-matrix + 1

  freq<-list()
  for (i in 1:nrow(p)){
    for (j in 1:ncol(p)){
      freq<- c(freq, p[i,j]/ sum(p[,j]))
    }
  }
  #Weight matrix
  weight<-list()
  for (i in freq){
    weight<-c(weight,log(i/.25))
  }
  
  weightmat<-matrix(unlist(weight), ncol=18, nrow =4,byrow=TRUE) 
  rownames(weightmat)<-c('a','c','g','t')
  View(weightmat)
  return(weightmat)
}

scoringmatrix<-pwm(pwmfile)
scoringmatrix


#write table
write.table(scoringmatrix, file = 'scoringmatrix.txt', quote = FALSE,sep='\t', col.names = NA)



##Q2
#Load the data
seq<- read.table('E_coli_K12_MG1655.400_50 3', as.is = T)
seq <- subset(seq, select = -c(V2,V4))
View(seq)

#####Generate sequences of length 18 (This code will take a while. Please be patient)
kmer18=dna(seq$V3)
kmer18=seq_split_kmer(kmer18, k = 18)
kmer18=as.data.frame(kmer18)

head(kmer18)
seqkmers<- kmer18 %>%
  mutate_all(as.character)
# Change sequences from upper case to lower case 
seqkmers<-data.frame(lapply(seqkmers, function(v) {
  if (is.character(v)) return(tolower(v))
  else return(v)
}))

colnames(seqkmers)= seq$V1 #add the gene ids as column names
head(seqkmers)

# Create a Function that scans the weights in the scoring matrix(Logo-odds Matrix) from q1 to identify binding sites in the data-frame
weight = function(kmerstr){
  x <- strsplit(x=kmerstr,split='')
  seq_score <- vector()
  for (i in 1:18){
    seq_score[i] <- scoringmatrix[x[[1]][i],i]
    
  }
  return(print(seq_score))
}

# apply function

#Reshape data
newcol= c(1:nrow(seqkmers))
seqkmers$newcol=newcol

reshaped_seqkmers<- melt(seqkmers, 
                         id.vars = c('newcol'),
                         variable_name = c('Seq_ID'),
                         value.name= c('Sequences'))

View(reshaped_seqkmers)
# apply function to add scores to the data frame
reshaped_seqkmers$scores<- lapply(reshaped_seqkmers$Sequences, weight) #Takes a while
View(reshaped_seqkmers)
####

#Create a new column(totalscore) to the dataframe which is the sum of all the scores
reshaped_seqkmers$totalscores<- lapply(reshaped_seqkmers$scores, function(a) sum(as.numeric(a)))
View(reshaped_seqkmers)
###
#remove the newcol colume from the dataframe
reshaped_seqkmers <- subset(reshaped_seqkmers, select = -c(newcol))
View(reshaped_seqkmers)

#Obtain the highest score for each gene id
topscoresforallgeneids<-setDT(reshaped_seqkmers)[, .SD[which.max(totalscores)], by= variable]
topscoresforallgeneids$totalscores<-as.numeric(topscoresforallgeneids$totalscores)

View(topscoresforallgeneids)
#####
# Another way to view top scores per gene-ID
reshaped_seqkmers$totalscores<-as.numeric(reshaped_seqkmers$totalscores)
Z<- reshaped_seqkmers %>%
  group_by(variable) %>%
  summarize(max.totalscores = max(totalscores))
View(Z)

#Top30
topscoresforallgeneids_arranged <- topscoresforallgeneids %>% arrange(desc(totalscores))
view(topscoresforallgeneids_arranged)

##top30
top30motifs<- head(topscoresforallgeneids_arranged, n=30)
view(top30motifs)

##Write table
top30motifs2<- subset(top30motifs, select = -c(scores))
write.table(top30motifs2, file = 'top30motifs.txt', quote = FALSE,sep='\t', col.names = NA)



####Reference:
##The results in this assignment was made with guidance from codes found in the website: 
##https://davetang.org/muse/2013/10/01/position-weight-matrix/ 
## and various references to stackoverflow




