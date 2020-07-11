#----------------------Get extension of the input files-----------------
getExtension <- function(file){ 
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
} 
#---------------------Converting user specified input files to bedGraph----------------
user_specified_converting_input_formats <- function(args)
{
  
  path <- paste(args[1],"/Converted_inputs",sep="")
  dir.create(path)
  #------------- Replicate 1 for training the model----------
  train_rep1_format <- getExtension(args[3])
  if(train_rep1_format=="bedGraph"){
    system(paste("cp" ,args[3] ,"rep1_train_model.bedGraph"))
  }else if(train_rep1_format!="bedGraph"){
    if(train_rep1_format=="bigWig"){
      system(paste("./bigWigToBedGraph " ,args[3] ," rep1_train_model.bedGraph -chrom=",as.name(args[5]),sep=""))
    }else if(train_rep1_format=="bam"){
      system(paste("bedtools genomecov -ibam" ,args[3], "-bga > rep1_train_model.bedGraph"))
    }
  }
  #------------- Replicate 2 for training the model----------
  train_rep2_format <- getExtension(args[4])
  if(train_rep2_format=="bedGraph"){
    system(paste("cp" ,args[4] ,"rep2_train_model.bedGraph"))
  }else if(train_rep2_format!="bedGraph"){
    if(train_rep2_format=="bigWig"){
      system(paste("./bigWigToBedGraph " ,args[4] ," rep2_train_model.bedGraph -chrom=",as.name(args[5]),sep=""))
    }else if(train_rep2_format=="bam"){
      system(paste("bedtools genomecov -ibam" ,args[4], "-bga > rep2_train_model.bedGraph"))
    }
  }
  #------------- Signals to be variance stabilized-----------
  test_rep1_format <- getExtension(args[6])
  if(test_rep1_format=="bedGraph"){
    system(paste("cp" ,args[6] ,"instabilized_signals.bedGraph"))
  }else if(test_rep1_format!="bedGraph"){
    if(test_rep1_format=="bigWig"){
      system(paste("./bigWigToBedGraph " ,args[6] ," instabilized_signals.bedGraph -chrom=",as.name(args[7]),sep=""))
    }else if(test_rep1_format=="bam"){
      system(paste("bedtools genomecov -ibam" ,args[6], "-bga > instabilized_signals.bedGraph"))
    }
  }
  system(paste("mv rep1_train_model.bedGraph ","'",path,"'/",sep=""))
  system(paste("mv rep2_train_model.bedGraph ","'",path,"'/",sep=""))
  system(paste("mv instabilized_signals.bedGraph ","'",path,"'/",sep=""))
}

#-----------------------Converting instabilized user provided signals to bedGraph----------
default_converting_input_formats <- function(args)
{
  
  path <- paste(args[1],"/Converted_inputs",sep="")
  dir.create(path)
  #------------- Signals to be variance stabilized-----------
  test_rep1_format <- getExtension(args[6])
  if(test_rep1_format=="bedGraph"){
    system(paste("cp" ,args[6] ,"instabilized_signals.bedGraph"))
  }else if(test_rep1_format!="bedGraph"){
    if(test_rep1_format=="bigWig"){
      system(paste("./bigWigToBedGraph " ,args[6] ," instabilized_signals.bedGraph -chrom=",as.name(args[7]),sep=""))
    }else if(test_rep1_format=="bam"){
      system(paste("bedtools genomecov -ibam" ,args[6], "-bga > instabilized_signals.bedGraph"))
    }
  }
  
  system(paste("mv instabilized_signals.bedGraph ","'",path,"'/",sep=""))
}
#-------------------------------

Weighted_Mean_1_over_sigma<-function(rep1,rep2,distance,bin,alpha,ordered_scores)
{
  
  L=ceiling(length(ordered_scores[,1])/bin)
  l=bin
  Y=matrix(nrow=L,ncol = 1)
  a=matrix(nrow=L,ncol = 1)
  z=matrix(nrow=((2*distance)+1),ncol = 1)
  weighted_z=matrix(nrow=((2*distance)+1),ncol=1)
  n=(2*distance)+1
  mean_1_over_sigma=matrix(nrow=L,ncol=2)
  #########################
  for(i in 1:(L-1))
  {
    Y[i,1]=sum(ordered_scores[(((i-1)*l)+1):(i*l),2])
    a[i,1]=sum((ordered_scores[(((i-1)*l)+1):(i*l),2])^2)
  }
  Y[L,1]=sum(ordered_scores[(((L-1)*l)+1):length(ordered_scores[,1]),2])
  a[L,1]=sum((ordered_scores[(((L-1)*l)+1):length(ordered_scores[,1]),2])^2)
  ##########################
  for (j in 1:n)
  {
    z[j,1] <- (1/(alpha^(abs(distance+1-j)*l)))
    weighted_z[j,1] <- l*(1/(alpha^(abs(distance+1-j)*l)))
  }
  ##########################
  for(i in 1:L)
  {
    output=Preprocess_Weighted_Mean_1_over_sigma(Y,z,weighted_z,a,i,distance)
    mean_1_over_sigma[i,1]=output[1]
    mean_1_over_sigma[i,2]=output[2]
  }
  mean_1_over_sigma
}
###############################
Preprocess_Weighted_Mean_1_over_sigma<- function(weighted_sum,weights,sum_of_weights,sum_square,center_bin,width)
{
  
  epsilon=0.0001
  y=weighted_sum
  W=weighted_sum
  L=length(y)
  z=weights
  weighted_z=sum_of_weights
  k=center_bin
  w=width
  a=sum_square
  A=sum_square
  
  lower=max(1,(k-w))
  upper=min(L,(k+w))
  xx=upper-lower+1
  center_z=ceiling(length(z)/2)
  lower_z=(center_z-(k-lower))
  counter=lower_z
  upper_z=(center_z+(upper-k))
  for(i in lower:upper)
  {
    y[i,1]=y[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_mean=sum(y[lower:upper])/sum(weighted_z[lower_z:upper_z])
  counter=lower_z
  for(i in lower:upper)
  {
    a[i,1]=a[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_var=sum(a[lower:upper])/sum(weighted_z[lower_z:upper_z])
  weighted_var=weighted_var -(weighted_mean^2)#best
  weighted_var=weighted_var+epsilon
  mean_1_over_sigma_values <- c(weighted_mean,1/sqrt(weighted_var))
  mean_1_over_sigma_values
}

######################
VSS_train_user_specified_replicates <- function(args)
{
  
  Input_path <- args[1]
  main_path=paste(Input_path,"/Converted_inputs",sep="")
  setwd(main_path)
  
  whole_rep1=read.table("rep1_train_model.bedGraph")
  whole_rep2=read.table("rep2_train_model.bedGraph")
  
  #train_chr <- 22
  #chromosomes <- c(21)
  #chromosomes <- args[4:length(args)]
  #chromosomes <- c()
  # counter <- 1
  # for(chr in args[4]:args[length(args)])
  #   {
  #   print(chr)
  #     x <- strsplit(chr,"chr")
  #     chromosomes[counter] <- as.integer(x[[1]][2])
  #     counter <- counter+1
  #   }
  x <- strsplit(args[5],"chr")
  chromosomes <- as.integer(x[[1]][2])
  
  path=paste(Input_path,"/trained_models",sep="")
  dir.create(path)  
  ############ Get the training data for calculating mean-variance curve ############
  for(chr in chromosomes)
  {
    print(chr)
    train_chr <- chr
    rep1<- subset(whole_rep1, whole_rep1[,1]==(paste("chr",train_chr,sep="")))
    rep2<- subset(whole_rep2, whole_rep2[,1]==(paste("chr",train_chr,sep="")))
    n_row=min(rep1[length(rep1[,1]),3],rep2[length(rep2[,1]),3])
    rep1_score=matrix(nrow=n_row,ncol=1) #Scores for all genomic positions
    rep2_score=matrix(nrow=n_row,ncol=1)
    rep11=subset(rep1, rep1[,3]<=n_row)
    rep22=subset(rep2, rep2[,3]<=n_row)
    l_rep11=length(rep11[,3])
    l_rep22=length(rep22[,3])
    
    ####################
    if(rep11[l_rep11,3]!=n_row)
      if(rep1[l_rep11+1,2]<n_row)
        
      {
        rep11[l_rep11+1,]=rep1[l_rep11+1,]
        rep11[l_rep11+1,3]=n_row
      }
    rep1=rep11
    ###################
    if(rep22[l_rep22,3]!=n_row)
      if(rep2[l_rep22+1,2]<n_row)
        
      {
        rep22[l_rep22+1,]=rep2[l_rep22+1,]
        rep22[l_rep22+1,3]=n_row
      }
    rep2=rep22
    ##################
    
    
    for(i in 1:length(rep1[,4]))
    {
      rep1_score[rep1[i,2]:rep1[i,3],1]<-rep1[i,4]
      
    }
    for(i in 1:length(rep2[,4]))
    {
      rep2_score[rep2[i,2]:rep2[i,3],1]<-rep2[i,4]
      
    }
    a=c(rep1_score,rep2_score)
    b=c(rep2_score,rep1_score)
    scores=matrix(nrow=length(a),ncol=2)
    scores[,1]  <- a
    scores[,2]  <- b
    keep_rows <- morder(scores, 1, decreasing = FALSE)
    ordered_scores <- scores[keep_rows,]
    ordered_scores <- data.frame(ordered_scores)
    L=length(ordered_scores[,1])
    #----------------
    zero_signals=subset(ordered_scores,ordered_scores[,1]==0)
    non_zero_signals=subset(ordered_scores,ordered_scores[,1]!=0)
    threshold <- 0.25*L
    
    if(length(zero_signals[,1])>=threshold){
      
      #---------mean-var curve smoothing parameters ---------
      beta_values=1000
      bin_sizes=100000
      alpha_values=2^(1/as.integer(beta_values))
      width_values=ceiling(log(1/0.01)/log(alpha_values))
      #------------------------------------------------------
      
      mean_ordered <- mean(ordered_scores[,1])
      less_than_mean <- subset(ordered_scores,ordered_scores[,1]<=mean_ordered)
      higher_than_mean <- subset(ordered_scores,ordered_scores[,1]>mean_ordered)
      Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(rep1,rep2,width_values,as.integer(bin_sizes),alpha_values,higher_than_mean)
      Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- mean(less_than_mean[,2])
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd(less_than_mean[,2])
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)+1),1] <- 0
      Mean_1_over_sigma[(nrow(Mean_1_over_sigma)),2] <- 1/sd(less_than_mean[,2])
      Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
      
    }else if(length(zero_signals[,1])<threshold){
      #---------mean-var curve smoothing parameters ---------
      beta_values=1e+07
      bin_sizes=1000
      alpha_values=2^(1/as.integer(beta_values))
      width_values=ceiling(log(1/0.01)/log(alpha_values))
      Mean_1_over_sigma <- Weighted_Mean_1_over_sigma(rep1,rep2,width_values,as.integer(bin_sizes),alpha_values,ordered_scores)
      Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    }
    
    write.table(Mean_1_over_sigma, file=paste(path,"/chr",chr,"_trained_model.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
    
  }
}
#############################
VSS_transform_user_specified <- function(args)
{
  
  Input_path <- args[1]
  x <- strsplit(args[7],"chr")
  chromosomes <- as.integer(x[[1]][2])
  
  main_path=paste(Input_path,"/Converted_inputs",sep="")
  setwd(main_path)
  instabilized_replicate_signals=read.table("instabilized_signals.bedGraph")
  
  path=paste(Input_path,"/Output",sep="")
  dir.create(path)  
  Replicate1_scores=c()
  for(chr in chromosomes)
  {
    Mean_1_over_sigma <- read.table(paste(Input_path,"/trained_models","/chr",chr,"_trained_model",'.bedGraph',sep=""))
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    test_rep1<- subset(instabilized_replicate_signals, instabilized_replicate_signals[,1]==(paste("chr",chr,sep="")))
    replicate1_scores <- calculating_vss_signals(Mean_1_over_sigma,test_rep1)
    Replicate1_scores <- rbind(Replicate1_scores,replicate1_scores)
  }
  
  write.table(Replicate1_scores, file=paste(path,"/Variance_stabilized_signals.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  
}

#######################
VSS_transform_default <- function(args)
{
  
  Input_path <- args[1]
  Source_path <- args[8]
  x <- strsplit(args[7],"chr")
  chromosomes <- as.integer(x[[1]][2])
  
  main_path=paste(Input_path,"/Converted_inputs",sep="")
  setwd(main_path)
  instabilized_replicate_signals=read.table("instabilized_signals.bedGraph")
  
  path=paste(Input_path,"/Output",sep="")
  dir.create(path)  
  Replicate1_scores=c()
  for(chr in chromosomes)
  {
    #Mean_1_over_sigma <- read.table(paste(Input_path,"/default_trained_models","/chr",chr,"_trained_model",'.bedGraph',sep=""))
    Mean_1_over_sigma <- read.table(paste(Source_path,"/bin/default_trained_models","/chr",chr,"_trained_model",'.bedGraph',sep=""))
    Mean_1_over_sigma <- data.frame(Mean_1_over_sigma)
    test_rep1<- subset(instabilized_replicate_signals, instabilized_replicate_signals[,1]==(paste("chr",chr,sep="")))
    replicate1_scores <- calculating_vss_signals(Mean_1_over_sigma,test_rep1)
    Replicate1_scores <- rbind(Replicate1_scores,replicate1_scores)
  }
  write.table(Replicate1_scores, file=paste(path,"/Variance_stabilized_signals.bedGraph",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
}

##########################
calculating_vss_signals <- function(Mean_1_over_sigma,rep1)
{
  
  replicate1_scores <- rep1
  replicate1_scores <- replicate1_scores[order(replicate1_scores[,4]),]
  unique_scores <- unique(replicate1_scores[,4])
  
  Mean_vs_sigma <- Mean_1_over_sigma
  Mean_vs_sigma <- data.frame(Mean_vs_sigma)
  Mean_vs_sigma[,1] <- Mean_1_over_sigma[,1]
  Mean_vs_sigma[,2] <- 1/Mean_1_over_sigma[,2]
  
  #-----------spline -----------
  fit_1 <- smooth.spline(Mean_vs_sigma[,1], (Mean_vs_sigma[,2]), cv=TRUE)
  predicted_sigma <- stats:::predict.smooth.spline(fit_1,replicate1_scores[,4])$y
  
  #-----------------------
  predicted_score=matrix(nrow=length(predicted_sigma),ncol=2)
  predicted_score[,1]=replicate1_scores[,4]
  predicted_score[,2]=predicted_sigma
  
  if(min(predicted_score[,2])<0){
    less <- subset(predicted_score,predicted_score[,1]<min(Mean_vs_sigma[,1]))
    if(min(less[,2])<0){
      less[,2] <- min(Mean_vs_sigma[,2])
    }
    
    middle <- subset(predicted_score,predicted_score[,1]>=min(Mean_vs_sigma[,1]) & predicted_score[,1]<=max(Mean_vs_sigma[,1]))
    high <- subset(predicted_score,predicted_score[,1]>max(Mean_vs_sigma[,1]))
    if(min(high[,2])<0){
      high[,2] <- max(Mean_vs_sigma[,2])
    }
    
    predicted_score <- rbind(less,middle,high)
  }
  predicted_score[predicted_score[,2]<0,2] <- mean(Mean_vs_sigma[,2])
  predicted_score[,2]=1/predicted_score[,2]
  
  x=predicted_score[,1]
  y=predicted_score[,2]
  ordered_data <- predicted_score[order(x,y),]
  ordered_data=ordered_data[!duplicated(ordered_data),]
  ordered_data=data.frame(ordered_data)
  
  for(i in 1:length(unique_scores))
  {
    
    scores_list <- which(replicate1_scores[,4]==unique_scores[i])
    start <- scores_list[1]
    end <- scores_list[length(scores_list)]
    
    sub <- subset(ordered_data, ordered_data[,1]<=unique_scores[i]) #recently
    score <-trapz(sub[,1],sub[,2])
    replicate1_scores[start:end,5] <- score
    
  }
  replicate1_scores<- replicate1_scores[order(replicate1_scores[,2]),]
  replicate1_scores[,4] <- NULL
  replicate1_scores
}



#-----------Preprocessing----------------
if(args[2]=="user_specified"){
  setwd(args[1])
  library(bigmemory)
  library(pracma)
  user_specified_converting_input_formats(args)
  VSS_train_user_specified_replicates(args)
  VSS_transform_user_specified(args)
}else if(args[2]=="default"){
  setwd(args[1])
  library(bigmemory)
  library(pracma)
  default_converting_input_formats(args)
  VSS_transform_default(args)
}
#-------------------------------------

