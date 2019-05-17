# step1: SD
require(stringr)

## Author:  Zhezhen Wang
## Email: zhezhen@uchicago.edu
## last update:  4/1/2019
## Acknowledgement:  National Institutes of Health  R21LM012619  (Xinan H Yang)

#' The length of a string (in characters).
#'
#' @param df a count matrix with unique loci row names:ID X column names:samples 
#' @param samplesL a list of characters with stages as names
#' @param cutoff numeric value, if < 1 automaticlly goes to select top x% transcripts 
#' if > 1 goes to select by x-fold more than the selected method (which is either the reference, other stages or pervious stage)
#' default 0.01
#' @param method select from 'reference','other', or 'previous'
#' for 'reference', the reference has to be the first 
#' for 'previous', make sure sampleL is in the right order from benign to malign
#' default uses 'other' 
#' @return numeric vector giving number of characters in each element of the
#'   character vector.  Missing strings have missing length.
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

sd_selection = function(df,samplesL,cutoff = 0.01,method = 'other'){
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  tmp = names(samplesL)
  samplesL = lapply(samplesL,as.character)
  test2 = sapply(tmp, function(x) apply(df[,samplesL[[x]]],1,sd,na.rm = T))
  
  if(method == 'reference'){
    ref = as.character(samplesL[[1]])
    sdref = apply(df[,ref],1,sd,na.rm = T)
    sds = lapply(tmp,function(x) test2[,x]/sdref)
    names(sds) = tmp
    
  }else if(method == 'other'){
    othersample = lapply(1: length(samplesL), function(x) do.call(c,samplesL[-x]))
    names(othersample) = tmp
    sdother = sapply(tmp, function(x) apply(df[,as.character(othersample[[x]])],1,sd,na.rm = T))
    
    sds = lapply(tmp,function(x) test2[,x]/sdother[,x])
    names(sds) = tmp
    
  }else if(method == 'previous'){
    warning('Using method "previous", make sure sampleL is in the right order')
    sds = lapply(2:ncol(test2),function(x) test2[,x]/test2[,x-1])
    tmp = tmp[-1]
    names(sds) = tmp
    
  }else{
    stop("method need to be selected from 'reference','other','previous'")
  }
  
  if(cutoff<1){
    topdf = nrow(df)*cutoff
    sdtop = lapply(tmp,function(x) names(sds[[x]][order(sds[[x]],decreasing = T)[1:topdf]]))
  }else{
    sdtop = lapply(tmp,function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }
  
  if(method == 'reference') tmp = tmp[-1]
  names(sdtop) = tmp
  subdf = lapply(tmp,function(x) df[,as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf), function(x) subset(subdf[[x]],row.names(subdf[[x]]) %in% sdtop[[x]]))
  names(subm)  = tmp
  return(subm)
}


#' simulation of sd selection
#' 
#' @description simulation of sd selection
#'
#' @param sampleL a list of characters with stages as names
#' @param lociL a named list of matrix for each states can be obtained from sd_selection
#' @param method select from 'reference','other', or 'previous'
#' for 'reference', the reference has to be the first 
#' for 'previous', make sure sampleL is in the right order from benign to malign
#' default uses 'other'  
#' @param B a integer, number of times to run this simulation
#' @param percent 
#' @param cutoff numeric value, if < 1 automaticlly goes to select top x% transcripts 
#' if > 1 goes to select by x-fold more than the selected method (which is either the reference, other stages or pervious stage)
#' default 0.01
#' @return a list of igraph object
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

sd_selection.simulation = function(samplesL,lociL,locibk,method = 'other',B=100,percent=0.8,cutoff=0.01){
   N.random = lapply(1:length(lociL), function(x) matrix(0, nrow = length(locibk),ncol=B))
   for(i in 1:length(N.random)){
     row.names(N.random[[i]]) = locibk
   }
   
   n = lengths(sampleL)
   k = n*percent
#   X <- nrow(counts)*top
#  Y <- sapply(lociL,nrow)
  
  for(i in c(1:B)) {
    random_sample = sapply(1:length(k),function(x) sample(1:n[[x]],k[[x]]))  # replace=F by default
    selected_counts = sapply(1:length(samplesL), function(x) lociL[[x]][,random_sample[[x]]])
    test2 = sapply(selected_counts, function(x) apply(x,1,sd,na.rm = T))
    tmp = names(lociL)
    colnames(test2) = tmp
    
  if(method == 'reference'){
    ref = selected_counts[[1]]
    sdref = apply(ref,1,sd,na.rm = T)
    sds = lapply(tmp,function(x) test2[,x]/sdref[,x])
    names(sds) = tmp
    
  }else if(method == 'other'){
    othersample = lapply(1: length(samplesL), function(x) do.call(c,samplesL[-x]))
    names(othersample) = tmp
    sdother = sapply(tmp, function(x) apply(df[,as.character(othersample[[x]])],1,sd,na.rm = T))
    
    sds = lapply(tmp,function(x) test2[,x]/sdother[,x])
    names(sds) = tmp
    
  }else if(method == 'previous'){
    warning('Using method "previous", make sure sampleL is in the right order')
    sds = lapply(2:ncol(test2),function(x) test2[,x]/test2[,x-1])
    names(sds) = tmp
    
  }else{
    stop("method need to be selected from 'reference','other','previous'")
  }
  
  if(cutoff<1){
    topdf = nrow(df)*cutoff
    sdtop = lapply(tmp,function(x) names(sds[[x]][order(sds[[x]],decreasing = T)[1:topdf]]))
  }else{
    sdtop = lapply(tmp,function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }
  names(sdtop) = tmp
    names(N.random) = tmp
    for(j in tmp){
      N.random[[j]][sdtop[[j]],i] = 1
    }
  }
  return(N.random)
} 


# sd_selection.simulation = function(counts,sdother,stage,cli,method = c('reference','other','previous'),B=100,percent=0.8,top=0.01){ 
#   #counts: a matrix of numvers
#   #B: an interger of simulation times
#   #percent: a numeric of sub-group to run simulation
#   #stage: a character of one of the colnames in the 'cli' data.frame
#   #cli: a data.frame with a column named 'group_xy' and the first column contains the colnames of counts
#   #top: a numeric, when usd, to be the top percentage
#   
#   if(!all(cli$Row.names %in% colnames(counts))) stop('check cli: must has a column named "Row.names" which contain the colnames of counts')
#   
#   #N.random = matrix(nrow=nrow(counts),ncol=B)
#   N.random = matrix(0, nrow=nrow(counts),ncol=B)
#   row.names(N.random) = row.names(counts)
#   
#   stage_counts = counts[,as.character(subset(cli,group == stage)[,'Row.names'])] 
#   
#   n = ncol(stage_counts)
#   k = n*percent
#   X <- nrow(counts)*top
#   Y <- ncol(stage_counts)
#   
#   for(i in c(1:B)) {
#     random_sample = sample(1:Y,k)  # replace=F by default
#     selected_counts = stage_counts[,random_sample]
#     test2 = apply(selected_counts,1,sd,na.rm = T) #
#     tmp3 = test2/sdother
#     tmp4 = tmp3[order(tmp3,decreasing = T)[1:X]]
#     N.random[names(tmp4),i] = 1
#   }
#   return(N.random)
# }


# step2: get network 
require(stringr)
require(psych)
require(igraph)

#' get a igraph object based on Pearson Correlation Coefficiency(PCC)
#' 
#' @description get a igraph object based on Pearson Correlation Coefficiency(PCC)
#'
#' @param optimal a list of count matrix
#' @param p a numeric cutoff
#' @return a list of igraph object
#' @export
#' @examples
#' 
#' 
#' @importFrom stringr psych igraph
#' @author Zhezhen Wang

# first version in F:\ZheZhen\GSE49711\code\DNB_PCCnetwork_top1.R
# add in p=0.05 param 4/4/2019
getNetwork = function(optimal,p = 0.05){
  rL = lapply(optimal,function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
  names(rL) = names(optimal)
  pL = lapply(optimal,function(x) corr.test(t(x),adjust = 'fdr',ci=F)$p)
  if(is.null(names(rL))) stop('give names to the input list')
  
  igraphL = list()
  for(i in names(rL)){
    test = rL[[i]]
    test.p = pL[[i]]
    test[lower.tri(test,diag = T)] = NA
    #test.p[lower.tri(test,diag = T)] = 1
    tmp = lapply(1:nrow(test),function(x) test[x,test.p[x,]<p])
    tmp_name = lapply(1:nrow(test),function(x) which(test.p[x,]<p))
    idx = which(lengths(tmp_name)==1)
    for(j in idx){
      names(tmp[[j]]) = names(tmp_name[[j]])
    }
    names(tmp) = row.names(test)
    edges = stack(do.call(c,tmp))
    edges = subset(edges, !is.na(values))
    tmp2 = subset(edges,grepl('\\.[1-9,A-z]\\.',ind))
    if(nrow(tmp2)!=0){
      tmp2$node1 = paste0(str_split_fixed(tmp2$ind,'\\.',3)[,1],'.',str_split_fixed(tmp2$ind,'\\.',3)[,2])
      tmp2$node2 = str_split_fixed(tmp2$ind,'\\.',3)[,3]
    }
    edges = subset(edges,!grepl('\\.[1-9,A-z]\\.',ind))
    edges$node1 = str_split_fixed(edges$ind,'\\.',2)[,1]
    edges$node2 = str_split_fixed(edges$ind,'\\.',2)[,2]
    edges = rbind(edges,tmp2)
    dim(edges) #[1] 1270    4
    dim(edges) #[1] 583   4
    edges = edges[,c('node1','node2','values')]
    edges$weight = abs(edges$values) # added in 1/8/2019
    #colnames(edges) = c('node1','node2','weight') # added in 12/18/2018
    
    nodes = data.frame(unique(c(edges$node1,edges$node2)))
    print(paste0(i,':',nrow(nodes),' nodes')) #[1] 48    1
    routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
    igraphL[[i]] = routes_igraph
  }
  return(igraphL)
}

## probably do not need anymore 4/4/ 2019
getNetwork_nocutoff = function(optimal){
  rL = lapply(optimal,function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
  
  igraphL = list()
  for(i in names(rL)){
    test = rL[[i]]
    test[lower.tri(test,diag = T)] = NA
    tmp = lapply(1:nrow(test),function(x) test[x,!is.na(test[x,])])
    names(tmp[[length(tmp)-1]]) = colnames(test)[ncol(test)]
    names(tmp) = row.names(test)
    edges = stack(do.call(c,tmp))
    edges = subset(edges, !is.na(values))
    tmp2 = subset(edges,grepl('\\.[1-9,A-z]\\.',ind))
    if(nrow(tmp2)!=0){
      tmp2$node1 = paste0(str_split_fixed(tmp2$ind,'\\.',3)[,1],'.',str_split_fixed(tmp2$ind,'\\.',3)[,2])
      tmp2$node2 = str_split_fixed(tmp2$ind,'\\.',3)[,3]
    }
    edges = subset(edges,!grepl('\\.[1-9,A-z]\\.',ind))
    edges$node1 = str_split_fixed(edges$ind,'\\.',2)[,1]
    edges$node2 = str_split_fixed(edges$ind,'\\.',2)[,2]
    edges = rbind(edges,tmp2)
    #dim(edges) #[1] 1270    4
    #dim(edges) #[1] 583   4
    edges = edges[,c('node1','node2','values')]
    edges$weight = abs(edges$values) # added in 1/8/2019
    #colnames(edges) = c('node1','node2','weight') # added in 12/18/2018
    
    nodes = data.frame(unique(c(edges$node1,edges$node2)))
    print(paste0(i,':',nrow(nodes))) #[1] 48    1
    routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
    igraphL[[i]] = routes_igraph
  }
  return(igraphL)
}

# first version in F:\ZheZhen\GSE49711\code\DNB_PCCnetwork_top1.R
getCluster = function(igraphL,steps = 4){
  #   if(is.na(step)){
  #     step = sapply(igraphL,function(x) nrow(as_data_frame(x, what="vertices"))*2)
  #     print(step)
  #   }else 
  if(length(steps) == 1 & steps %% 1 == 0){               ##grepl("^[1-9]{1,}$", step) only works for 1 digit
       steps = rep(steps,length(igraphL))
     }else if(length(steps) != 1 | length(steps) != length(igraphL)){
       stop('check step: must be postive integer(s) of length 1 or length of igraphL')
     }
  groups = list()
  for(i in 1:length(igraphL)){
    if(nrow(as_data_frame(igraphL[[i]])) != 0){
      groups[[i]] = cluster_walktrap(igraphL[[i]],weight = abs(E(igraphL[[i]])$weight),steps = steps[i])
    }else{
      groups[[i]] = NA
    }
  }
  #groups = lapply(1:length(igraphL), function(x) cluster_walktrap(igraphL[[x]],weight = abs(E(igraphL[[x]])$weight),steps = steps[x])) # changed weight to abs(PCC) 12/18/2018
  names(groups) = names(igraphL)
  return(groups)
}

#' get clusters of nodes by clustering methods
#' 
#' get clusters of nodes by clustering methods
#'
#' @param igraphL a list of count matrix or a list of igraph object
#' @param method a character need to be select from 'rw', 'hcm', 'km', 'pam', or 'natrual'. Uses 'rw' by default
#' for 'hcm', using complete. for 'km' using 'euclidean'.
#' @param cutoff a numeric value, default NULL
#' @param countsL a list of count matrix
#' @return A list of vertex id vectors that can be used in plot function in igraph package as 'mark.groups' parameter
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

getCluster_methods = function(igraphL, method = 'rw',cutoff = NULL){
  if(method == 'rw'){
    if(all(sapply(igraphL,class) != 'igraph')) 
      stop('random walk clustering needs a list of igraph object which can be obtained using getNetwork')
    if(!is.null(cutoff)) if(cutoff%%1 !=0) warning('Please provide a integer as "cutoff" for the cluster method random walk')
    if(is.null(cutoff)) cutoff = 4
    groups = getCluster(igraphL,cutoff)
  }else if(method == 'hcm'){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame'))) 
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r) 
    groups = lapply(1:length(testL), function(x) hclust(dist(testL[[x]]), method = "complete"))
  }else if(method %in% c('km','pam')){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame'))) 
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
    groups = lapply(1:length(testL), function(x) KMedoids(testL[[x]],3,distance = 'euclidean'))
  }else if(method == 'natrual'){
    if(all(sapply(igraphL,class) != 'igraph')) 
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    groups = lapply(1:length(igraphL), function(x) components(igraphL[[x]])$membership)
    
  }else(stop('please select from "rw", "hcm","km", "pam", "natrual" as method'))
  return(groups)
}


# step3: calulate CI
require(psych)
require(igraph)

#' get CI score
#' 
#' @description get CI score
#'
#' @param groups A named list of vertex id vectors
#' @param countsL a named list of numeric count matrix which can be obtained through sd_selection
#' @param plot a boolean to decide if plot a bar plot of CI scores or not. Default TRUE
#' @param adjust.size a boolean to decide if CI score should be adjust by size or not. Default FALSE
#' @param ylim a vector needed if the output barplots need to be on the same y scale 
#' @param nr the number of rows to plot
#' @param nc the number of column to plot, default length(groups)
#' @param order the order of the barplot. default NULL which is using the input list order
#' @return a list of CI score and their components
#' @export
#' @examples
#' 
#' 
#' @importFrom psych igraph
#' @author Zhezhen Wang

# first version in F:\ZheZhen\GSE49711\code\DNB_PCCnetwork_top1.R
getCI = function(groups,countsL,plot = TRUE,adjust.size = FALSE,ylim = NULL,nr=1,nc = NULL,order = NULL){
  if(all(is.na(groups))){
    warning('no loci in any of the state in the list given, please rerun getCluster_methods with a larger cutoff or provide a list of loci')
  }else{
    if(all(sapply(groups,class) =="communities")){
      membersL = lapply(groups,membership)
    }else if(any(is.na(groups)) & any(sapply(groups,class) =="communities")){
      removed = groups[is.na(groups)]
      groups = groups[!is.na(groups)]
      membersL = lapply(groups,membership)
    }else{
      membersL = groups
    }
    CIl = PCCol = PCCl = sdl =list()
    names(membersL) = names(groups) # probably need to be changed to names(groups) instead of names(countsL)
    if(is.null(names(groups))) warning('No names provided for "groups"')
    if(is.null(names(countsL))) warning('No names provided for "countsL"')
    if(is.null(nc)) nc = length(groups)
    if(plot) par(mfrow =c(nr,nc))
    
    if(is.null(order)){
      loop = names(membersL)
    }else{
      if(any(!order %in% names(membersL))) stop('make sure all names provided in "order" is included in names of "countsL"')
      loop = order
    }
    
    for(i in loop){
      test = membersL[[i]]
      if(all(is.na(test))){
        CI = sdL = PCC = PCCo = NA
      }else{
        test.counts = countsL[[i]]
        m = lapply(1:max(test),function(x) subset(test.counts, row.names(test.counts) %in% names(test[test==x])))
        comple = lapply(1:max(test),function(x) subset(test.counts, !row.names(test.counts) %in% names(test[test==x])))
        names(m) = names(comple) = letters[1:max(test)]
        if(length(m)>26) names(m) = names(comple) = paste0(letters[1:max(test)],1:max(test))
        
        #PCCo = mapply(function(x,y) abs(cor(t(x),t(y))),comple,m)
        PCCo = lapply(names(m), function(x) abs(cor(t(comple[[x]]),t(m[[x]]))))
        PCCo_avg = sapply(PCCo,mean) 
        
        PCC = lapply(m,function(x) abs(cor(t(x))))
        #for(j in 1:length(PCC)){
        #  PCC[[j]][lower.tri(PCC[[j]],diag = T)] = NA
        #}
        PCC_avg = sapply(PCC,function(x) (sum(x,na.rm = T)-nrow(x))/(nrow(x)^2-nrow(x)))
        #PCC_avg = sapply(PCC,function(x) mean(x,na.rm = T))
        sdL = lapply(m, function(x) apply(x,1,sd))
        if(adjust.size){
          CI = mapply(function(x,y,z,w) mean(x)*(y/z)*sqrt(nrow(w)), sdL,PCC_avg,PCCo_avg,m)
        }else{
          CI = mapply(function(x,y,z) mean(x)*(y/z), sdL,PCC_avg,PCCo_avg)
        } 
        
        if(plot) {
          tn = ifelse(is.null(order),names(membersL)[i],i)
          cex = ifelse(length(CI)>20,0.7,1) ## added 1/9/2019
          ## 3/1/2019 changed legend 'n=' to #letters=
          barplot(CI,legend = paste0('#',names(m),'=',sapply(m,nrow)),col = rainbow(length(m), alpha = 0.3),
                  main = tn,ylab = 'CI',xlab='modules',args.legend = list(cex = cex),ylim =ylim) 
        }
      }
      CIl[[i]] = CI
      sdl[[i]] = sdL
      PCCl[[i]] = PCC
      PCCol[[i]] = PCCo
    }
    if(is.null(order)){
      names(CIl) = names(PCCol) = names(PCCl) = names(sdl) = names(groups)
    }else{
      membersL = membersL[order]
    }
    return(list(members = membersL,CI = CIl,sd = sdl,PCC=PCCl,PCCo=PCCol))
  }
}

#' get the index and members with the maximum CI score
#' 
#' @description get the index and members with the maximum CI score
#'
#' @param membersL a list of characters, the first element of output from function getCI
#' @param CIl alist of numerics, the secned element of output from function getCI
#' @param minsize a numeric value of minimum size allowed for a cluster
#' @return a list of index and members
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

getMaxCImember = function(membersL,CIl,minsize = 1){
  listn = names(membersL)
  if(!minsize <1){
    minsize = minsize-1
    CIl = lapply(1:length(membersL),function(x) ifelse(table(membersL[[x]])>minsize,CIl[[x]],NA))
    module_keep = lapply(1:length(membersL), function(x) names(table(membersL[[x]])[table(membersL[[x]])>(minsize-1)]))
    membersL = lapply(1:length(membersL),function(x) membersL[[x]][membersL[[x]] %in% module_keep[[x]]])
  }else(stop('please provide a minimum size for the cluster, which should be integer that is larger than 0'))
  
  idx = sapply(CIl,which.max)
  maxCI = sapply(1:length(idx),function(x) names(membersL[[x]][membersL[[x]] == idx[x]]))
  names(maxCI) = listn
  names(idx) = listn
  return(list(idx = idx,members = maxCI))
}

#' get the index and members with the maximum CI score
#' 
#' @description get the index and members with the maximum CI score
#' 
#' @param membersL a list of characters or numbers, any element of the result from function getCI.
#' @param idx of clusters 
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

getMaxStats = function(membersL,idx){
  member_max = sapply(names(idx),function(x) membersL[[x]][idx[[x]]])
  names(member_max) = names(idx)
  member_max = sapply(member_max,mean)
  return(member_max)
}

#' plot the maximized CI score for each stage
#' 
#' @description plot a line plot with maximum CI score for each stage
#' 
#' @param maxCIms a list of 2 elements, 1st one is a stage-named idex vector,2nd one is a list of members in each stages; can be obtained by run function getMaxCImember
#' @param membersL a list of charactersand numbers, can be obtained from function getCI
#' @param las the direction of x-labels. default 0. labels are parallel (=0) or perpendicular(=2) to axis
#' @param order the ordered stages names
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

plotMaxCI = function(maxCIms, membersL, las = 0,order = NULL){
  if(is.null(maxCIms[[1]]) | is.null(membersL[[1]])) stop('Please give names for the second element of the "maxCIms" and "membersL"')
  if(is.null(order)){
    CI = sapply(names(maxCIms[[1]]), function(x) membersL$CI[[x]][maxCIms[[1]][[x]]])
  }else{
    if(any(!order %in% names(maxCIms[[2]]))) stop('make sure all names in "order" are in names of the 2nd element of maxCIms')
    if(any(!names(maxCIms[[2]]) %in% order)) warning('not every state in "simulation" is plotted, make sure "order" is complete')
    CI = sapply(order, function(x) membersL$CI[[x]][maxCIms[[1]][[x]]])
  }
  if(class(CI) == 'list'){
    warning('changing NA CI score(s) to 0')
    idx = sapply(CI,function(x) length(x)==0)
    CI[idx] = 0
    CI = do.call(c,CI)
  }
  matplot(CI,type = 'l',ylab = 'CI(m|r)',axes=F)
  if(is.null(order)){
    len = lengths(maxCIms[[2]])
  }else{
    len = sapply(order,function(x) length(maxCIms[[2]][[x]]))
  }
  text(1:length(maxCIms[[1]]),CI+0.01,paste0('(',len,')'))
  
  axis(2)
  # customize x-axis
  if(is.null(order)){
    stages = names(maxCIms[[1]])
  }else{
    stages = order
  }
  axis(side=1,at=1:length(maxCIms[[1]]),labels=stages,las = las)
}

getCI_inner = function(members,countsL,adjust.size){
  #len = ifelse(class(members) == "numeric",members,length(members))
  random_id = sample(nrow(countsL[[1]]),members)
  randomL = lapply(names(countsL), function(x) countsL[[x]][random_id,])
  comple = lapply(names(countsL), function(x) subset(countsL[[x]],!row.names(countsL[[x]]) %in% row.names(randomL[[x]])))
  names(randomL) = names(comple) = names(countsL)
  PCCo = lapply(names(countsL), function(x) abs(cor(t(comple[[x]]),t(randomL[[x]]))))
  PCCo_avg = sapply(PCCo,function(x) mean(x,na.rm = T)) 
  PCC = lapply(randomL,function(x) abs(cor(t(x))))
  PCC_avg = sapply(PCC,function(x) (sum(x,na.rm = T)-nrow(x))/(nrow(x)^2-nrow(x)))
  sdL = lapply(randomL, function(x) apply(x,1,sd))
  
  if(adjust.size){
    CI = mapply(function(x,y,z,w) mean(x)*(y/z)*sqrt(nrow(w)), sdL,PCC_avg,PCCo_avg,m)
  }else{
    CI = mapply(function(x,y,z) mean(x)*(y/z), sdL,PCC_avg,PCCo_avg)
  } 
  return(CI)
}

#' get Ic score
#' 
#' @description get Ic score see (ref)
#' 
#' @param len a number of length of DNS 
#' @param countsL a named list of numeric count matrix for each states 
#' @param adjust.size a boolean to decide if CI score should be adjust by size or not. Default FALSE
#' @param B number of times to run this simulation, default 1000
#' @param order the ordered stages names
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

simulationCI = function(len,samplesL,df,adjust.size = FALSE,B=1000){
  if(is.null(names(samplesL))) stop('please provide names for list countsL')
  countsL = lapply(samplesL, function(x) df[,x])
  if(is.null(names(countsL))) names(countsL) = names(samplesL)
  m = sapply(1:B, function(x) getCI_inner(len,countsL,adjust.size))
  row.names(m) = names(countsL)
  return(m)
}

### can be removed ##########
simulationCI_1st = function(len,countsL,adjust.size = FALSE,B=1000){
  #if(is.null(names(samplesL))) stop('please provide names for list countsL')
  #countsL = lapply(samplesL, function(x) df[,x])
  if(is.null(names(countsL))) names(countsL) = names(samplesL)
  m = sapply(1:B, function(x) getCI_inner(len,countsL,adjust.size))
  row.names(m) = names(countsL)
  return(m)
}

#' simulation of loci to calculate the CI score
#' 
#' @description simulation of loci to calculate the Ic score
#' 
#' @param CI 
#' @param simulation a matrix state * number of simulated times, can be obtained from simulation_loci
#' @param order the ordered stages names
#' @param B run time of simulation, default 1000 
#' @return a line plot of simulated 
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

plot_CI_Simulation = function(CI,simulation,las = 0,order = NULL,ylim = NULL,B=100){
  #toplot = simulation
  if(!is.null(order)){
    if(any(!order %in% row.names(simulation))) stop('make sure "simulation" has row.names which are in "order"')
    if(any(!row.names(simulation) %in% order)) warning('not every state in "simulation" is plotted, make sure "order" is complete')
    simulation = simulation[order,]
  }
  maxCI = max(simulation,CI,na.rm = T)
  tmp = c(min(simulation),maxCI)
  #if(is.null(ylim)) ylim = ifelse(min(simulation)<maxCI, tmp, rev(tmp)) # ifelse only returns the 1st element
  if(is.null(ylim)){
    if(min(simulation)<maxCI){
      ylim = tmp
    }else{
      ylim = rev(tmp)
    }
  }  
  #matplot(simulation,type = 'l',col = c(rep('grey',B),'red'),ylab = 'CI(m|r)',axes=F,ylim=ylim)
  matplot(simulation,type = 'l',col = rep('grey',B),ylab = 'CI(m|r)',axes=F,ylim=ylim)
  
  x = which.max(CI)
#   if(is.null(order)){
#     #len = length(maxCIms[[2]][[x]])
#   }else{
#     #len = length(maxCIms[[2]][[x]])
#     x = which(order == names(CI)[x])
#     CI = CI[order]
#   }
#   text(x,maxCI+0.01,paste0('(',len,')'))
  
  if(!is.null(order)){
    if(is.null(names(CI))) stop('make sure "CI" is named using names in "order"')
    x = which(order == names(CI)[x])
   # CI = CI[order]
  }else if(length(CI) != nrow(simulation)) x = x+1

  axis(2)
  # customize x-axis
  if(is.null(order)){
    stages = row.names(simulation)
  }else{
    stages = order
  }
  axis(side=1,at=1:nrow(simulation),labels=stages,las = las)
  points(x,maxCI,col = 'red',pch = 16)
}


#' get Ic score
#' 
#' @description get Ic score see (ref)
#' 
#' @param counts a count matrix with unique loci row names:ID X column names:samples
#' @param sampleL a list of samples seperate by stages. The order of this output is determined by the order of this list 
#' @param genes a DNS loci IDs in counts
#' @param output select from 'Ic','PCCg', or 'PCCs'. Uses 'Ic' by default. 
#' 'PCCg' is the PCC between genes (numerator) and 'PCCs' is PCC between samples (denominator) 
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

getIc = function(counts,sampleL,genes,output = 'Ic'){
  subsetC = subset(counts, row.names(counts) %in% genes)
  subsetC = lapply(sampleL,function(x) subsetC[,x])
  
  # subsetC = lapply(maxCIms[[2]],function(x) subset(counts,row.names(counts) %in% x))
  # names(subsetC) = names(maxCIms[[2]])
  # sampleL = lapply(names(maxCIms[[2]]),function(x) as.character(subset(cli,group_xy == x)$Row.names))
  # names(sampleL) = names(maxCIms[[2]])
  # subsetC = lapply(names(sampleL),function(x) subsetC[[x]][,sampleL[[x]]])
  # sapply(subsetC,dim)
  # #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
  # # [1,]   75  142   46  119   15   46   12  145
  # # [2,]  168   23   72   29   19   11   56  120
  # 
  PCCg = lapply(subsetC,function(x) abs(cor(t(x))))
  for(i in 1:length(PCCg)) PCCg[[i]][lower.tri(PCCg[[i]],diag = T)] = NA
  PCCg = sapply(PCCg,function(x) mean(x,na.rm = T))
  #print(PCCg)
  
  PCCs = lapply(subsetC,function(x) cor(x))
  for(i in 1:length(PCCs)) PCCs[[i]][lower.tri(PCCs[[i]],diag = T)] = NA
  PCCs = sapply(PCCs,function(x) mean(x,na.rm = T))
  #print(PCCs)
  toplot = PCCg/PCCs
  names(toplot) = names(PCCg) = names(PCCs) = names(sampleL)
  if(output == 'Ic'){
    return(toplot)
  }else if(output == 'PCCg'){
    return(PCCg)
  }else if(output == 'PCCs'){
    return(PCCs)
  }
}

#' plot the Ic score for each stage
#' 
#' @description plot a line plot with Ic score for each stage
#' 
#' @param Ic a vector, names needed. if order is not assigned then plot by the order of this vector
#' @param las the direction of x-labels. default 0. labels are parallel (=0) or perpendicular(=2) to axis
#' @param order the ordered stages names
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

plotIc = function(Ic,las = 0,order = NULL){
  if(!is.null(order)){
    if(any(!order %in% names(Ic))) stop('make sure "Ic" is named using names in "order"')
    if(any(!names(Ic) %in% order)) warning('not every state in "Ic" is plotted, make sure "order" is complete')
    Ic = Ic[order]
  }
  matplot(Ic,type = 'l',ylab = 'Ic',axes=F)
  axis(2)
  stages = names(Ic)
  axis(side=1,at=1:length(Ic),labels=stages,las = las)
}

#' simulation of loci to calculate the Ic score
#' 
#' @description simulation of loci to calculate the Ic score
#' 
#' @param obs.x anumeric value, length of identified DNS
#' @param bk a vector of loci in the background
#' @param sampleL a list of samples seperate by stages. The order of this output is determined by the order of this list
#' @param counts 
#' @param B number of times to run this simulation, default 1000
#' @return a list of stats
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

# commented out 4/30/2019
# simulation_loci = function(x,bk,sampleL,counts,saveN,B = 1000){
#   obs.x = ifelse(length(x),x)
#   random = sapply(1:B, function(x) sample(bk,obs.x))
#   tmp= sapply(1:B, function(x) getIc(counts,sampleL,random[,x],output = 'Ic'))
#   return(tmp)
# }

simulation_Ic_loci = function(obs.x,bk,sampleL,counts,B = 1000){
  random = sapply(1:B, function(x) sample(bk,obs.x))
  tmp= sapply(1:B, function(x) getIc(counts,sampleL,random[,x],output = 'Ic'))
  row.names(tmp) = names(sampleL)
  return(tmp)
}


#' simulation of loci to calculate the Ic score
#' 
#' @description simulation of loci to calculate the Ic score
#' 
#' @param Ic score obtained in 
#' @param simulation a matrix state * number of simulated times, can be obtained from simulation_loci
#' @param B run time of simulation, default 1000 
#' @param order the ordered stages names
#' @return a line plot of simulated 
#' @export
#' @examples
#' 
#' 
#' @author Zhezhen Wang

plot_Ic_Simulation = function(Ic,simulation,las = 0,ylim = NULL,B=1000,order = NULL){
  if(!identical(row.names(simulation),names(Ic))) stop("Check if simulation and Ic is in the same order, make sure that row names of simulation is identical to names of Ic")
  toplot = cbind(simulation,Ic)
  if(!is.null(order)){
    if(any(!names(Ic) %in% order)) stop('make sure "Ic" is named using names in "order"')
    toplot = toplot[order,]
  }
  matplot(toplot,type = 'l',col = c(rep('grey',B),'red'),ylab = 'Ic',axes=F,ylim=ylim)
  axis(2)
  # customize x-axis
  stages = row.names(toplot)
  axis(side=1,at=1:length(Ic),labels=stages,las = las)
}

# plotdensity = function(tmp,obs){
#   #obs= obs[idx,idx]
#   plot(density(tmp))
#   abline(v = obs,col = 'red',lty = 2)
# }

#' simulation of samples to calculate the Ic score
#' 
#' @description simulation of samples to calculate the Ic score
#' 
#' @param Ic a vector, names needed. plot ordered by this order
#' @param bk
#' @param genes
#' @param B number of times to run this simulation, default 1000 
#' @param title default 'simulation of samples'
#' @return a numeric vector and a density plot
#' @export 
#' @examples
#' 
#' 
#' @author Zhezhen Wang

plot_simulation_sample = function(counts,nobs,obs,bk,genes,B = 1000,title = 'simulation of samples'){
  sampleL = lapply(1:B, function(x) sample(bk,nobs))
  tmp= sapply(1:B, function(x) getIc(counts,sampleL[x],genes,output = 'Ic'))
  p_v = length(tmp[tmp>obs])/B
  den = density(tmp)
  plot(den,main = title)
  abline(v = obs,col = 'red',lty = 2)
  x = max(den$x) - 0.2*diff(range(den$x))
  text(x,max(den$y),paste0('p-value=',p_v))
  return(tmp)
}

