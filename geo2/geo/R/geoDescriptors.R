require(jsonlite)
require(RCurl)
require(pmml)

#source("https://bioconductor.org/biocLite.R")
#biocLite(c('org.Hs.eg.db','GSEABase','GOstats'))
require(org.Hs.eg.db)
require(GSEABase)
require(GOstats)
require(Rgraphviz)
require(blockcluster)
require(Matrix)
require(vegan)

#1st step: read-in data

#dat2<- fromJSON('training_go_request.json')#datasetp
#dat3<- fromJSON('prediction_goM_request.json')#dataset_pred
#dat4<- fromJSON('training_go_requestH.json')
##dat5<- fromJSON('prediction_goM_requestH.json')#dataset_pred
#datg<- fromJSON('geoData.json')
#datBIO<- list(dataset=datg,parameters=dat2$parameters)
#datBIOh<- list(dataset=datg,parameters=dat4$parameters)
#sink("goData_request1.json")
#cat(toJSON(datBIOh))
#sink()

#datBIO<- fromJSON('goData_request.json')
#datBIOh<- fromJSON('goData_request1.json')

#dat1<- dat2$dataset#$dataEntry[,2]#dataEntries
#dat1p<- dat3$dataset#$dataEntry[,2]#dataEntries
#dat1m<- dat3$rawModel#?
#dat1i<- dat3$additionalInfo#
#dat4h<- dat4$dataset
##dat4param<- dat4$parameters

#save(datBIO,file='datBIO.rda')
#save(datBIO,file='datBIOh.rda')

#save(dat1m,file='dat1m.rda')
#save(dat1i,file='dat1i.rda')
#save(dat4h,file='dat4h.rda')

#x1<- generate.biclust.model(dat2$dataset,dat2$predictionFeature,dat2$parameters)
#dat3N<- list(dataset=list(datasetURI=character(),dataEntry=corona.dat$dataEntry),
#rawModel=x1$rawModel,additionalInfo=
#list(predictedFeatures=x1$additionalInfo$predictedFeatures,summaryStat=data.frame('mean')))
##base64(x1$rawModel,mode='raw')
#sink("prediction_goM_request.json")
#cat(toJSON(dat3N))
#sink()

#y1<- generate.hierar.model(dat4$dataset,dat4$predictionFeature,dat4$parameters)
#dat3NN<- list(dataset=list(datasetURI=character(),dataEntry=corona.dat$dataEntry),
#rawModel=y1$rawModel,additionalInfo=data.frame('mean'))#base64(x1$rawModel,mode='raw')
#sink("prediction_goM_requestH.json")
#cat(toJSON(dat3NN))
#sink()


##
#2 functions that include all steps for GO.descr
#1st: generate.descr.biclust, doing biclustering using 'blockcluster' library
#2nd: generate.descr.hierar, doing hierarchical clustering using 'vegan' library




generate.biclust.model<- function(dataset,parameters){
  #dataset:= list of 2 objects -
  #datasetURI:= character sring, code name of dataset
  #dataEntry:= data.frame with 2 columns,
  #1st:name of compound,2nd:data.frame with values (colnames are feature names)
  #rawModel:= ?
  #returns in additionalInfo the system names of properties (proteins)
  #additionalInfo:= FOR LM MODEL what are the features used for setting the model,
  #data frame with 2 columns: 1st  ModelCoef giving the dummy coefficient names produced for independent
  #features in the model, and 2nd RealFeatureNames which are ambit's feature names


  prot.data<- dataset$dataEntry$values
  #parameters<- dat2$parameters
  #prot.data<- dataset$dataEntry[,2]#d1$x.mat
  #prot.n<- colnames(prot.data)
  prot.n<- dataset$features$name
  ##parameters<- d1$model#key='UNIPROT',onto=c('GO','MF'),pvalCutoff=0.05,nclust=c(5,4),FUN=mean)

  #prot.data:=data frame with proteomics data, columns are dependent feature (1st) and proteins (remaining)
  #we only need the column protein names, a character vector of UNIPROT/SwissProt ids
  #key:= character sring for protein names id, default value is 'UNIPROT'
  #onto:= character vector showing the ontology and sub.ontologies used. Default value is c('GO','MF')
  #pvalCutoff:= numeric value for hypergeometric p-values cutoff. Default value is 0.05.
  #nclust:= numeric vector indicating the number of clusters for GOs(x axis) and proteins (y axis). Default value c(5,4).
  #FUN:= string, R function to summarize vector's groups

  #returns clustering classes for proteins to produce GO descriptors


  #steps:
  #a. produce binary go.table (kegg.table)
  #b. map protein ids (SwissProt/UNIPROT) to entrez ids & gene symbols
  #c. filter go.table using hypergeometric stats and pvalCutoff
  #d. cluster go.tab to find cluster ids for each protein
  #e. summarize prot.data to produce new descriptors using FUN


  #names1<- lapply(strsplit(prot.n,'/'),unlist)
  #names2<- unlist(lapply(names1,function(x)(return(x[length(x)]))))
  #rm(names1)

  #colnames(prot.data)<- names2
  prot.ids<- prot.n#names2

  unip.go<- select(org.Hs.eg.db, columns=c("GO"),
                   keys= prot.ids, keytype=parameters$key)#"UNIPROT")#g stands for go

  unip.go11<- as.data.frame(unip.go[1:dim(unip.go)[1],])
  go.idx<- duplicated(unip.go11[,1:2])
  unip.go11<- unip.go11[which(go.idx!=TRUE),1:2]

  go.tab1<- table(unip.go11$GO,unip.go11$UNIPROT)#matrix of 0-1s

  rm(unip.go,unip.go11,go.idx)

  #step a. construct binary GO x UNIPROT table
  rb.go.tab1<- go.tab1


  #step b. map protein ids (SwissProt/UNIPROT) to entrez ids & gene symbols
  #gids, gnam

  print("before unip.gid")
  unip.gid<- select(org.Hs.eg.db, columns=c("ENTREZID",'SYMBOL'),
                    keys= prot.ids, keytype=parameters$key)#"UNIPROT")
  gids<- unique(unip.gid$ENTREZID)
  gids<- gids[which(gids!='NA')]
  gnam<- unique(unip.gid$SYMBOL)
  gnam<- gnam[which(gnam!='NA')]


  #step c. filter out non-significant genes
  print("toTable")
  frame = toTable(org.Hs.egGO)# MAYBE LATER include in parameters
  #frame.om<- toTable(org.Hs.egOMIM)

  print("goFrameData")
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)

  goFrame=GOFrame(goframeData,organism="Homo sapiens")# MAYBE LATER include in parameters
  goAllFrame=GOAllFrame(goFrame)


  print("gsc")
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

  universe = Lkeys(org.Hs.egGO)
  my.ontol.choice<- c("MF","BP","CC")
  my.onto.res<- vector("list", length(my.ontol.choice))
  names(my.onto.res)<- my.ontol.choice

  print("for...")
  for(i in 1:length(my.ontol.choice)){
    params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                 geneSetCollection=gsc,
                                 geneIds = gids,
                                 universeGeneIds = universe,
                                 ontology = my.ontol.choice[i],
                                 pvalueCutoff = parameters$pvalCutoff,
                                 conditional = FALSE,
                                 testDirection = "over")
    over <- hyperGTest(params)
    #sum.go<- summary(over,htmlLinks=TRUE, pvalue=0.05)#perilamvanei ola ta GO me p-val<0.05
    sum.go <- GOstats::summary(over,htmlLinks=TRUE)
    #plot(goDag(Over))
    #dat.html<- htmlReport(Over)
    #description(Over)
    name1<- my.ontol.choice[i]
    my.onto.res[[i]]<- sum.go

  }

  imp.gos<- my.onto.res
  imp.gos<- get(parameters$onto[2],pos = imp.gos, envir = as.environment(imp.gos))#'MF'
  imp.gos.n<- paste(paste(parameters$onto[1],parameters$onto[2],sep=''),'ID',sep='')
  imp.gos.ids<- get(imp.gos.n, pos = imp.gos, envir = as.environment(imp.gos))#imp.gos$MF$GOMFID#!!! can I use onto[1] instead?


  #step d. cluster proteins
  print("step d")
  if(length(imp.gos.ids)!=0){
    sig.idx<- which(imp.gos.ids %in% rownames(rb.go.tab1))
    rb.go.tab1<- rb.go.tab1[sig.idx,]
    rm(sig.idx)
  }

  set.seed(1)
  newstrategy <- coclusterStrategy (nbtry =5, nbxem =10 , algo =' XCEMStrategy ',stopcriteria='Likelihood')
  clust.go<- cocluster(unclass(rb.go.tab1),'binary',nbcocluster=parameters$nclust)
  #NA KOITAKSW AN SUNEXIZEI KAI VGAZEI DIAFORETIKA CLUSTERS ANA RUN

  #clust.go.rowclass<- clust.go@rowclass
  clust.go.colclass<- clust.go@colclass
  #clust.go.ICL<- clust.go@ICLvalue

  rm(newstrategy, clust.go)

  #pred.descr<- function(dataset,rawModel,additionalInfo){
  #dataset:= list of 2 objects -
  #datasetURI:= character sring, code name of dataset
  #dataEntry:= data.frame with 2 columns,
  #1st:name of compound,2nd:data.frame with values (colnames are feature names)
  #rawModel:= numeric vector showing cluster memberships for proteins
  #returns in additionalInfo the system names of properties (proteins)
  #additionalInfo:= FOR LM MODEL what are the features used for setting the model,
  #data frame with 2 columns: 1st  ModelCoef giving the dummy coefficient names produced for independent
  #features in the model, and 2nd RealFeatureNames which are ambit's feature names


  #	d1<- read.in.json.for.pred(dataset,rawModel,additionalInfo)
  #	prot.data<- d1$x.mat
  clust.classes<- clust.go.colclass#d1$model#
  #adInfo<- d1$additionalInfo


  #prot.data:=data frame with proteomics data, columns are proteins

  #returns matrix with descriptors (one per class label in class.cprot.labels produced by biclustering for y axis (proteins))



  #summarize prot.data to produce new descriptors using FUN and clust.classes

  print("summarize prot.data big if!!!!")
  if(length(clust.go.colclass)!=0){
    class.cprot.labels<- clust.go.colclass


    #e. summarize tabel to produce descriprors using FUN
    class.go<- sort(unique(class.cprot.labels))
    if(class.go[1]==0){clust.names<- paste('clust.go.sep',class.go+1,sep='')}
    else{clust.names<- paste('clust.go.sep',class.go,sep='')}
    clust.descr<- matrix(0,dim(prot.data)[1],length(class.go))
    colnames(clust.descr)<- clust.names
    rownames(clust.descr)<- rownames(prot.data)

    for(i in 1:length(class.go)){
      idx.class<- which(class.cprot.labels==class.go[i])
      if(length(idx.class)>1){
        prot.f<- prot.data[,idx.class]
        clust.descr[1:dim(clust.descr)[1],i]<- apply(prot.f,1,
                                                     function(x)(get(parameters$FUN[1])(x,na.rm=T)))
        #FUN)#function(x)quantile(x,probs=0.75))
        #function(x){FUN(which(x==class.cprot.labels[i]))})
      }else{
        clust.descr[1:dim(clust.descr)[1],i]<- prot.data[,idx.class]
      }
    }

  }else{
    stop('Empty clusters specified.')
  }

  print("clust.descr")
  dataset$dataEntry$values <- as.data.frame(clust.descr)
  #feature <- dataset$features
  dataset$features<- as.data.frame(matrix(0,dim(clust.descr)[2],2))
  colnames(dataset$features)<- c('name','uri')
  dataset$features$uri <- colnames(clust.descr)

  #	for(i in 1:dim(clust.descr)[1]){
  #		w1<- data.frame(t(clust.descr[i,]))
  #		colnames(w1)<- clust.names
  #		if(i==1){p7.1<- list(unbox(w1))
  #		}else{
  #			p7.1[[i]]<- unbox(w1)
  #		}
  #	}
  #	p7.2<- list(predictions=p7.1)

  #return(p7.2)#as.data.frame(clust.descr))
  responseDataset <- list(responseDataset=dataset)
  return(responseDataset)#as.data.frame(clust.descr))
}

#try1<- generate.biclust.model(datBIO$dataset,datBIO$parameters)



generate.hierar.model<- function(dataset,parameters){
  #dataset:= list of 2 objects -
  #datasetURI:= character sring, code name of dataset
  #dataEntry:= data.frame with 2 columns,
  #1st:name of compound,2nd:data.frame with values (colnames are feature names)
  #rawModel:= ?
  #returns in additionalInfo the system names of properties (proteins)
  #additionalInfo:= FOR LM MODEL what are the features used for setting the model,
  #data frame with 2 columns: 1st  ModelCoef giving the dummy coefficient names produced for independent
  #features in the model, and 2nd RealFeatureNames which are ambit's feature names


  prot.data<- dataset$dataEntry$values
  prot.n<- dataset$features$name
  #parameters<- d1$model#key='UNIPROT',onto=c('GO','MF'),pvalCutoff=0.05,nclust=c(5,4),FUN=mean)

  #prot.data:=data frame with proteomics data, columns are dependent feature (1st) and proteins (remaining)
  #we only need the column protein names, a character vector of UNIPROT/SwissProt ids
  #key:= character sring for protein names id, default value is 'UNIPROT'
  #onto:= character vector showing the ontology and sub.ontologies used. Default value is c('GO','MF')
  #pvalCutoff:= numeric value for hypergeometric p-values cutoff. Default value is 0.05.
  #distMethod:= distance method, could be "manhattan", "euclidean", "canberra", "bray", "kulczynski",
  #"jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao"
  #or "mahalanobis"
  #hclustMethod:=hierarchical clustering method, the agglomeration method to be used.
  #This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete",
  #"average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  #nORh:= either numeric or character giving number of clusters or a function to define height respectively.
  #FUN:= string, R function to summarize vector's groups

  #returns clustering classes for proteins to produce GO descriptors


  #steps:
  #a. produce binary go.table (kegg.table)
  #b. map protein ids (SwissProt/UNIPROT) to entrez ids & gene symbols
  #c. filter go.table using hypergeometric stats and pvalCutoff
  #d. cluster go.tab to find cluster ids for each protein
  #e. summarize prot.data to produce new descriptors using FUN



  prot.ids<- prot.n#names2

  unip.go<- select(org.Hs.eg.db, columns=c("GO"),
                   keys= prot.ids, keytype=parameters$key)#"UNIPROT")#g stands for go

  unip.go11<- as.data.frame(unip.go[1:dim(unip.go)[1],])
  go.idx<- duplicated(unip.go11[,1:2])
  unip.go11<- unip.go11[which(go.idx!=TRUE),1:2]

  go.tab1<- table(unip.go11$GO,unip.go11$UNIPROT)#matrix of 0-1s

  rm(unip.go,unip.go11,go.idx)

  #step a. construct binary GO x UNIPROT table
  rb.go.tab1<- go.tab1


  #step b. map protein ids (SwissProt/UNIPROT) to entrez ids & gene symbols
  #gids, gnam
  unip.gid<- select(org.Hs.eg.db, columns=c("ENTREZID",'SYMBOL'),
                    keys= prot.ids, keytype=parameters$key)#"UNIPROT")
  gids<- unique(unip.gid$ENTREZID)
  gids<- gids[which(gids!='NA')]
  gnam<- unique(unip.gid$SYMBOL)
  gnam<- gnam[which(gnam!='NA')]


  #step c. filter out non-significant genes
  frame = toTable(org.Hs.egGO)# MAYBE LATER include in parameters
  #frame.om<- toTable(org.Hs.egOMIM)
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)

  goFrame=GOFrame(goframeData,organism="Homo sapiens")# MAYBE LATER include in parameters
  goAllFrame=GOAllFrame(goFrame)

  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

  universe = Lkeys(org.Hs.egGO)
  my.ontol.choice<- c("MF","BP","CC")
  my.onto.res<- vector("list", length(my.ontol.choice))
  names(my.onto.res)<- my.ontol.choice
  for(i in 1:length(my.ontol.choice)){
    params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                 geneSetCollection=gsc,
                                 geneIds = gids,
                                 universeGeneIds = universe,
                                 ontology = my.ontol.choice[i],
                                 pvalueCutoff = parameters$pvalCutoff,
                                 conditional = FALSE,
                                 testDirection = "over")

    Over <- hyperGTest(params)
    sum.go<- GOstats::summary(Over,htmlLinks=TRUE)#perilamvanei ola ta GO me p-val<0.05

    #plot(goDag(Over))
    #dat.html<- htmlReport(Over)

    description(Over)
    name1<- my.ontol.choice[i]
    my.onto.res[[i]]<- sum.go

  }


  imp.gos<- my.onto.res
  imp.gos<- get(parameters$onto[2],imp.gos)#'MF'
  imp.gos.n<- paste(paste(parameters$onto[1],parameters$onto[2],sep=''),'ID',sep='')
  imp.gos.ids<- get(imp.gos.n,imp.gos)#imp.gos$MF$GOMFID#!!! can I use onto[1] instead?


  #step d. cluster proteins
  if(length(imp.gos.ids)!=0){
    sig.idx<- which(imp.gos.ids %in% rownames(rb.go.tab1))
    rb.go.tab1<- rb.go.tab1[sig.idx,]
    rm(sig.idx)
  }

  set.seed(15)
  d <- vegdist(t(rb.go.tab1),method = parameters$distMethod,binary=TRUE)#"euclidean"#manhattan

  fit <- hclust(d, method=parameters$hclustMethod)#"ward.D2")
  q1<- parameters$nORh
  if(is.numeric(q1)==TRUE){
    class.cprot.Hlabels<- as.vector(cutree(fit,k=q1))
  }else{
    class.cprot.Hlabels<- as.vector(cutree(fit,h=get(q1)(fit$height)))#mean
  }

  rm(d,fit)

  clust.go.colclass<- class.cprot.Hlabels

  ####
  ####
  #		clust.go.colclass.ser<- serialize(clust.go.colclass,connection=NULL)
  #	#p1.pmml<- pmml(p1)
  #	#p1.pmml<- toString(p1.pmml)
  #	p1.n<- paste('clust.go.sep',1:length(unique(clust.go.colclass)),sep='')
  #
  #	m1.ser.list<- list(rawModel=clust.go.colclass.ser,
  #	predictedFeatures=p1.n,
  #	pmmlModel=NULL,independentFeatures=prot.n,additionalInfo=list(predictedFeatures=p1.n))
  #    return(m1.ser.list)

  clust.classes<- clust.go.colclass#d1$model#
  #adInfo<- d1$additionalInfo


  #prot.data:=data frame with proteomics data, columns are proteins

  #returns matrix with descriptors (one per class label in class.cprot.labels produced by biclustering for y axis (proteins))



  #summarize prot.data to produce new descriptors using FUN and clust.classes


  if(length(clust.go.colclass)!=0){
    class.cprot.labels<- clust.go.colclass


    #e. summarize tabel to produce descriprors using FUN
    class.go<- sort(unique(class.cprot.labels))
    if(class.go[1]==0){clust.names<- paste('clust.go.sep',class.go+1,sep='')}
    else{clust.names<- paste('clust.go.sep',class.go,sep='')}
    clust.descr<- matrix(0,dim(prot.data)[1],length(class.go))
    colnames(clust.descr)<- clust.names
    rownames(clust.descr)<- rownames(prot.data)

    for(i in 1:length(class.go)){
      idx.class<- which(class.cprot.labels==class.go[i])
      if(length(idx.class)>1){
        prot.f<- prot.data[,idx.class]
        clust.descr[1:dim(clust.descr)[1],i]<- apply(prot.f,1,
                                                     function(x)(get(parameters$FUN[1])(x,na.rm=T)))
        #FUN)#function(x)quantile(x,probs=0.75))
        #function(x){FUN(which(x==class.cprot.labels[i]))})
      }else{
        clust.descr[1:dim(clust.descr)[1],i]<- prot.data[,idx.class]
      }
    }

  }else{
    stop('Empty clusters specified.')
  }

  dataset$dataEntry$values <- as.data.frame(clust.descr)
  #feature <- dataset$features
  dataset$features<- as.data.frame(matrix(0,dim(clust.descr)[2],2))
  colnames(dataset$features)<- c('name','uri')
  dataset$features$uri <- colnames(clust.descr)

  #	for(i in 1:dim(clust.descr)[1]){
  #		w1<- data.frame(t(clust.descr[i,]))
  #		colnames(w1)<- clust.names
  #		if(i==1){p7.1<- list(unbox(w1))
  #		}else{
  #			p7.1[[i]]<- unbox(w1)
  #		}
  #	}
  #	p7.2<- list(predictions=p7.1)

  #return(p7.2)#as.data.frame(clust.descr))
  responseDataset <- list(responseDataset=dataset)
  return(responseDataset)#as.data.frame(clust.descr))

}
#generate.hierar.model(dat4$dataset,dat4$predictionFeature,dat4$parameters)
#try2<- generate.hierar.model(datBIOh$dataset,datBIOh$parameters)

read.in.json.for.pred<- function(dataset,rawModel,additionalInfo){
  #dataset:= list of 2 objects -
  #datasetURI:= character sring, code name of dataset
  #dataEntry:= data.frame with 2 columns,
  #1st:name of compound,2nd:data.frame with values (colnames are feature names)
  #rawModel:= serialized raw model for prediction - nly contains parameters for clustering
  #additionalInfo:= return FUN used for summary in prediction
  #FOR LM MODEL what are the features used for setting the model,
  #data frame with 2 columns: 1st  ModelCoef giving the dummy coefficient names produced for independent
  #features in the model, and 2nd RealFeatureNames which are ambit's feature names

  dat1.t<- dataset$dataEntry[,2]#$featureValues # data table
  dat1.n<- colnames(dat1.t)#dat1.n<- unlist(strsplit(colnames(dat1.t),'/'))
  names1<- lapply(strsplit(dat1.n,'/'),unlist)
  names2<- unlist(lapply(names1,function(x)(return(x[length(x)]))))
  rm(names1)

  colnames(dat1.t)<- names2

  dat1.m<- rawModel
  dat1.m<- base64Decode(dat1.m,'raw')
  dat1.model<- unserialize(dat1.m)

  return(list(x.mat=dat1.t,model=dat1.model,additionalInfo=additionalInfo))#data.frame(dat1.n)))

}

#read.in.json.for.pred(dat3$dataset,dat3$rawModel,dat3$additionalInfo)




pred.descr<- function(dataset,rawModel,additionalInfo){
  #dataset:= list of 2 objects -
  #datasetURI:= character sring, code name of dataset
  #dataEntry:= data.frame with 2 columns,
  #1st:name of compound,2nd:data.frame with values (colnames are feature names)
  #rawModel:= numeric vector showing cluster memberships for proteins
  #returns in additionalInfo the system names of properties (proteins)
  #additionalInfo:= FOR LM MODEL what are the features used for setting the model,
  #data frame with 2 columns: 1st  ModelCoef giving the dummy coefficient names produced for independent
  #features in the model, and 2nd RealFeatureNames which are ambit's feature names


  d1<- read.in.json.for.pred(dataset,rawModel,additionalInfo)
  prot.data<- d1$x.mat
  clust.classes<- d1$model#
  adInfo<- d1$additionalInfo


  #prot.data:=data frame with proteomics data, columns are proteins

  #returns matrix with descriptors (one per class label in class.cprot.labels produced by biclustering for y axis (proteins))



  #summarize prot.data to produce new descriptors using FUN and clust.classes


  if(length(clust.classes)!=0){
    class.cprot.labels<- clust.classes


    #e. summarize tabel to produce descriprors using FUN
    class.go<- sort(unique(class.cprot.labels))
    if(class.go[1]==0){clust.names<- paste('clust.go.sep',class.go+1,sep='')}
    else{clust.names<- paste('clust.go.sep',class.go,sep='')}
    clust.descr<- matrix(0,dim(prot.data)[1],length(class.go))
    colnames(clust.descr)<- clust.names
    rownames(clust.descr)<- rownames(prot.data)

    for(i in 1:length(class.go)){
      idx.class<- which(class.cprot.labels==class.go[i])
      if(length(idx.class)>1){
        prot.f<- prot.data[,idx.class]
        clust.descr[1:dim(clust.descr)[1],i]<- apply(prot.f,1,
                                                     function(x)(get(parameters$FUN[1])(x,na.rm=T)))
        #FUN)#function(x)quantile(x,probs=0.75))
        #function(x){FUN(which(x==class.cprot.labels[i]))})
      }else{
        clust.descr[1:dim(clust.descr)[1],i]<- prot.data[,idx.class]
      }
    }

  }else{
    stop('Empty clusters specified.')
  }

  for(i in 1:dim(clust.descr)[1]){
    w1<- data.frame(t(clust.descr[i,]))
    colnames(w1)<- clust.names
    if(i==1){p7.1<- list(unbox(w1))
    }else{
      p7.1[[i]]<- unbox(w1)
    }
  }
  p7.2<- list(predictions=p7.1)

  return(p7.2)#as.data.frame(clust.descr))
}



#new.descr<- pred.descr(dat3$dataset,dat3$rawModel,dat3$additionalInfo)





#package.skeleton(list=c('generate.biclust.model','generate.hierar.model','read.in.json.for.pred','pred.descr','datBIO','datBIOh'),name='GOdescrPred')
