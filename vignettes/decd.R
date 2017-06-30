## ---- eval=FALSE---------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("DEComplexDisease")

## ---- eval=FALSE---------------------------------------------------------
#  rm(list=ls())

## ----  eval=FALSE--------------------------------------------------------
#  >download.file("https://superb-sea2.dl.sourceforge.net/project/decd/decd_data.rda",
#                "decd_data.rda");
#  >load("decd_data.rda")
#  >exp.brca[1:5,1:5]
#           TCGA.D8.A27K.01A TCGA.BH.A0HU.01A TCGA.BH.A5IZ.01A TCGA.AR.A5QM.01A TCGA.BH.A1EV.01A
#  TSPAN6               6801             2139             3274             1268             1649
#  TNMD                    5               33                2               44                0
#  DPM1                 1772             2582             2895             1375             2073
#  SCYL3                4659             2115             1019             1490             2321
#  C1orf112             1203              884             1591              478              748
#  >head(cl.brca, 4)
#  TCGA.D8.A27K.01A TCGA.BH.A0HU.01A TCGA.BH.A5IZ.01A TCGA.AR.A5QM.01A
#                 1                1                1                1

## ---- eval=FALSE---------------------------------------------------------
#  hc=hclust(dist(t(exp.brca[, cl==0])))

## ---- eval=FALSE---------------------------------------------------------
#  plot(hc)
#  rect.hclust(hc, k= 4, border = 'red')
#  

## ---- eval=FALSE---------------------------------------------------------
#  group=cutree(hc, k=4)
#  exp=exp.brca[, !colnames(exp.brca) %in% names(group[ group != 2])]
#  cl=cl.brca[ colnames(exp) ]

## ---- eval=FALSE---------------------------------------------------------
#  p=pnbinom(x, size=1/disp, mu=mu, lower.tail = F)

## ---- eval=FALSE---------------------------------------------------------
#    z=(x-mu)/sd
#    p=pnorm(z,lower.tail=F)

## ---- eval=FALSE---------------------------------------------------------
#  deg=bi.deg(exp, cl, method="edger", cutoff=0.05, cores=4)

## ---- eval=FALSE, results='axis'-----------------------------------------
#  Plot(deg, ann=ann.er, show.genes=c("ESR1","FOXA1","GATA3","FOXC1"))

## ---- eval=FALSE---------------------------------------------------------
#  res.deg=deg.specific(deg, min.genes=50, min.patients=5, cores=4)
#  

## ---- eval=FALSE---------------------------------------------------------
#  res.deg.test=deg.specific(deg, test.patients=brca1.mutated.patients, min.genes=50,
#                            min.patients=8, cores=4)
#  

## ---- eval=FALSE, results='asis'-----------------------------------------
#  Plot(res.deg, ann=ann.er, show.genes=c("ESR1","FOXA1","GATA3","FOXC1"))

## ---- eval=FALSE, results='asis'-----------------------------------------
#  Plot(res.deg.test, ann=ann.er, show.genes=c("ESR1","FOXA1","GATA3","FOXC1"))

## ---- eval=FALSE---------------------------------------------------------
#  seed.mod1 = seed.module(deg, res.deg=res.specific, min.genes=100, min.patients=50,
#                          overlap=0.85, cores=4)

## ---- eval=FALSE---------------------------------------------------------
#  seed.mod2 = seed.module(deg, test.patients=brca1.mutated.patients, min.genes=100,
#                          min.patients=20, overlap=0.85,  cores=4)

## ---- eval=FALSE---------------------------------------------------------
#  Plot(seed.mod1, ann=er.ann, type="model", max.n=5)

## ---- eval=FALSE---------------------------------------------------------
#  cluster.mod1 <-  seed.module(seed.mod1, cores=4)

## ---- eval=FALSE---------------------------------------------------------
#  cluster.mod2 <- seed.module(seed.mod1, vote.seed=T, cores=4)

## ---- eval=FALSE---------------------------------------------------------
#  >sort(names(cluster.mod1), decreasing=T)
#    [1] "M99"             "M98"             "M97"             "M96"             "M95"
#    [6] "M94"             "M93"             "M92"             "M91"             "M90"
#   [11] "M9"              "M89"             "M88"             "M87"             "M86"
#   [16] "M85"             "M84"             "M83"             "M82"             "M81"
#   [21] "M80"             "M8"              "M79"             "M78"             "M77"
#   [26] "M76"             "M75"             "M74"             "M73"             "M72"
#   [31] "M71"             "M70"             "M7"              "M69"             "M68"
#   [36] "M67"             "M66"             "M65"             "M64"             "M63"
#   [41] "M62"             "M61"             "M60"             "M6"              "M59"
#   [46] "M58"             "M57"             "M56"             "M55"             "M54"
#   [51] "M53"             "M52"             "M51"             "M50"             "M5"
#   [56] "M49"             "M48"             "M47"             "M46"             "M45"
#   [61] "M44"             "M43"             "M42"             "M41"             "M40"
#   [66] "M4"              "M39"             "M38"             "M37"             "M36"
#   [71] "M35"             "M34"             "M33"             "M32"             "M31"
#   [76] "M30"             "M3"              "M29"             "M28"             "M27"
#   [81] "M26"             "M25"             "M24"             "M23"             "M22"
#   [86] "M21"             "M20"             "M2"              "M19"             "M18"
#   [91] "M17"             "M16"             "M15"             "M14"             "M130"
#   [96] "M13"             "M129"            "M128"            "M127"            "M126"
#  [101] "M125"            "M124"            "M123"            "M122"            "M121"
#  [106] "M120"            "M12"             "M119"            "M118"            "M117"
#  [111] "M116"            "M115"            "M114"            "M113"            "M112"
#  [116] "M111"            "M110"            "M11"             "M109"            "M108"
#  [121] "M107"            "M106"            "M105"            "M104"            "M103"
#  [126] "M102"            "M101"            "M100"            "M10"             "M1"
#  [131] "M0"              "decd.input"      "decd.clustering"))
#  >names(cluster.mod1[["decd.input"]])
#  [1] "genes"         "patients"      "overlap"       "deg"           "test.patients" "min.genes"
#  [7] "min.patients"  "vote.seed"     "model.method"
#  >names(cluster.mod1[["decd.clustering"]])
#  [1] "group"     "represent"
#  >names(cluster.mod1[["M1"]])
#  [1] "max.genes"      "max.patients"   "genes.removed"  "patients.added" "curve"
#  [6] "seed"           "model"

## ---- eval=FALSE---------------------------------------------------------
#  Plot(cluster.module, ann=er.ann, type="model", max.n=5)

## ---- eval=FALSE---------------------------------------------------------
#  brca1.module <- seed.module(deg[, brca1.mutated.patients], min.row=100, min.col=7,
#                             cutoff=0.8, cores=4)
#  Plot(brca1.module, ann=ann.er, max.n=5, type="model")
#  Plot(brca1.module, ann=ann.er, deg=deg, max.n=5, type="model")

## ---- eval=FALSE---------------------------------------------------------
#  module.overlap(cluster.mod1, max.n=5)

## ---- eval=FALSE---------------------------------------------------------
#  res.mod1 <- seed.module(deg[,er.pos],  min.genes=100, min.patients=50,cutoff=0.85, cores=4)
#  res.mod1 <- cluster.module(res.mod1)
#  #modules of ER+ samples
#  res.mod2 <- seed.module(deg[,er.neg],  min.genes=100, min.patients=50,cutoff=0.85, cores=4)
#  res.mod2 <- cluster.module(res.mod2)
#  #modules of ER- samples

## ---- eval=FALSE---------------------------------------------------------
#  module.compare(res.mod1, res.mod2, max.n1=10, max.n2=10)

## ---- eval=FALSE---------------------------------------------------------
#  >names(cluster.mod1[["M1"]][["curve"]])
#  [1] "no.gene"    "no.patient" "score"
#  >head(cluster.mod1[["M1"]][["curve"]][["no.gene"]])
#  [1] 6616 6503 6475 6365 6264 6225
#  >head(cluster.mod1[["M1"]][["curve"]][["no.patient"]])
#  [1] 10 11 12 13 14 15

## ---- eval=FALSE---------------------------------------------------------
#  module.curve(cluster.mod1, "M1")

## ---- eval=FALSE---------------------------------------------------------
#  x=c(100,300)
#  names(x)<-c("M1","M3")
#  new.cluster.mod1=module.modelling(cluster.mod1, keep.gene.num = x, method='slope.clustering',
#                                  cores=4)
#  #here, only "M1" and "M3" are modified
#  new.cluster.mod1=module.modelling(cluster.mod1, keep.gene.num = 150)
#  # here, all the modules are modified
#  module.curve(new.cluster.mod1, "M1")

## ---- eval=FALSE---------------------------------------------------------
#  module.screen(cluster.mod1, feature.patients=brca1.mutated.patients)
#  #search modules
#  
#  module.screen(seed.mod1, feature.patients=brca1.mutated.patients,
#                method="fisher.test")

## ---- eval=FALSE---------------------------------------------------------
#  genes=cluster.mod1[["M1"]][["genes.removed"]]
#  # split the genes into overlapped windows
#  nas<-as.character(seq(1, length(genes)-400, by=50));
#  
#  #functional annotation
#  library(clusterProfiler)
#  go.res=lapply(nas, function(x){
#  	ges=genes[x:(x+400)]
#  	eg = bitr(ges, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#  	ego <- enrichGO(gene          = unique(as.vector(eg$ENTREZID)),
#  		            OrgDb         = org.Hs.eg.db,
#  		            ont           = "BP",
#  		            pAdjustMethod = "BH",
#  		            pvalueCutoff  = 1,
#  		            qvalueCutoff  = 1,
#  		    		    readable      = TRUE)
#  	
#  
#  	result=ego@result;
#  	rr=as.vector(result$GeneRatio);
#  	rr=as.numeric(gsub("\\/\\d+","",rr, perl=T))
#  	names(rr)<-as.vector(result$Description);
#  	out=data.frame(rr=rr, p=as.vector(result$pvalue))
#  	return(out)
#  })
#  names(go.res)<-as.character(nas)
#  
#  # collect the GO terms
#  gos=vector()
#  for(x in nas){
#  	gos=append(gos, row.names(go.res[[x]]));
#  }
#  gos=unique(gos)
#  
#  #collect the gene number
#  res1=matrix(ncol=length(nas), nrow=length(gos));
#  colnames(res1)<-nas;
#  row.names(res1)<-gos;
#  for(x in nas){
#  	res1[,x]=as.vector(go.res[[x]][gos,]$rr);
#  }
#  
#  #collect the p-value
#  res2=matrix(ncol=length(nas), nrow=length(gos));
#  colnames(res2)<-nas;
#  row.names(res2)<-gos;
#  for(x in nas){
#  	res2[,x]=as.vector(go.res[[x]][gos,]$p);
#  }
#  
#  res1=apply(res1, c(1,2), function(x) if(is.na(x)) 0 else x)
#  res2=apply(res2, c(1,2), function(x) if(is.na(x)) 1 else x)
#  
#  #one example
#  test.go="response to estrogen";
#  
#  #location of the removed genes
#  location=nas +200;
#  
#  df=data.frame(location=res[test,], x=res1[test.go,]] p=res2[test.go,])
#  
#  library(ggplot2)
#  qplot(x, cc, data=df, color=-log10(p), geom=c("point", "smooth"),
#        xlab="Windows location",ylab="Genes with GO annotation")
#  + scale_colour_gradientn(colours=rainbow(3)[c(3,2,1)])
#  + annotate(geom="text", x=2000, y=17, label="GO:0043627:response to estrogen",
#             color="green",size=4)

