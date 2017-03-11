#FINDSUB test Data
#testing of findsub function

#EDITED: 3.5.17



require(igraph) #load igraph package. need iGraph package
require(testthat)
folder<-"C:/Users/HPDM1/Documents/CanadaUofT/4thYear/BCB420/ekplektoR/R/"
#note: you will have to change folder and filename depending on changes to the function's filename
#and file location
filename<-paste(folder,"FINDSUBv1_2.R",sep="")

source(filename) #sources the script encoding FINDSUB function
HGNCsymb<-c("ABC1","DEF2","GHI4","JKL3","LMN5","OPQ6","RST7","UVW8","XYZ9")
VertScores<-c(10,2,4,6,7,8,2,5,9)

BigNetVertices<-data.frame(HGNC_symbol=HGNCsymb,Gene_Score=VertScores)
From<-c("ABC1","ABC1","DEF2","JKL3","DEF2","DEF2","JKL3","GHI4","JKL3","GHI4","GHI4","GHI4","LMN5",
        "XYZ9","LMN5","OPQ6","OPQ6","UVW8","OPQ6","RST7","UVW8","XYZ9")
To<-c("DEF2","JKL3","ABC1","ABC1","JKL3","GHI4","DEF2","DEF2","GHI4","JKL3","LMN5","XYZ9","GHI4",
      "GHI4","OPQ6","LMN5","UVW8","OPQ6","RST7","OPQ6","XYZ9","UVW8")
edgeID<-character(length=length(From))

for (i in 1:length(edgeID)) {
    edgeID[i]<-paste("BNE_",as.character(i),sep="")
}



#Influence vector created to assign influence levels to directed edges
Influence<-numeric(length(From))
Influence[1:length(Influence)]=5
Influence[(From=="GHI4")&(To=="LMN5")]=1
Influence[(From=="LMN5")&(To=="GHI4")]=2
Influence[(From=="GHI4")&(To=="XYZ9")]=1.5
Influence[(From=="XYZ9")&(To=="GHI4")]=2.5
Influence[(From=="OPQ6")&(To=="UVW8")]=2
Influence[(From=="UVW8")&(To=="OPQ6")]=3
Influence[(From=="UVW8")&(To=="XYZ9")]=3.5
Influence[(From=="XYZ9")&(To=="UVW8")]=3.5
#With a threshold influence of 3.5, LMN5 OPQ6 RST7 should form a network,
#ABC1,DEF2,GHI4,JKL3 should form a network
#When this is lowered to 3, LMN5,OPQ6,RST7,UVW8,and XYZ9 should form a network
#ABC1 DEF2 JKL3 GHI4 from same subnetwork
#Lower further to 2, and the 'discovered subnetwork' will include all vertices, but not all edges
#Lower influence threshold (delta) to 1, and every edge included

BigNetEdges<-data.frame(from=From,to=To,Influence=Influence,edgeID=edgeID) #edge data frame created for EGG
EGG<-graph_from_data_frame(BigNetEdges,directed=TRUE,vertices=BigNetVertices) #setup EGG graph
igraph::graph_attr(EGG, "EGGversion") <- "1.0"  #set EGG version attribute

as_data_frame(EGG,what="edges")
as_data_frame(EGG,what="vertices")

#[END]
