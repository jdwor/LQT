#' @title Plot subject-level summaries
#' @description This function
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils read.csv write.csv
#' @importFrom ggplot2 ggplot scale_fill_manual theme_bw
#' @importFrom ggplot2 xlab ylab theme element_blank ggtitle
#' @importFrom ggplot2 guides guide_legend scale_fill_gradient2
#' @importFrom nationalparkcolors park_palette
#'
#' @return A data.frame giving the percent of streamlines in each tract that were disconnected by the lesions (column "Discon")
#' and the associated tract names (column "Tract").
#'
#' @export

plot_subject_summaries<-function(cfg, type=c("")){

  # Parcel damage plot
  if("parcel.damage"%in%type | "all"%in%type){
    pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage")
    parc.dam=read.csv(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                             "_percent_damage.csv"),header=T)[,-1]
    parc.dam$Parcel=gsub(" \\(Right\\)","_R",parc.dam$Parcel)
    parc.dam$Parcel=gsub(" \\(Left\\)","_L",parc.dam$Parcel)
    parc.dam$Parcel=gsub("Lenticular nucleus, p",
                         "Lent_Nuc_P",parc.dam$Parcel)
    parc.dam$Parcel=gsub("RH_","R_",parc.dam$Parcel)
    parc.dam$Parcel=gsub("LH_","L_",parc.dam$Parcel)
    top10=sort(parc.dam$PercentDamage,decreasing=T,index.return=T)$ix[1:10]
    subdf=parc.dam[top10,]
    subdf$Parcel=factor(subdf$Parcel,levels=rev(subdf$Parcel))
    subdf$ParcelGroup=as.character(subdf$ParcelGroup)
    cols=c(park_palette("GeneralGrant")[c(2,3,4,5,1,6,8)],
           park_palette("CraterLake")[1:3])

    ggplot(subdf,aes(x=Parcel,y=PercentDamage,fill=ParcelGroup))+
      geom_bar(stat="identity",alpha=.9,col="black")+coord_flip()+
      ggtitle("Parcel Damage")+
      scale_fill_manual(values=cols)+
      theme_bw()+ylab("% Damage")+xlab(NULL)+
      theme(plot.title=element_text(size=12,
                                    face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=7),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9),
            legend.text =element_text(size=8),
            legend.position='bottom',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

  }

  # Tract disconnection plot
  if("tract.discon"%in%type | "all"%in%type){
    td.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Tract_Disconnection")
    trac.dis=read.csv(paste0(td.path,"/",cfg$pat_id,
                             "_percent_discon_tracts.csv"),header=T)[,-1]
    order=sort(trac.dis$Discon,decreasing=T,index.return=T)$ix
    subdf=trac.dis[order,]; subdf=subdf[subdf$Discon>cfg$tract_sdc_thresh,]
    subdf$Tract=factor(subdf$Tract,levels=rev(subdf$Tract))
    subdf$Pathway=as.character(subdf$Pathway)
    cols=park_palette("Yosemite")

    ggplot(subdf,aes(x=Tract,y=Discon,fill=Pathway))+
      geom_bar(stat="identity",alpha=.9,col="black")+coord_flip()+
      ggtitle("Tract Disconnection")+
      scale_fill_manual(values=cols)+
      theme_bw()+ylab("% Disconnection")+xlab(NULL)+
      theme(plot.title=element_text(size=12,
                                    face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9),
            legend.text =element_text(size=8),
            legend.position='bottom',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  }

  # Parcel disconnection plot
  if("parcel.discon"%in%type | "all"%in%type){
    pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage")
    parc.dam=read.csv(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                             "_percent_damage.csv"),header=T)[,-1]
    pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Disconnection")
    load(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                "_percent_parcel_mats.RData"))
    at.path=paste0(cfg$out_path,"/Atlas")
    at.file=list.files(at.path,pattern="connectivity\\.RData$")
    load(paste0(at.path,"/",at.file))

    at_con=connectivity; rm(connectivity,atlas)
    at_con[pct_sdc_matrix==0]=0

    ord1=sort(parc.dam$PercentDamage,index.return=T)$ix
    mat1=pct_sdc_matrix[ord1,ord1]; names1=parc.dam$ParcelGroup[ord1]
    ord2=sort(parc.dam$ParcelGroup[ord1],index.return=T)$ix
    mat2=mat1[ord2,ord2]; names2=names1[ord2]
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

    names2lag=c(names2[-1],NA)
    diffsrow=nrow(mat2)-which(names2lag!=names2)+.5
    diffscol=which(names2lag!=names2)+.5
    rownames(mat1)=names1
    colnames(mat1)=NULL
    rownames(mat2)=names2
    colnames(mat2)=NULL
    corrplot(mat2/100,is.corr=F,tl.pos="l",
             method="color",tl.cex=.5,col=rev(col2(21)))
    segments(diffscol, rep(nrow(mat1)+.5,length(diffs)),
             diffscol, rep(0.5,length(diffs)), lwd=1)
    segments(rep(nrow(mat1)+.5,length(diffs)), diffsrow,
             rep(0.5,length(diffs)), diffsrow, lwd=1)

    g = graph_from_adjacency_matrix(at_con,mode="undirected",
                                    diag=F,weighted=T)
    g2 = graph_from_adjacency_matrix(pct_sdc_matrix,mode="undirected",
                                    diag=F,weighted=T)

    E(g)$color=park_palette("ArcticGates")[6]; cols=E(g)$color
    weights=E(g2)$weight/100
    #weights=rep(1,length(E(g2)$weight))
    editcols=function(x,cols,weights){
      return(adjustcolor(cols[x],alpha.f=weights[x]))
    }
    E(g)$color=unlist(lapply(1:length(cols),editcols,cols,weights))

    plot.igraph(g,vertex.label.cex=0.01,
                vertex.label.family="Helvetica",
                vertex.label.font=2,
                vertex.label.color="black",
                vertex.color=as.factor(parc.dam$ParcelGroup),
                vertex.size=6*(parc.dam$PercentDamage+.01)^(1/6),
                edge.width=(E(g)$weight)^(1/1.5)/100,
                xlim=c(-.9,.9),ylim=c(-.9,.9),
                layout=cfg$parcel_coords[,c(1,2)])

    ## Check on .edge files, something is off there
    ## Also what the hell is the 'name' vector in the conn files?

  }


}
