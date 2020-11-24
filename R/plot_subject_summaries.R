#' @title Plot subject-level summaries
#' @description This function
#' @param cfg a pre-made cfg structure (as list object).
#'
#' @importFrom neurobase readnii writenii
#' @importFrom utils read.csv write.csv
#' @importFrom ggplot2 ggplot scale_fill_manual theme_bw geom_point
#' @importFrom ggplot2 xlab ylab theme element_blank ggtitle geom_segment
#' @importFrom ggplot2 guides guide_legend scale_fill_manual coord_fixed
#' @importFrom ggplot2 scale_identity theme_bw geom_bar coord_flip
#' @importFrom nationalparkcolors park_palette
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom igraph graph_from_adjacency_matrix E
#' @importFrom network as.matrix.network.edgelist network
#'
#' @return A data.frame giving the percent of streamlines in each tract that were disconnected by the lesions (column "Discon")
#' and the associated tract names (column "Tract").
#'
#' @export

plot_subject_summaries<-function(cfg, type=c("")){

  cols=c(park_palette("GeneralGrant")[c(2,3,4,1,8)],
         park_palette("CraterLake")[c(-c(4,6))],
         park_palette("Zion")[c(4,3,1,5,2)],
         park_palette("DeathValley"),
         park_palette("BlueRidgePkwy"))

  # Parcel damage plot
  if(type=="parcel.damage"){
    pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage")
    parc.dam=read.csv(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                             "_percent_damage.csv"),header=T)[,-1]
    parc.dam$Parcel=gsub(" \\(Right\\)","_R",parc.dam$Parcel)
    parc.dam$Parcel=gsub(" \\(Left\\)","_L",parc.dam$Parcel)
    parc.dam$Parcel=gsub("Lenticular nucleus, p",
                         "Lent_Nuc_P",parc.dam$Parcel)
    parc.dam$Parcel=gsub("RH_","R_",parc.dam$Parcel)
    parc.dam$Parcel=gsub("LH_","L_",parc.dam$Parcel)
    parc.dam$Color=cols[as.numeric(as.factor(parc.dam$ParcelGroup))]
    top10=sort(parc.dam$PercentDamage,decreasing=T,index.return=T)$ix[1:10]
    subdf=parc.dam[top10,]
    subdf$Parcel=factor(subdf$Parcel,levels=rev(subdf$Parcel))
    subdf$ParcelGroup=factor(subdf$ParcelGroup,
                             levels=unique(subdf$ParcelGroup))

    p1=ggplot(subdf,aes(x=Parcel,y=PercentDamage,fill=Color))+
      geom_bar(stat="identity",alpha=.9,col="black")+coord_flip()+
      ggtitle("Parcel Damage")+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=subdf$ParcelGroup,
                          breaks=subdf$Color)+
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
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    return(p1)

  }else if(type=="tract.discon"){
    td.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Tract_Disconnection")
    trac.dis=read.csv(paste0(td.path,"/",cfg$pat_id,
                             "_percent_discon_tracts.csv"),header=T)[,-1]
    order=sort(trac.dis$Discon,decreasing=T,index.return=T)$ix
    subdf2=trac.dis[order,]; subdf2=subdf2[subdf2$Discon>cfg$tract_sdc_thresh,]
    subdf2$Tract=factor(subdf2$Tract,levels=rev(subdf2$Tract))
    subdf2$Pathway=as.character(subdf2$Pathway)
    cols2=park_palette("Yosemite")

    p2=ggplot(subdf2,aes(x=Tract,y=Discon,fill=Pathway))+
      geom_bar(stat="identity",alpha=.9,col="black")+
      coord_flip()+
      ggtitle("Tract Disconnection")+
      scale_fill_manual(values=cols2)+
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
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    return(p2)
  }else if(type=="parcel.discon"){
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
    pat_sdc=at_con*pct_sdc_matrix/100

    communities=unique(parc.dam$ParcelGroup)
    comm_sdc=matrix(nrow=length(communities),
                    ncol=length(communities))
    for(i in 1:length(communities)){
      for(j in 1:i){
        comm_i=communities[i]
        comm_j=communities[j]
        inds_i=parc.dam$ParcelGroup==comm_i
        inds_j=parc.dam$ParcelGroup==comm_j
        disc_ij=sum(pat_sdc[inds_i,inds_j])/sum(at_con[inds_i,inds_j])
        if(is.na(disc_ij)){disc_ij=0}
        comm_sdc[i,j]=disc_ij
        comm_sdc[j,i]=disc_ij
      }
    }
    comm_sdc
    rownames(comm_sdc)=communities
    colnames(comm_sdc)=communities
    corr <- reshape2::melt(comm_sdc, na.rm = TRUE)
    colnames(corr) <- c("Var1", "Var2", "value")
    corr$value=corr$value*100
    corr$size=corr$value/10
    corr=corr[corr$value>0,]
    corr$color=park_palette("ArcticGates")[6]

    c=ggplot(corr, mapping=aes_string(x="Var1",y="Var2",fill="color"))+
      geom_point(color="black",shape=21,
                 aes_string(size="size",alpha="value"))+
      coord_fixed(ratio=1)+
      scale_size_identity()+guides(alpha = FALSE) +
      scale_fill_identity()+theme_bw()+xlab(NULL)+ylab(NULL)+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_text(size=7,angle=30,hjust=1),
            axis.text.y=element_text(size=7,angle=30,hjust=1),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    ndata=data.frame(Parcel=parc.dam$Parcel,
                     X=cfg$parcel_coords[,1],
                     Y=cfg$parcel_coords[,2],
                     Z=cfg$parcel_coords[,3],
                     Comm=parc.dam$ParcelGroup,
                     Damage=parc.dam$PercentDamage)
    ndata$Size=5*(ndata$Damage+.01)^(1/6)
    ndata$Color=cols[as.numeric(as.factor(ndata$Comm))]
    xyz = ndata[,2:4]

    at_con2=at_con; at_con2[pct_sdc_matrix==0]=0
    g_at = graph_from_adjacency_matrix(at_con2,mode="undirected",diag=F,weighted=T)
    g_pat = graph_from_adjacency_matrix(pct_sdc_matrix,mode="undirected",diag=F,weighted=T)
    net = network(pct_sdc_matrix, directed=F,)
    edges = as.matrix.network.edgelist(net)
    edges = data.frame(xyz[edges[, 1], ], xyz[edges[, 2], ])
    names(edges) = c("x1", "y1", "z1", "x2", "y2", "z2")
    edges$Size=(E(g_at)$weight)^(1/1.5)/200
    edges$Damage=E(g_pat)$weight/100

    g1=ggplot(ndata,aes(x=X,y=Y))+coord_fixed(ratio=1)+
      geom_segment(data=edges,aes(x=x1,y=y1,xend=x2,yend=y2,
                                  size=Size,alpha=Damage),
                   color=park_palette("ArcticGates")[6])+
      geom_point(data=ndata,aes(fill=Color),size=2.5,
                 color="black",shape=21)+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=ndata$Comm,breaks=ndata$Color)+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    g3=ggplot(ndata,aes(x=X,y=Z))+coord_fixed(ratio=1)+
      geom_segment(data=edges,aes(x=x1,y=z1,xend=x2,yend=z2,
                                  size=Size,alpha=Damage),
                   color=park_palette("ArcticGates")[6])+
      geom_point(data=ndata,aes(fill=Color),size=2.5,
                 color="black",shape=21)+
      scale_fill_identity()+scale_alpha(guide="none")+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    g2=ggplot(ndata,aes(x=Y,y=Z))+coord_fixed(ratio=1)+
      geom_segment(data=edges,aes(x=y1,y=z1,xend=y2,yend=z2,
                                  size=Size,alpha=Damage),
                   color=park_palette("ArcticGates")[6])+
      geom_point(data=ndata,aes(fill=Color),size=2.5,
                 color="black",shape=21)+
      scale_fill_identity()+scale_alpha(guide="none")+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    p3=c+g2+g1+g3+
      plot_layout(byrow=F,nrow=2,guides="collect")+
      plot_annotation(title="Parcel disconnection",
                      subtitle="Visualization of inter-parcel and intra/inter-community disconnectivity") &
      theme(plot.title=element_text(size=12,face='bold',
                                    colour = '#3B3B3B',
                                    hjust=.05,vjust=1),
            plot.subtitle = element_text(size=10.5,
                                         colour = '#3B3B3B',
                                         hjust=.14,vjust=1))
    return(p3)


    ## Also what is the 'name' vector in the conn files?

  }


}
