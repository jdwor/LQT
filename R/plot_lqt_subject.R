#' @title Plot subject-level summaries
#' @description This function constructs subject-level visualizations of lesion-based damage and disconnection
#' @param cfg a pre-made cfg structure (as list object).
#' @param subject either a string giving a subject ID, or an integer giving the subject index.
#' @param type a string specifying the desired plot type. Options are "parcel.damage", which gives a barplot of percent parcel damage;
#' "tract.discon", which gives a barplot of percent tract-level disconnection;
#' "parcel.discon", which gives a barplot and brain network plots showing parcel-level disconnection;
#' and "parcel.sspl", which gives a barplot and brain network plots showing parcel-level structural shortest path length increases.
#'
#' @importFrom utils read.csv
#' @importFrom ggplot2 ggplot scale_fill_manual theme_bw geom_point
#' @importFrom ggplot2 xlab ylab theme element_blank ggtitle geom_segment
#' @importFrom ggplot2 guides guide_legend scale_fill_manual coord_fixed
#' @importFrom ggplot2 theme_bw geom_bar coord_flip scale_alpha aes
#' @importFrom ggplot2 scale_size_identity scale_fill_identity
#' @importFrom ggplot2 scale_alpha_identity scale_alpha_continuous
#' @importFrom ggplot2 aes element_text scale_alpha_manual geom_text
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom igraph graph_from_adjacency_matrix E
#' @importFrom network as.matrix.network.edgelist network
#' @importFrom reshape2 melt
#'
#' @return A plot object visualizing either  parcel damage (type = 'parcel.damage'),
#' tract disconnection (type = 'tract.discon'),
#' inter/intra-parcel structural disconnection (type = 'parcel.discon'),
#' or changes in inter-parcel structural shortest path lengths (type = 'parcel.sspl').
#'
#' @export

plot_lqt_subject<-function(cfg, subject=1, type=NULL){

  if(is.null(type)){
    stop("Please specify desired plot 'type'. Options are 'parcel.damage', 'tract.discon',
         'parcel.discon', 'parcel.sspl', 'group.to.group', and 'parcel.to.parcel'.")
  }

  if(is.null(cfg$pat_id)){
    pat_ids = unlist(lapply(cfg,`[[`,1))
    ind=ifelse(is.character(subject),
               which(pat_ids==subject),subject)
    cfg=cfg[[ind]]
  }

  cols=c("#F3AE6D","#516888","#C9DACA","#FBE697","#802729",
         "#7DCCD3","#4E7147","#BE9C9D","#376597","#DBA662",
         "#C4878C","#F1E3B6","#469BEC","#6D882B","#C9FAFF",
         "#B23539","#FAB57C","#F7E790","#73652D","#E79498",
         "#514289","#EC8FA3","#FCBA65","#FAECCF","#8D7F99",
         "#8C9D57","#163343")
  ag6="#A12A19"

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

    ggplot(subdf,aes(x=Parcel,y=PercentDamage,fill=Color))+
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
  }else if(type=="tract.discon"){
    td.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Tract_Disconnection")
    trac.dis=read.csv(paste0(td.path,"/",cfg$pat_id,
                             "_percent_discon_tracts.csv"),header=T)[,-1]
    order=sort(trac.dis$Discon,decreasing=T,index.return=T)$ix
    subdf2=trac.dis[order,]; subdf2=subdf2[subdf2$Discon>cfg$tract_sdc_thresh,]
    subdf2$Tract=factor(subdf2$Tract,levels=rev(subdf2$Tract))
    subdf2$Pathway=as.character(subdf2$Pathway)
    cols2=c("#9FC2B2","#DFDED3","#A49A69","#3F5B66","#869144")

    ggplot(subdf2,aes(x=Tract,y=Discon,fill=Pathway))+
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
  }else if(type=="parcel.discon"){
    if(dir.exists(paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage"))){
      pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage")
      parc.dam=read.csv(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                               "_percent_damage.csv"),header=T)[,-1]
    }else{
      parc.dam=data.frame(Parcel=cfg$node_label,ParcelGroup=cfg$node_group,
                          PercentDamage=0)
    }
    parc.dam$Parcel=gsub(" \\(Right\\)","_R",parc.dam$Parcel)
    parc.dam$Parcel=gsub(" \\(Left\\)","_L",parc.dam$Parcel)
    parc.dam$Parcel=gsub("Lenticular nucleus, p",
                         "Lent_Nuc_P",parc.dam$Parcel)
    parc.dam$Parcel=gsub("RH_","R_",parc.dam$Parcel)
    parc.dam$Parcel=gsub("LH_","L_",parc.dam$Parcel)
    parc.dam$Color=cols[as.numeric(as.factor(parc.dam$ParcelGroup))]

    pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Disconnection")
    load(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                "_percent_parcel_mats.RData"))
    at.path=paste0(cfg$out_path,"/Atlas")
    load(paste0(at.path,'/atlas_',cfg$file_suffix,'_connectivity.RData'))

    at_con=connectivity; rm(connectivity)
    pat_sdc=at_con*pct_sdc_matrix/100

    parc.dam$ParcelDiscon=apply(pat_sdc,1,sum)/apply(at_con,1,sum)
    top10=sort(parc.dam$ParcelDiscon,decreasing=T,index.return=T)$ix[1:10]
    subdf=parc.dam[top10,]
    subdf$Parcel=factor(subdf$Parcel,levels=rev(subdf$Parcel))
    subdf$ParcelGroup=factor(subdf$ParcelGroup,
                             levels=unique(subdf$ParcelGroup))

    c=ggplot(subdf,aes(x=Parcel,y=ParcelDiscon*100,fill=Color))+
      geom_bar(stat="identity",alpha=.9,col="black")+coord_flip()+
      ggtitle("Parcel Disconnection")+
      geom_text(aes(x = Parcel, y=0.2, label=Parcel),
                hjust=0,size=3,fontface="bold")+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=subdf$ParcelGroup,
                          breaks=subdf$Color)+
      theme_bw()+ylab("Overall % Disconnection")+xlab(NULL)+
      theme(plot.title=element_text(size=12,
                                    face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_text(size=8),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9),
            legend.text =element_text(size=8),
            legend.position='n',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    ndata=data.frame(Parcel=parc.dam$Parcel,
                     X=cfg$parcel_coords[,1],
                     Y=cfg$parcel_coords[,2],
                     Z=cfg$parcel_coords[,3],
                     Comm=parc.dam$ParcelGroup,
                     Discon=parc.dam$ParcelDiscon)
    ndata$Size=5*(ndata$Discon+.01)^(1/6)
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
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_alpha_continuous(name="Disconnection",
                             breaks=c(.25,.50,.75,1.00),
                             labels=c("25%","50%","75%","100%"))+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=ndata$Comm,breaks=ndata$Color)+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
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
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_fill_identity()+scale_alpha(guide="none")+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
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
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_fill_identity()+scale_alpha(guide="none")+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    c+g2+g1+g3+
      plot_layout(byrow=F,nrow=2,guides="collect")+
      theme(plot.title=element_text(size=12,face='bold',
                                    colour = '#3B3B3B',
                                    hjust=.2,vjust=1),
            plot.subtitle = element_text(size=10.5,
                                         colour = '#3B3B3B',
                                         hjust=.14,vjust=1))
  }else if(type=="parcel.sspl"){
    if(dir.exists(paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage"))){
      pd.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_Damage")
      parc.dam=read.csv(paste0(pd.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                               "_percent_damage.csv"),header=T)[,-1]
    }else{
      parc.dam=data.frame(Parcel=cfg$node_label,ParcelGroup=cfg$node_group,
                          PercentDamage=0)
    }
    parc.dam$Parcel=gsub(" \\(Right\\)","_R",parc.dam$Parcel)
    parc.dam$Parcel=gsub(" \\(Left\\)","_L",parc.dam$Parcel)
    parc.dam$Parcel=gsub("Lenticular nucleus, p",
                         "Lent_Nuc_P",parc.dam$Parcel)
    parc.dam$Parcel=gsub("RH_","R_",parc.dam$Parcel)
    parc.dam$Parcel=gsub("LH_","L_",parc.dam$Parcel)
    parc.dam$Color=cols[as.numeric(as.factor(parc.dam$ParcelGroup))]

    ps.path=paste0(cfg$out_path,"/",cfg$pat_id,"/Parcel_SSPL")
    load(paste0(ps.path,"/",cfg$pat_id,"_",cfg$file_suffix,
                "_SSPL_matrices.RData"))
    perc_d_sspl=delta_sspl_matrix/(sspl_matrix-delta_sspl_matrix)
    perc_d_sspl[is.na(perc_d_sspl)]=0

    parc.dam$ParcelSSPL=apply(perc_d_sspl,1,mean)
    top10=sort(parc.dam$ParcelSSPL,decreasing=T,index.return=T)$ix[1:10]
    subdf=parc.dam[top10,]
    subdf$Parcel=factor(subdf$Parcel,levels=rev(subdf$Parcel))
    subdf$ParcelGroup=factor(subdf$ParcelGroup,
                             levels=unique(subdf$ParcelGroup))

    c=ggplot(subdf,aes(x=Parcel,y=ParcelSSPL*100,fill=Color))+
      geom_bar(stat="identity",alpha=.9,col="black")+coord_flip()+
      ggtitle("Parcel SSPL Increase")+
      geom_text(aes(x = Parcel, y=0.2, label=Parcel),
                hjust=0,size=3,fontface="bold")+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=subdf$ParcelGroup,
                          breaks=subdf$Color)+
      theme_bw()+ylab("Avg. % SSPL Increase")+xlab(NULL)+
      theme(plot.title=element_text(size=12,
                                    face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_text(size=8),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9),
            legend.text =element_text(size=8),
            legend.position='n',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    parc.sspl=apply(perc_d_sspl,1,mean)
    ndata=data.frame(Parcel=parc.dam$Parcel,
                     X=cfg$parcel_coords[,1],
                     Y=cfg$parcel_coords[,2],
                     Z=cfg$parcel_coords[,3],
                     Comm=parc.dam$ParcelGroup,
                     Damage=parc.sspl*100)
    ndata$Color=cols[as.numeric(as.factor(ndata$Comm))]
    xyz = ndata[,2:4]

    g_del = graph_from_adjacency_matrix(delta_sspl_matrix,mode="undirected",diag=F,weighted=T)
    g_pdel = graph_from_adjacency_matrix(perc_d_sspl,mode="undirected",diag=F,weighted=T)
    net = network(perc_d_sspl, directed=F,)
    edges = as.matrix.network.edgelist(net)
    edges = data.frame(xyz[edges[, 1], ], xyz[edges[, 2], ])
    names(edges) = c("x1", "y1", "z1", "x2", "y2", "z2")
    edges$Size=(E(g_del)$weight)/10
    edges$Damage=E(g_pdel)$weight*100

    ndata$Size=log(ndata$Damage+.5)*3
    g1=ggplot(ndata,aes(x=X,y=Y))+coord_fixed(ratio=1)+
      geom_segment(data=edges,aes(x=x1,y=y1,xend=x2,yend=y2,
                                  size=Size,alpha=Damage),
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_fill_identity("Parcel Group",guide="legend",
                          labels=ndata$Comm,breaks=ndata$Color)+
      scale_size_identity()+theme_bw()+
      scale_alpha_continuous(name="SSPL Change (%)",
                             range=c(0.1,1))+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
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
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_fill_identity()+
      scale_alpha_continuous(name="SSPL Change (%)",
                             range=c(0.1,1))+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
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
                   color=ag6)+
      geom_point(data=ndata,aes(fill=Color),size=2,
                 color="black",shape=21)+
      scale_fill_identity()+
      scale_alpha_continuous(name="SSPL Change (%)",
                             range=c(0.1,1))+
      scale_size_identity()+theme_bw()+
      theme(plot.title=element_text(size=12,face = 'bold',
                                    colour = '#3B3B3B'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=9,angle=0,vjust=.5),
            axis.title.x=element_text(size=9),
            legend.title =element_text(size=9,vjust=4),
            legend.text =element_text(size=8),
            legend.position='right',
            axis.ticks=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    c+g2+g1+g3+
      plot_layout(byrow=F,nrow=2,guides="collect")+
      theme(plot.title=element_text(size=12,face='bold',
                                    colour = '#3B3B3B',
                                    hjust=.2,vjust=1),
            plot.subtitle = element_text(size=10.5,
                                         colour = '#3B3B3B',
                                         hjust=.14,vjust=1))
  }else if(type=="group.to.group"){

  }else if(type=="parcel.to.parcel"){

  }
}
