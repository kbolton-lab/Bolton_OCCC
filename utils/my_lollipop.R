library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

my_lollipop <- function(dat.genetics.unique, current.gene, diseases=NULL, protein_domain, protein_length, cutoff.hotspot=20, colors.mutation_type = NULL) {

    ### FIX: filter out the splice site ####
    dat.genetics.unique <- dat.genetics.unique %>% mutate(PROTEIN_POS = ifelse(PROTEIN_CHANGE == 0|is.na(PROTEIN_CHANGE), NA, PROTEIN_POS))

    ## 1. Protein
    protein_domain <- protein_domain %>% mutate(width = interpro_end  - interpro_start + 1,
                                                protein_length = protein_length) %>%
                                         group_by(interpro_description) %>% 
                                         mutate(height = 1 / (10*n())) %>% 
                                         ungroup()

    # TODO: FIX COLOR
    colourCount <- length(unique(protein_domain$interpro_description))
    p.protein <- ggplot(protein_domain %>% arrange(width)) + 
            geom_segment(aes(y=0,yend=0,x=1,xend=protein_length)) +  # 1. first plot the full protein
            geom_rect(aes(xmin=interpro_start,
                                  xmax=interpro_end,
                                  ymin=-height,
                                  ymax=height,
                                  fill=interpro_description),
                                  alpha= 0.3) + # 2. plot domains 
#             geom_text(data = protein_domain %>% distinct(interpro_description, .keep_all = TRUE) %>% mutate(interpro_mid = (interpro_end + interpro_start)/2), aes(x=interpro_mid,
#                                   y=height / 2 ,
#                                   label = interpro_description)) + 
            scale_fill_manual(values = getPalette(colourCount))+
            theme_void() +
            theme(legend.position = "bottom") +
            guides(fill=guide_legend(nrow=ifelse(colourCount < 3, 1, colourCount %/% 2))) +
            xlim(0, protein_length+10)

    ## 2. lollipop
    dat.genetics.unique.gene <- dat.genetics.unique %>% filter(gene == current.gene, !is.na(PROTEIN_POS))
    total_patients <- length(unique(dat.genetics.unique$ID))
    dat.score = dat.genetics.unique.gene %>%
        group_by(hit, PROTEIN_CHANGE, mutation_type) %>%
        summarise(
            gene = unique(gene),
            n_patients = length(unique(ID)),
            PROTEIN_POS = unique(PROTEIN_POS)
        )  %>%
#         mutate(score = ifelse(hit == '1', n_patients/total_patients, - n_patients/total_patients)) %>%
        mutate(score = ifelse(hit == '1', n_patients, - n_patients)) %>%
        ungroup() %>%
        dplyr::select(PROTEIN_CHANGE, gene, score, PROTEIN_POS, mutation_type)

    p.lollipop = ggplot(dat.score) + 
        geom_segment(
            aes(x=PROTEIN_POS, xend=PROTEIN_POS, y=0, yend=score),
            alpha=0.5, na.rm=T
        ) +
        geom_point(
            aes(x=PROTEIN_POS, y=score, color=mutation_type),
            alpha=0.7, na.rm=T
        ) + 
        ylab('Number of samples') + xlab("") + xlim(0, protein_length + 10) +
        theme(legend.position='top') + 
        scale_colour_manual(values=colors.mutation_type)
    
    p.final = ggarrange(
        p.lollipop, 
        p.protein,
        ncol=1, 
        nrow=2, 
        common.legend=F, 
        heights=c(4 * 4, 8), 
        align="v")

    p.final <- annotate_figure(p.final, top = text_grob(current.gene, color="black", face="bold", size=14))

    return(p.final)
}
