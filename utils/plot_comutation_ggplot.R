plot_comutation <- function(dat.integrate, 
                            out.dir="results/", 
                            logOR.limits = c(-3, 3),
                            pair.range = c(-1,0,5,15,25,50,100,200),
                            cache=T, 
                            Reordering=NULL)
{

        # Loading
        require(RColorBrewer)

        ########################
        # TESTING
        ########################

        # dat.binary <- dat.final %>% filter(DISEASE == "AML") %>% select(4:ncol(dat.final))
        # dat.binary[is.na(dat.binary)] <- 0
        # dat.binary[dat.binary>1] <- 1

        # dat.integrate <- as.matrix(dat.binary)
        # out.dir <- "results/"
        # cache <- F
        # Reordering <- NULL
        # logOR.limits <- c(-3,3)
        # pair.range <- c(-1,0,5,10,20,50,100,200)

        ## Reorder dat.integrate
        if (!is.null(Reordering))
                dat.integrate <- dat.integrate %>% mutate( variable = factor(variable, levels=Reordering$order))

        ########################
        # TESTING default values
        ########################
        # TODO: Verify if logOR.limits are integer values

        # TODO: Work on pair.range

        # TODO: Work on Reordering/Grouping


        ########################
        # CALCULATING OR/Pvals/Pairs
        ########################

        if (!cache)
        {
                system(paste0("mkdir -p ", out.dir,"/cache/"))
                logPInt <- sapply(1:ncol(dat.integrate), function(i) sapply(1:ncol(dat.integrate), function(j) {f<- try(fisher.test(dat.integrate[,i], dat.integrate[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
                odds <- sapply(1:ncol(dat.integrate), function(i) sapply(1:ncol(dat.integrate), function(j) {f<- try(fisher.test(table(dat.integrate[,i], dat.integrate[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
                pairs <- sapply(1:ncol(dat.integrate), function(i) colMeans(dat.integrate * dat.integrate[,i], na.rm=TRUE))

                # for caching purposes
                save(logPInt, file=paste0(out.dir,"/cache/logPInt.RData"))
                save(odds, file=paste0(out.dir,"/cache/odds.RData"))

        } else
        {
                logPInt <- get(load(paste0(out.dir,"/cache/logPInt.RData")))
                odds <- get(load(paste0(out.dir,"/cache/odds.RData")))
        }

        ########################
        # Formatting
        ########################

        # names
        colnames(pairs) <- rownames(pairs) <- colnames(dat.integrate)
        colnames(odds) <- rownames(odds) <- colnames(dat.integrate)
        colnames(logPInt) <- rownames(logPInt) <- colnames(dat.integrate)

        # formatting diag
        diag(logPInt) <- 0
        diag(odds) <- 1

        # Thresholding the odd ratios
        odds[odds<10^logOR.limits[1]] <- 10^(logOR.limits[1]-1)
        odds[odds>10^logOR.limits[2]] <- 10^(logOR.limits[2]+1)

        # putting non-significant OR at 1 for plotting
        odds[10^-abs(logPInt) > 0.05] = 1

        # logOdds
        logOdds=log10(odds)
        diag(logPInt) <- NA

        ###########################
        ## FORMATTING
        ###########################

        PInt.m <- reshape2::melt(logPInt) %>% tbl_df()
        pairs.m <- reshape2::melt(pairs) %>% tbl_df()
        logOdds.m <- reshape2::melt(logOdds) %>% tbl_df()

        if (is.null(Reordering))
        {
                order.variables <- levels(PInt.m$Var1)
        } else
        {
                order.variables <- Reordering
        }

        ## Working on range
        # i. Breaks for pairs 
        # TODO: improve pair range and labelling
        max.val <- pairs.m %>% .$value %>% max(., na.rm=T)
        if ((max.val *nrow(dat.integrate) + 100)> max(pair.range))
                pair.range <- c(pair.range, (max.val *nrow(dat.integrate) + 100))
        pair.labels <- formatC(pair.range[-1])

        # ii. Breaks for OR
        OR.breaks <- c((logOR.limits[1]-1):0-.Machine$double.eps,0:(logOR.limits[2]+1))
        # TODO: improve formatting
        OR.labels <- trimws(unique(formatC(10^OR.breaks, format = "g", digits=2)))

        ## Formatting data.frames
        pairs.upper <- pairs.m %>% 
                mutate(Var1=factor(Var1, levels=order.variables),
                       Var2=factor(Var2, levels=order.variables),
                       value=cut(value*nrow(dat.integrate), label=pair.labels, breaks=pair.range)) %>%
                filter(as.numeric(Var1) >= as.numeric(Var2))

        logOdds.lower <- logOdds.m %>% 
                mutate(Var1=factor(Var1, levels=order.variables),
                       Var2=factor(Var2, levels=order.variables),
                       value=cut(value,breaks = OR.breaks, label= OR.labels, include.lowest=TRUE))%>%
                filter(as.numeric(Var1) < as.numeric(Var2))

        PInt.lower <- PInt.m %>% 
                mutate(Var1=factor(Var1, levels=order.variables),
                       Var2=factor(Var2, levels=order.variables),
                       Pval=10^-abs(value)) %>%
                filter(as.numeric(Var1) < as.numeric(Var2))
        PInt.lower$adj.pval <- p.adjust(PInt.lower$Pval, method="fdr")

        ##################
        ## Color scale
        ##################

        # pairs colorscale
        colors.pairs <- brewer.pal(length(pair.labels), "Purples")
        names(colors.pairs) <- pair.labels

        # odds ratio colorscale
        colors.odds <- c(brewer.pal(length(OR.labels), "RdBu"),"white")
        names(colors.odds) <- c(setdiff(OR.labels," 10"), "  1") # TODO: rework with formatting
        colors.odds <- colors.odds[OR.labels] # Reorder

        # Add dummy values
        logOdds.lower <- bind_rows( logOdds.lower, data.frame(Var1=NA, Var2=NA, value=levels(logOdds.lower$value)))
        pairs.upper <- bind_rows( pairs.upper, data.frame(Var1=NA, Var2=NA, value=pairs.upper$value))

        ##################
        ## PLOT
        ##################

        ## Plot
        pp <- ggplot() + 
                geom_tile(data=pairs.upper, aes(x=Var1, y=Var2, fill=value), colour="white", show.legend=F) +
                geom_point(data=pairs.upper, aes(x=NA, y=NA, color=value)) + # dummy plots for legend purposes
                geom_tile(data=logOdds.lower, aes(x=Var1, y=Var2, fill=value), colour="white", show.legend=F) +
                geom_point(data=logOdds.lower, aes(x=NA, y=NA, alpha=value)) + # dummy plots for legend purposes
                geom_point(data=PInt.lower %>% filter(adj.pval <= 0.05), aes(x=Var1, y=Var2), shape=8, colour="salmon", show.legend=F) +
                geom_point(data=PInt.lower , aes(x=NA, y=NA, shape="Q < 0.05")) + # dummy plots for legend purposes
                scale_fill_manual(values=c(colors.pairs, colors.odds)) + 
                scale_x_discrete(limits=order.variables) +
                scale_y_discrete(limits=order.variables) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
                xlab("") + ylab("") +
                guides(color = guide_legend(title = "Frequency", order = 1, 
                                            override.aes = list(shape = 15, size = 5, color = colors.pairs) ),
                       alpha = guide_legend(title = "Odds ratio", order = 2, 
                                            override.aes = list(shape = 15, size = 5, color = colors.odds, alpha=1) ),
                       shape = guide_legend(title = "Significance", order = 3,
                                            override.aes = list(shape = 8, size = 3, color = "salmon"))) +
                theme(legend.key = element_rect(fill = "white") )


        return(pp)

}
