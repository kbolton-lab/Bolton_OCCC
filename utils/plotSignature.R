toPyr <- function(nctds) {
  nctds.pyr <- nctds
  nctds.pyr[nctds=='A'] <- 'T'
  nctds.pyr[nctds=='T'] <- 'A'
  nctds.pyr[nctds=='G'] <- 'C'
  nctds.pyr[nctds=='C'] <- 'G'
  as.character(nctds.pyr)
}


plotSignature <-function(subs.df, mut.order, main='', plot = FALSE) {
  oldw <- getOption("warn")
  options(warn = -1)

  options(warn = oldw)
  
  subs.df$bb <- substr(subs.df$triplets,1,1)
  subs.df$ba <- substr(subs.df$triplets,3,3)
  subs.df$isPyr <- subs.df$REF %in% c('C', 'T')
  subs.df$candi <- ''
  subs.df$candi <- paste0(subs.df$bb,'[', subs.df$REF,'>',subs.df$ALT,']', subs.df$ba )
  subs.df$candi[!subs.df$isPyr] <- paste0(toPyr(subs.df$ba[!subs.df$isPyr]),
                                          '[', toPyr(subs.df$REF[!subs.df$isPyr]),'>',
                                          toPyr(subs.df$ALT[!subs.df$isPyr]),']', 
                                          toPyr(subs.df$bb[!subs.df$isPyr]))
  
  b <- table(subs.df$candi)
  b <- b[ as.character(mut.order)]		
  b[is.na(b)] <- 0	
  if (plot) {
    barplot(b, col=c(rep('deepskyblue',16), rep('gray19',16), rep('red2', 16), rep('darkgrey', 16), rep('olivedrab1', 16), rep('mistyrose',16)),
            main=paste(nrow(subs.df$candi),'', main), las=2, border=NA, cex.axis=0.5, cex.names=0.4) + theme_bw()
      
  } else {
    return(b)
  }
}
