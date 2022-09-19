require(GenVisR)
require(biomaRt)

filter_residue <- function (residueSeq) {
        if (any(residueSeq %in% c("OPAL", "OCHRE", "AMBER"))) {
                stopRes <- c("OPAL", "OCHRE", "AMBER")
                residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
        }

        return(residueSeq)
        
}

get_protein_info <- function (gene) {

        # TODO Recover good version
        ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= "hsapiens_gene_ensembl", ensemblRedirect = F)
        cds.infos <- getBM(attributes=c("coding", "cds_length","ensembl_transcript_id"), filters = "hgnc_symbol", values = gene, mart= ensembl)

        # Removing non available sequences
        cds.infos <- cds.infos[!is.na(cds.infos$cds_length),]
        # Removing sequences that are not of length /3
        cds.infos <- cds.infos[(cds.infos$cds_length%%3==0),]

        if (nrow(cds.infos)!=0)
        {
                # Keeping longest transcript
                longest.transcript.list <- cds.infos$ensembl_transcript_id[(cds.infos$cds_length==max(cds.infos$cds_length))]
                # If several keep the first one by default
                transcript.taken <- longest.transcript.list[1] 
                # Recover sequence
                codingSeq <- cds.infos$coding[cds.infos$ensembl_transcript_id==transcript.taken]
                # Transform sequence into AA
                residueSeq <- GenVisR:::lolliplot_DNAconv(codingSeq, to = "residue")
                # Filter out non AA codes
                residueSeq <- filter_residue(residueSeq)
                # Get length
                protein_length <- length(residueSeq)
                # Recover domains
                protein_domain <- getBM(attributes=c("interpro_description", "interpro_start","interpro_end"), filters=c("ensembl_transcript_id"), values=transcript.taken, mart=ensembl)
                if(all(is.na(protein_domain$interpro_description))) {
                        protein_domain <- data.frame(interpro_description=current.gene, interpro_start=1, interpro_end=protein_length)
                }
        } else
        {
                return(NULL)
        }

        return(list(protein_length=protein_length,protein_domain=protein_domain))
        
}
