library('maftools')

#get maf files from mutation dataframe

mdf_to_maf <- function(mut_df, clin_df) {
  
  indel.pass.maf <- mut_df %>%
    filter(variant_type != 'SNV') %>%
    dplyr::rename('Hugo_Symbol' = 'VAG_GENE', 'Chromosome' = 'CHR', 'Tumor_Sample_Barcode' = 'TARGET_NAME', 'Protein_Change' = 'VAG_PROTEIN_CHANGE') %>%
    mutate(Reference_Allele = case_when(VAG_VT == 'Del' ~ substr(REF, start = 2, stop = 1000000L),
                                        VAG_VT == 'Ins' ~ '-',
                                        TRUE ~ REF),
           Tumor_Seq_Allele2 = case_when(VAG_VT == 'Del' ~ '-',
                                         VAG_VT == 'Ins' ~ substr(ALT, start = 2, stop = 1000000L),
                                         TRUE ~ ALT),
           Start_Position = case_when(VAG_VT == 'Del' ~ START + 1L,
                                      VAG_VT == 'Ins' ~ START,
                                      TRUE ~ START),
           End_Position = case_when(VAG_VT == 'Del' ~ END,
                                    VAG_VT == 'Ins' ~ START + 1L,
                                    TRUE ~ END),
           Variant_Classification = case_when(VAG_EFFECT == 'frameshift_variant' & VAG_VT == 'Del' ~ 'Frame_Shift_Del',
                                              VAG_EFFECT == 'frameshift_variant' & VAG_VT == 'Ins' ~ 'Frame_Shift_Ins',
                                              VAG_EFFECT == 'splice_site_variant' ~ 'Splice_Site',
                                              VAG_EFFECT == '5_prime_UTR_variant'~ "5'UTR",
                                              VAG_EFFECT == 'inframe_codon_loss' ~ 'In_Frame_Del',
                                              VAG_EFFECT == 'inframe_codon_gain' ~ 'In_Frame_Ins',
                                              VAG_EFFECT =='inframe_variant' ~ 'In_Frame_Del',
                                              VAG_EFFECT =='stop_gained' ~ 'Nonsense_Mutation',  #to check
                                              VAG_EFFECT == 'complex_change_in_transcript' ~ 'In_Frame_Del', 
                                              VAG_EFFECT == '2KB_upstream_variant' ~ "Translation_Start_Site",
                                              VAG_EFFECT == 'frameshift_variant' & VAG_VT == 'Complex' ~ 'Frame_Shift_Del'), #?
           Variant_Type = case_when(VAG_VT == 'Del' ~ 'DEL',
                                    VAG_VT == 'Ins' ~ 'INS',
                                    TRUE ~ NA_character_)) %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type, Tumor_Sample_Barcode, TARGET_VAF_MEAN, Protein_Change)
  
  #Format SNV data
  #MAF variant classification : Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR (See Notes Section #1) , Intron, RNA, Targeted_Region
  # Should we agregate consecutive substitutions ?
  
  snv.pass.maf <- mut_df %>%
    filter(variant_type == 'SNV') %>%
    dplyr::rename('Hugo_Symbol' = 'VAG_GENE', 'Chromosome' = 'CHR', 'Reference_Allele'= 'REF', 'Tumor_Seq_Allele2' = 'ALT', 'Start_Position' = 'START', 'End_Position' = 'END','Tumor_Sample_Barcode' = 'TARGET_NAME', 'Protein_Change' = 'VAG_PROTEIN_CHANGE') %>%
    mutate(Variant_Classification = case_when(VAG_EFFECT == 'stop_gained' ~ 'Nonsense_Mutation',
                                              VAG_EFFECT == 'non_synonymous_codon' ~ 'Missense_Mutation',
                                              VAG_EFFECT == 'splice_site_variant' ~ 'Splice_Site',
                                              VAG_EFFECT == '5_prime_UTR_variant'~ "5'UTR",
                                              VAG_EFFECT == 'stop_retained_variant' ~ 'Nonstop_Mutation',
                                              VAG_EFFECT == 'initiator_codon_change' ~ 'Translation_Start_Site', 
                                              VAG_EFFECT =='synonymous_codon' ~ 'Silent',
                                              VAG_EFFECT == 'extended_intronic_splice_region_variant' ~ 'Splice_Site', 
                                              VAG_EFFECT == '2KB_upstream_variant' ~ "Translation_Start_Site"), 
           Variant_Type = 'SNP') %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type, Tumor_Sample_Barcode, TARGET_VAF_MEAN, Protein_Change)
  
  #merge indel and snv data
  pass.maf <- rbind(indel.pass.maf, snv.pass.maf)
  
  laml.clin <- clin_df %>%
    dplyr::rename('Tumor_Sample_Barcode' = 'leukgen_sample_id')
  
  #get maf df
  laml = read.maf(maf = pass.maf, clinicalData = laml.clin, removeDuplicatedVariants = FALSE)
  
  return(list(maf = pass.maf, clin = laml.clin, laml = laml))
  
}