#' @import stringr
#' @import GenomicRanges
#' @import IRanges
#' @import seqinr
#' @title Makes a reference file for Salmon
#' @description This function creates a reference file that will be used by Salmon to create an index from a RepMask file.
#' The user can enter a RepMask file without deleting co-transcribed or overlapping repeats with the RepMask argument, or enter a RepMask file without co-transcribed but overlapping repeats with the RepMask.clean argument, or a file free of co-transcribed or overlapping repeats with the RepMask.ovlp.clean argument. When the file contains co-transcribed repeats, it must indicate rm.cotrans = T and when the file contains overlaps it must indicate overlapping = T.
#' @param overlapping Indicates whether the RepMask file contains overlapping repetitions (TRUE) or not (FALSE). When the RepMask file contains overlapping repetitions, the ovlp.res() function will be used to solve them and the resolution criteria must be indicated (higher score (HS), longer length (LE) or lower Kimura distances (LD))
#' @param rule A numerical vector respectively indicating the minimum percentage of identity, the percentage of the minimum length of a repeat with respect to the length of the transcript, and the length in base pairs of the repeats to be analyzed. Example: c(80, 80, 80)
#' @param trme transcriptome in fasta format
#' @param RepMask RepeatMasker output file. If rm.cotrans = F it is assumed that the file does not contain cotranscribed repeats. If overlapping = F it is assumed that the file does not contain overlapping.
#' @param rm.cotrnas logical vector indicating whether co-transcribed repeats should be removed
#' @param align .align file
#' @param over.res Indicates the method by which the repetition overlap will be resolved.
#' HS: higher score, bases are assigned to the element with the highest score
#' LS: longer element, bases are assigned to the longest element
#' LD: lower divergence, bases are assigned to the element with the least divergence.
#' in all cases both elements have the same characteristics, the bases are assigned to the first element.
#' @param anot annotation file in outfmt6 format. It is necessary when the option rm.cotrans = T
#' @param gff3 gff3 file. It is necessary when the option rm.cotrans = T
#' @param stranded logical vector indicating if the library is strand specific
#' @param cleanTEsProt logical vector indicating whether the search for TEs-related proteins should be carried out (e.g.
#' transposases, integrases, env, reverse transcriptase, etc.). We recommend that users use a curated annotations file,
#' in which these genes have been excluded; therefore the default option is F. When T is selected, a search is performed
#' against a database obtained from UniProt, so we recommend that the annotations file have this format for the subject
#' sequence id (e.g. "CO1A2_MOUSE"/"sp|Q01149|CO1A2_MOUSE"/"tr|H9GLU4|H9GLU4_ANOCA")
#' @param featureSum Returns statistics related to the characteristics of the transcripts. Requires a gff3 file. If TRUE, returns a list of the
#' @param outdir Output directory
#' @param by The column by which the repeats will be classified
#' @param best If best = T, the repetition with the highest percentage of coverage in the transcript will be selected (by default).
#' @param ignore.aln.pos The RepeatMasker alignments file may have discrepancies in the repeats positions with respect to the output file. If you selected over.res = "LD", then you can choose whether to take into account the positions of the alignment file or to take the average per repeats class (default).
#' @param threads Number of cores to use in the processing. By default threads = 1
#' @export
mk.reference <- function(RepMask,overlapping=F, by=c("namRep","classRep", "class", "supFam", "Fam"), rule=c(80,80,80), best=T, trme, threads=1, outdir, ...){
if(overlapping==T){
    RM <- ovlp.res(RepMask=RepMask, threads=threads, outdir=outdir,...)
  }else{
    RM <- RepMask
  }
message("reading lengths from fasta file ...")
  RM$width <- (RM$EPMQuer - RM$SPMQuer)
  f <- seqinr::read.fasta(trme)
  f_df <- data.frame(seqName=seqinr::getName(f), seqLength=seqinr::getLength(f))

  RM <- RM[RM$namSeq%in%f_df$seqName,]

cl <- parallel::makeCluster(threads)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, list("f_df", "RM"), envir=environment())
suppressWarnings(
RM$seqLength <- foreach::foreach(i=1:nrow(RM),.combine = rbind) %dopar% {
  f_df$seqLength[as.character(f_df$seqName)==as.character(RM$namSeq[i])]
}
)

message("building references ...")
refSeqs <- RM[RM$PersubM<(100-rule[1]) & RM$width>rule[3],c("namSeq", by, "width")]
suppressWarnings(
Ref.salmon <- foreach::foreach (i=unique(refSeqs[,1]), .combine = rbind ) %dopar% {
  Ref <- data.frame(
    seqNam=vector(),
    repNam=vector(),
    percentRep=numeric())
n <- 1
for(j in unique(refSeqs[refSeqs$namSeq==i,2])){
  perRep <- (sum(RM[RM$namSeq==i & RM[,by]== j, "width"])/mean(RM[RM$namSeq==i & RM[,by]== j, "seqLength"])*100)
    if(perRep>rule[3]){
      Ref[n,1] <- i
      Ref[n,2] <- j
      Ref[n,3] <- perRep
      n <- n+1
      }
    }
  Ref
}
)
parallel::stopCluster(cl)
foreach::registerDoSEQ()

message("writing files ...")
  if(best==T){
  Ref.salmon <- Ref.salmon[order(Ref.salmon$percentRep, decreasing = T),]
  Ref.salmon <- Ref.salmon[!duplicated(Ref.salmon$seqNam),]
}

Ref.salmon <-Ref.salmon[,c(1,2)]

write.table(Ref.salmon,paste0(outdir,"/references.csv"), col.names = F, row.names = F, quote = F, sep = ";")
decoys <- unique(f_df$seqName[f_df$seqName%!in%Ref.salmon$seqNam])
write.table(decoys,paste0(outdir,"/decoy.txt"), col.names = F, row.names = F, quote = F)

refseqs <- f[names(f)%in%unique(Ref.salmon$seqNam)]
decoyseqs <- f[names(f)%in%unique(decoys)]
seqinr::write.fasta(sequences = refseqs, names = names(refseqs),file.out = paste0(outdir,"/trmeSalmon.fasta"))
seqinr::write.fasta(sequences = decoyseqs, names = names(decoyseqs),file.out = paste0(outdir,"/trmeSalmon.fasta"),open = "a")

message(paste("The reference.csv and decoy.txt files are in", outdir, "directory"))

  Ref.salmon
}


