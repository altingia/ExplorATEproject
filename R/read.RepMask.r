#' @import readr
#' @title Read RepeatMasker output
#' @description Read a RepeatMasker output
#' @param file RepeatMasker output file
#' @details This function loads the RepeatMasker output file and rename columns as:
#'\describe{
#'\item{scoreSW}{Smith-Waterman score of the match, usually complexity adjusted. The SW scores are not always directly comparable. Sometimes the complexity adjustment has been turned off, and a variety of scoring-matrices are used.}
#'\item{PersubM}{ % substitutions in matching region compared to the consensus}
#'\item{PerBasDel}{ % of bases opposite a gap in the query sequence (deleted bp)}
#'\item{PerBasIns}{ % of bases opposite a gap in the repeat consensus (inserted bp)}
#'\item{namSeq}{ name of query sequence}
#'\item{SPMQuer}{ starting position of match in query sequence}
#'\item{EPMQuer}{ ending position of match in query sequence}
#'\item{NbasAfEQuer}{ no. of bases in query sequence past the ending position of match}
#'\item{st}{ match is with the Complement of the consensus sequence in the database}
#'\item{namRep}{ name of the matching interspersed repeat}
#'\item{classRep}{ the class of the repeat, in this case a DNA transposon fossil of the MER2 group (see below for list and references)}
#'\item{NbasCompRep}{ no. of bases in (complement of) the repeat consensus sequence prior to beginning of the match}
#'\item{SPMRepdb}{ starting position of match in database sequence (using top-strand numbering)}
#'\item{EPMRepdb}{ ending position of match in database sequence}
#'\item{nu}{num index}
#'\item{HSM}{higher-scoring match whose domain partly (<80%) includes the domain of this match.}
#'}
#' @export
read.RepMask <- function(file) {as.data.frame(readr::read_table2(file, skip = 3,
                                                                 col_names = c("scoreSW", "PersubM","PerBasDel","PerBasIns","namSeq","SPMQuer","EPMQuer","NbasAfEQuer", "st",
                                                                               "namRep","classRep","NbasCompRep","SPMRepdb","EPMRepdb","nu","HSM") ))
}





