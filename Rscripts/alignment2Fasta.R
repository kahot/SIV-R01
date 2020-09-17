alignment2Fasta <- function(alignment, filename) {
    sink(filename)
    
    n <- length(rownames(alignment))
    for(i in seq(1, n)) {
        cat(paste0('>', rownames(alignment)[i]))
        cat('\n')
        the.sequence <- toString(unmasked(alignment)[[i]])
        cat(the.sequence)
        cat('\n')  
    }
    
    sink(NULL)
}