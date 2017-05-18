get.mutect2.data <- function(path='.', pattern='hg19_multianno.txt$') {
  files <- list.files(path=path, pattern=pattern, recursive=TRUE, full.names=TRUE)
  data <- data.frame()
  for(i in 1:length(files)) {
    tmp.data <- try(
      read.table(
        file = files[i],
        header = FALSE,
        as.is = TRUE,
        sep = '\t',
        quote = "\"",
        skip = 1
        ),
      silent = TRUE
      )

    # the "try" block will return a class of "try-error"
    if (class(tmp.data) == 'try-error') {
      next()
      }
    else {
      data <- rbind(tmp.data, data)
      }
    }
  colnames(data) <- ShlienLab.Core.SSM::get.mutect2.annotated.header()
  return(data)
  }
