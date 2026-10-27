args <- commandArgs(trailingOnly=TRUE)
for (p in args) {
  invisible(parse(file=p))
  cat('SYNTAX_CHECK PASS', p, '\n')
}
