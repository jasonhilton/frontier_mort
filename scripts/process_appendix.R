
args <- commandArgs(trailingOnly = T)
prefix <- args[1]

if(is.na(prefix)){
  prefix <- ".."
}



rmarkdown::render(input = "paper/appendix.Rmd",
                  output_format = bookdown::pdf_document2(
                    toc=FALSE,
                    fig_caption=T,
                    keep_tex=T
                  ),
                  output_dir="paper",
                  params=list(prefix=prefix,
                              linear_ts="20191013_210236",
                              quad_ts="20191014_121020",
                              ind_ts="20191014_092729")
)
