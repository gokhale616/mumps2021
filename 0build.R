library(wrapr)
library(knitr)
## tinytex requires local install (fast, once per machine)
## see https://yihui.org/tinytex/ 
library(tinytex)
## %.>%
library(wrapr)

# this scripts produces a new directory to store processed data 
if(dir.create("./output") == TRUE) {
  message("'./output' already exists at root! Proceeeding to populate")
} else {
  dir.exists("./output")  
  message("'./output' created at root! Proceeeding to populate")
}


## convenience function: build path
mk_path <- function(x, dir='output/') paste0(dir,x)

## knit: Rnw to tex (returns file path)
## pdflatex: tex to pdf
(   knit(input = 'figs.Rnw', output=mk_path('figs.tex'))
  %.>% tinytex::pdflatex(
    file=., pdf_file=mk_path('figs.pdf')
  )
)
##





if(FALSE) {
  
  seq(0, 1, by = 0.1)
  
  10^seq(0, 1, by = 0.1)
  
  10^seq(-1,0, by = 0.1)
  
  10^seq(1, 2, by = 0.1)
  
}
  
  