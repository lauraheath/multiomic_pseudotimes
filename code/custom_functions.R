library(uuid)
library(tidyverse) # tidy
library(magrittr) # tidy %<>%
library(monocle) # pseudotime
# authors : JK, RB
list_to_rmd <- function(
    lst,                      # Input list. May contain ggplot objects, data.frames and strings
    outfile='output.html',    # Output filename
    title="Output",           # Title of RMarkdown document
    use.plotly=F,             # Output ggplots with plotly?
    keep.tmp.source=F,         # Keep the .RMD script and .RDS file?
    fig.width=7, fig.height = 5 # RMD parameters
) {
  
  # temporary files
  file.rmd <- paste0(tools::file_path_sans_ext(outfile),'.RMD')
  file.rds <- paste0(tools::file_path_sans_ext(outfile),'.rds')
  
  #### helper functions
  out <- function(str){writeLines(str, h)}
  # out <- function(str)cat(str, "\n")
  
  writechunk <- function(code, params='') {
    if (nchar(params)>0) params=paste0(' ', params)
    out(sprintf("```{r fig.width=%f,fig.height=%f,%s}", fig.width, fig.height, params))
    out(code)
    out("```\n")
  }
  
  #### SAVE LIST TO TMP FILE
  save(lst, file=file.rds)
  
  #### BUILD RMD 
  
  # initialize output
  h <- file(file.rmd, open='wt')
  
  #### markdown header and first heading
  out(glue::glue('
---
title: {title}
output:
  html_document:
    toc: true
    toc_float: TRUE
    number_sections: true
---

    '))
  
  ## global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")
  
  ## chunk that loads libraries
  if (!use.plotly) {
    writechunk('# load libraries\nlibrary(tidyverse)\n')
  } else  {
    writechunk('# load libraries\nlibrary(tidyverse)\nlibrary("plotly")')
  }
  ## chunk that loads data
  writechunk(glue::glue('# load data\nload("{file.rds}")\n'))
  
  
  ## go through list
  for (i in 1:length(lst)) {
    # add name as header
    out(sprintf("# %s", names(lst)[i]))
    # check what's inside
    if ("data.frame" %in% class(lst[[i]])) {
      # DATA FRAME
      # write out datatable
      writechunk(glue::glue('
# extract result table
df<-lst[[{i}]]
# output
DT::datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 20, lengthMenu = c(10*(2^(0:3)), nrow(df)), autoWidth = TRUE, width = 1200, dom = "Bitlrp", buttons = c("copy", "csv", "excel", "pdf", "print")), class = "cell-border stripe", extensions = "Buttons")  %>% DT::formatStyle(columns = c(1:ncol(df)), fontSize = "80%", target= "row", lineHeight="80%")'),
                 params = "results='asis'")
      
    } else if ("ggplot" %in% class(lst[[i]])) {
      # GGPLOT
      # plot
      if (!use.plotly) {
        writechunk( glue::glue("lst[[{i}]]"))
      } else {
        writechunk( glue::glue("
plotlist = list(lst[[{i}]] %>% ggplotly())
htmltools::tagList(setNames(plotlist, NULL))
                             "), params='results="show"')
      }
      
    } else if (class(lst[[i]]) %in% "character") {
      # STRING
      # just output
      out(lst[[i]])
      out("") # need a newline
      
    } else {
      stop(sprintf("Don't know what to do with list entry of type: %s", paste0(class(lst[[i]]), collapse = ", ")))
    }
  }
  
  # close rmd file
  close(h)
  
  
  #### KNIT
  rmarkdown::render(file.rmd)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(file.rmd),'.html'), outfile)
  # clean up?
  if (!keep.tmp.source) {
    file.remove(file.rmd)
    file.remove(file.rds)
  }
  
}

get_plot_for_all_outcomes <- function(hsmm, # monocle object
                                      outcomes # outcomes to overlay on hsmm
                                      ){
  
  # overall state designations
  p <- plot_cell_trajectory(hsmm, color_by = "State",
                            show_branch_points=F,
                            use_color_gradient = F,cell_size = 1.5) + 
    ggplot2::scale_color_viridis_d() +
    ggplot2::labs(color="State")
  ### empty list of plots
  plot_list <- list(); k<-1; plot_list[[k]] <- p; k <- k+1
  #### for loop over outcomes
  for(i in c(1:nrow(outcomes))) {
    # ith outcome name 
    outcome_name <- outcomes$outcome[i]
    # plotting matrix with pseudotime and ith outcome
    plot_mat <-data.frame(Pseudotime=as.numeric(as.matrix(hsmm@phenoData@data$Pseudotime)),
                          outcome=as.numeric(as.matrix(hsmm@phenoData@data[[outcome_name]])))
    # correlation between the pseudotime and ith outcome
    pval <- cor.test(plot_mat$Pseudotime, plot_mat$outcome, method='kendall')$p.value
    # if ith outcome is a numeric outcome
    if(outcomes$outcome_type[i]=='numeric'){
      # outcome values
      y <- plot_mat$outcome
      # median for splitting data into two
      med_y <- median(y, na.rm=T) # median
      y[which(plot_mat$outcome <= med_y)] <- 1 # less than median
      y[which(plot_mat$outcome >= med_y)] <- 2 # more than median
      plot_mat$factor_outcome <- as.factor(as.matrix(y)) # factorize it
      
    } else {
      plot_mat$factor_outcome <- plot_mat$outcome # create a dummy outcome variable
      plot_mat$factor_outcome <- as.factor(as.matrix(y)) # factorize it
    }
    # plot the trajectory annotated with ith outcome
    g <- plot_cell_trajectory(hsmm,color_by = plot_mat$factor_outcome,
                              show_branch_points=F,use_color_gradient = F,cell_size = 1.5) +
      ggplot2::labs(color=outcome_name) +
      ggplot2::ggtitle(sprintf("Kendall correlation pvalue: %0.3g", pval))
    # add plot to list and increment counter
    plot_list[[k]] <- g ; k<- k+1
    # plot the boxplot
    g <- ggplot2::ggplot(plot_mat, 
                         aes(x=factor_outcome, y=scale(Pseudotime,center=F),
                             fill=factor_outcome)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme(axis.text=element_text(size=15), 
                     axis.title=element_text(size=15,face="bold"),
                     legend.text=element_text(size=15))  +
      ggplot2::labs(y="Pseudotime",x=outcome_name) +  theme_bw() +
      theme(legend.position = 'none') 
    # add plot to list and increment counter
    plot_list[[k]] <- g ; k<- k+1
  }
  return(plot_list)
}

get_trajectories <- function(D # summarized experiment
                             ){
  ## monocle data object with all metabolites
  hsmm <- newCellDataSet(cellData =  assay(D), 
                              phenoData = new("AnnotatedDataFrame", 
                              data = colData(D) %>% as.data.frame), 
                              featureData = new("AnnotatedDataFrame", 
                                                data = rowData(D) %>% as.data.frame), 
                              expressionFamily = VGAM::tobit())
  # create 2D projection --> use DDRTree method
  hsmm<- reduceDimension(hsmm, max_components=2, norm_method = 'none',
                               pseudo_expr = 0,relative_expr = T,
                               reduction_method = 'DDRTree')
  # order "cells" (individuals) and estimate pseudotime
  hsmm <- orderCells(hsmm) 
  return(hsmm)
}