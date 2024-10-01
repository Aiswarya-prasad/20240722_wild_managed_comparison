working_dir <- "/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20240722_wild_managed_comparison"
setwd(working_dir)

library(ggplot2)
library(readxl)
library(knitr)
library(RColorBrewer)
library(scales)
library(dplyr)
library(vegan)
library(ggsignif)
library(plotly)
library(htmlwidgets)
library(viridis)
library(hrbrthemes)
library(ggthemes)
library(ggrepel)
library(gridExtra)
library(emmeans)
library(microbiome)
library(ape)
library(phyloseq)
library(ComplexHeatmap)
library(ggnewscale)
library(VennDiagram)
library(ggforce)
library(circlize)
library(magrittr)
library(tidyverse)
library(corrplot)
library(ggkegg)
library(ggfx)
library(igraph)
library(tidygraph)
library(ggdendro)
library(dendextend)
library(ggpubr)
library(gggenes)

make_theme <- function(theme_name=theme_classic(), max_colors=0, palettefill="Pastel1", palettecolor="Dark2", modify_guide = T,
                        setFill=TRUE, setCol=TRUE,
                        guide_nrow=2, guide_nrow_byrow=TRUE, leg_pos="top", leg_size=12,
                        axis_x_title = 12, axis_y_title = 12,
                        x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                        y_angle=0 ,y_vj=0, y_hj=0, y_size=12){
  n_11 = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  n_12 = c("Paired", "Set3")
  n_8 = c("Accent", "Dark2", "Pastel2", "Set2")
  if (palettefill %in% n_12) {
    n_f = 12
  } else {
    if (palettefill %in% n_11) {
      n_f = 11
    } else {
      if (palettefill %in% n_8) {
        n_f  = 8
      } else {
        n_f = 9
      }
    }
  }
  if (palettecolor %in% n_12) {
    n_c = 12
  } else {
    if (palettecolor %in% n_11) {
      n_c = 11
    } else {
      if (palettecolor %in% n_8) {
        n_c  = 8
      } else {
        n_c = 9
      }
    }
  }
  getFill = colorRampPalette(brewer.pal(n_f, palettefill))
  getColor = colorRampPalette(brewer.pal(n_c, palettecolor))
  theme_params <- theme(axis.text.x = element_text(angle = x_angle,
    vjust = x_vj, hjust=x_hj,
    size = x_size),
    axis.text.y = element_text(angle = y_angle,
      vjust = y_vj, hjust=y_hj,
      size = y_size),
      axis.title.x = element_text(size=axis_x_title),
      axis.title.y = element_text(size=axis_y_title),
      legend.position=leg_pos,
      legend.text = element_text(size=leg_size)
    )
  if (modify_guide == T) {
    guide_params <- guides(fill = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  ),
                          col = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  )
                    )
  my_theme <- list(
                theme_name,
                theme_params,
                guide_params
              )
  } else {
    my_theme <- list(
                  theme_name,
                  theme_params
                )
  }
  if(setFill) {
    if (n_f < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_fill_manual(values = getFill(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_fill_brewer(palette=palettefill, na.value="grey")
                  )
    }
  }
  if(setCol) {
    if (n_c < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_color_manual(values = getColor(max_colors), na.value="grey")
                  )

    } else {
      my_theme <- list(
                    my_theme,
                    scale_color_brewer(palette=palettecolor, na.value="grey")
                  )
    }
  }
  return(my_theme)
}


sample_name_from_shortname <- function(shortname) {
  sample_name <- gsub("-", "_", shortname)
  return(paste0("DNA", sample_name))
}
