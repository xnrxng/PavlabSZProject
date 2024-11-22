# Author: Rui Xiang Yu
# Date: 2024 November 21
# Example usage: make all

all: load_data \
  summary_figures \
  countcells_genes \
  umaps \
  pseudobulk_dea \
  reports\sz-pavlab-report.html \
  reports\sz-pavlab-report.pdf


load_data: R/1-loaddata.R
  Rscript R/1-loaddata.R

summary_figures: R/2-summaryfigures.R
  Rscript R/2-summaryfigures.R

countcells_genes: R/3-countcellsandgenes.R
  Rscript R/3-countcellsandgenes.R

umaps: R/4-szagedistributionandumap.R
  Rscript R/4-szagedistributionandumap.R

pseudobulk_dea: R/6-pseudobulkanddea.R
  Rscript R/6-pseudobulkanddea.R

reports/sz-pavlab-report.html: results reports/sz-pavlab-report.qmd
	quarto render reports/sz-pavlab-report.qmd --to html

reports/sz-pavlab-report.pdf: results reports/sz-pavlab-report.qmd
	quarto render reports/sz-pavlab-report.qmd --to pdf

clean: load_data \
  summary_figures \
  countcells_genes \
  umaps \
  pseudobulk_dea \
  reports\sz-pavlab-report.html \
  reports\sz-pavlab-report.pdf
