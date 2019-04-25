.PHONY: all clean submission

all: spaceness.pdf

# get these by doing
#   echo $(grep includegr spaceness.tex | sed -e 's/.*{\([^}]*\)}.*/\1/')
FIGURES = figures/sampling_maps.pdf figures/pop_params.pdf figures/sfs_w_sumstats.pdf figures/stairwayplot_facet_rollmean.pdf figures/gwas_summary.pdf figures/sumstats_by_neighbors_allstats.pdf figures/fig_S1_sfs_grid_model_by_sampling.pdf figures/het_z_by_ind.pdf 

spaceness.pdf : spaceness.bib $(FIGURES)

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	-bibtex $<

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

%.fonts.pdf : %.pdf
	# this will embed all fonts
	inkscape $< --export-pdf=$@

