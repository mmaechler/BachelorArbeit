BachelorThesis-Engl.pdf: BachelorThesis-Engl.tex chapter1.tex chapter2.tex \
	chapter3.tex chapter4.tex References.bib App_llnorMmix.tex App_plots.tex
	pdflatex BachelorThesis-Engl.tex

chapter1.tex: chapter1.Rnw
	R CMD Sweave $^
chapter2.tex: chapter2.Rnw
	R CMD Sweave $^
chapter3.tex: chapter3.Rnw
	R CMD Sweave $^
chapter4.tex: chapter4.Rnw
	R CMD Sweave $^

App_llnorMmix.tex: App_llnorMmix.Rnw
	R CMD Sweave $^
App_plots.tex: App_plots.Rnw
	R CMD Sweave $^

.PHONY: bib
bib:
	make -B
	bibtex BachelorThesis-Engl
	pdflatex BachelorThesis-Engl.tex
	pdflatex BachelorThesis-Engl.tex

.PHONY: clean
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot,
