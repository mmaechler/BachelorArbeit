BachelorThesis-Engl.pdf: BachelorThesis-Engl.tex chapter1.tex References.bib App_llnorMmix.tex
	pdflatex BachelorThesis-Engl.tex

chapter1.tex: chapter1.Rnw
	R CMD Sweave chapter1.Rnw

App_llnorMmix.tex: App_llnorMmix.Rnw
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
