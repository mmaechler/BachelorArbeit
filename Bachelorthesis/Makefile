BachelorThesis-Engl.pdf : BachelorThesis-Engl.tex chapter1.tex References.bib
	pdflatex BachelorThesis-Engl.tex

chapter1.tex : chapter1.Rnw
	R CMD Sweave chapter1.Rnw

.PHONY: bib
bib:
	pdflatex BachelorThesis-Engl.tex
	bibtex BachelorThesis-Engl
	pdflatex BachelorThesis-Engl.tex
	pdflatex BachelorThesis-Engl.tex

.PHONY: clean
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot,
