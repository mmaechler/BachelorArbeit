BachelorThesis-Engl.pdf : BachelorThesis-Engl.tex chapter1.tex chapter1.Rnw
	pdflatex BachelorThesis-Engl.tex

BachelorThesis-Engl.tex : chapter1.tex chapter1.Rnw

chapter1.tex : chapter1.Rnw
	R CMD Sweave chapter1.Rnw
