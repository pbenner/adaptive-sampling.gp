TEX = presentation.tex
BIB = ../../literature.bib
TARGET = $(TEX:.tex=.pdf)
#LOGOS = minerva.pdf mislogo.pdf
FIGURES = correlations1.pdf correlations2.pdf correlations3.pdf correlations4.pdf  

all: $(TARGET)

$(TARGET): $(TEX) $(FIGURES) allreferences.bbl #$(LOGOS)
	lualatex $<
	lualatex $<


# list of all references, to be used to extract single bibitems in place:
allreferences.bbl: $(BIB)
	printf "\\\\bibdata{$(<:.bib=)}\n\\\\citation{*}\n\\\\bibstyle{plain}" > allreferences.aux
	bibtex allreferences
	rm allreferences.aux allreferences.blg


#### Figures:

correlations%.pdf: correlations.r
	R --vanilla < $^

%.pdf: %.r
	R --vanilla < $^

%.pdf : %.svg
	inkscape -A $@ $^

%.pdf: %.eps
	epstopdf $^

clean: 
	rm -f $(TARGET) $(TEX:.tex=.aux) $(TEX:.tex=.log) $(TEX:.tex=.nav) $(TEX:.tex=.snm) $(TEX:.tex=.out) $(TEX:.tex=.toc) 
	
update:
	svn update
	touch *.eps
	touch *.pdf
