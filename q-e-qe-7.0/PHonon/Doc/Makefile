VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)

LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf developer_man.pdf dfpt_tetra.pdf dfpt_self_energy.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  all
all:  defs pdf html
pdf: $(PDFS)
html: user_guide developer_man

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/ developer_man/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf input_xx.xsl
	- rm -rf ../../Doc/INPUT_PH.*

user_guide: user_guide.pdf
	- rm -rf user_guide/
	latex2html \
                -t "User's Guide for the PHonon package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                user_guide.tex
	cd user_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in user_guide/user_guide.html"
	@echo "***"
	@echo ""


developer_man: developer_man.pdf
	rm -rf developer_man/
	latex2html \
                -t "User's Guide for the PHonon package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                developer_man.tex
	cd developer_man; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** Developer manual created in developer_man/developer_man.html"
	@echo "***"
	@echo ""


defs: input_xx.xsl INPUT_PH.txt INPUT_PH.html INPUT_DYNMAT.txt INPUT_DYNMAT.html \
	INPUT_MATDYN.txt INPUT_MATDYN.html \
	INPUT_Q2R.txt INPUT_Q2R.html \
	INPUT_POSTAHC.txt INPUT_POSTAHC.html link_on_main_doc

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

INPUT_PH.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_PH.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_DYNMAT.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_DYNMAT.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_MATDYN.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_MATDYN.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_Q2R.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_Q2R.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_POSTAHC.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_POSTAHC.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

link_on_main_doc:
	-( cd ../../Doc ; \
	ln -fs ../PHonon/Doc/INPUT_PH.html . ; \
	ln -fs ../PHonon/Doc/INPUT_PH.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_PH.txt . ; \
	ln -fs ../PHonon/Doc/INPUT_DYNMAT.html . ; \
	ln -fs ../PHonon/Doc/INPUT_DYNMAT.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_DYNMAT.txt . ; \
	ln -fs ../PHonon/Doc/INPUT_MATDYN.html . ; \
	ln -fs ../PHonon/Doc/INPUT_MATDYN.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_MATDYN.txt . ; \
	ln -fs ../PHonon/Doc/INPUT_Q2R.html . ; \
	ln -fs ../PHonon/Doc/INPUT_Q2R.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_Q2R.txt . ; \
	ln -fs ../PHonon/Doc/INPUT_POSTAHC.html . ; \
	ln -fs ../PHonon/Doc/INPUT_POSTAHC.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_POSTAHC.txt . )
