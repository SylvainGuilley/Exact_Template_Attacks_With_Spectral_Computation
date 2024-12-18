#!/bin/sh

#-------------------------------------------------------------------------------
# Filename:    fig2pdf.sh
# Description: Bourne shell script to compile an XFIG file with LaTeX commands
#              into a PDF document
# Author:      Sylvain GUILLEY (sylvain.guilley@TELECOM-ParisTech.fr)
# Remarks:      - Cross-references with the main document are allowed.
#               - It is usefull to have the following alias for xfig:
#                 xfig -specialtext -latexfonts -startlatexFont default
#-------------------------------------------------------------------------------

# Usage:
if [ -z $1 ]; then
  echo "Usage: $0 <file>"
  echo "$0 produces <file>.pdf from <file>.fig, where <file>.fig contains arbitrary drawings and LaTeX mathematical formulas"
  exit 1;
fi;

doc=$1

# Creation of multi-threaded proof working copies:
# Also prevents from clobbering an unlucky file "driver.tex"
sed "s/\<temp\>/$doc.$$/g" .driver.tex >driver.$$.tex
cp $doc.fig $doc.$$.fig
# Creation of the postscript layer:
fig2dev -L pstex $doc.$$.fig > $doc.$$.pstex
# Creation of the LaTeX layer (adjusted over the PS):
fig2dev -L pstex_t -p $doc.$$.pstex $doc.$$.fig > $doc.$$.pstex_t
# Generating a blank (apart from the graphic) DVI document:
latex driver.$$.tex >/dev/null
# Copying the "\label" defined in ../ in the driver auxiliary file:
# grep -h "newlabel" ../*.aux >> driver.$$.aux
# grep -h "citation" ../*.aux >> driver.$$.aux
# grep -h "bibcite"  ../*.aux >> driver.$$.aux
# Converting the DVI to EPS
dvips -Ppdf -G0 -E driver.$$.dvi -o $doc.eps 2>&1 > /dev/null
epstopdf $doc.eps; # With gs, we would have needed "-dEPSCrop"
# Removing temporay files
rm -f \
	$doc.eps \
  $doc.$$.fig $doc.$$.pstex $doc.$$.pstex_t \
  driver.$$.tex driver.$$.dvi driver.$$.aux driver.$$.log
