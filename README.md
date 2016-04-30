# Visual Sample Plan and Prior Information: What do we Need to Know to Find UXO?
## Kenneth A Flagg, Montana State University

Welcome to the online home of my MS writing project, featuring code and data files. This repository contains everything you need to reproduce my analyses and recreate my sldies and paper without actually running the simulations. It also includes all the code so you can run the simulations yourself. I have run everything successfully on both Windows 10 and Debian.

This project is a simulation study investigating the influence of prior information inputs on the results produced by the unexploded ordnance target area delineation features of [Visual Sample Plan](http://vsp.pnnl.gov). This is an initial exploratory analysis using fictitious sites. I am making all of my code publically available to support discussions about my conclusions and ways to improve similar simulation studies in the future. You are welcome to send any questions, comments, or concerns to kenneth.flagg@msu.montana.edu.

You can find the full paper [here](writeup/flagg_writeup.pdf) and an outline of the 15-minute presentation [here](presentation/flagg_presentation.md).

# Software Needed to Build the Presentation and Paper

* [TeX Live](https://www.tug.org/texlive/) or another LaTeX system
* [Pandoc](http://www.pandoc.org/)
* [R](http://www.r-project.org)
* R packages:
    * [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html)
    * [knitr](http://www.yihui.name/knitr/)
    * [extrafont](https://cran.r-project.org/web/packages/extrafont/index.html)
    * [rmarkdown](http://rmarkdown.rstudio.com/)
    * [revealjs](https://cran.r-project.org/web/packages/revealjs/index.html)

If you use [RStudio](https://www.rstudio.com/), you probably have all of the above already set up except for spatstat, extrafont, and revealjs.

I use extrafont to make R use a clone of the Computer Modern font for the plots in the paper. After installing the extrafont package, you need to run the R command
```r
font_install('fontcm')
```
to install the appropriate font package.

# Building the Presentation

To build the presentation slides, open an R terminal in the `presentation` directory and run
```r
rmarkdown::render('flagg_presentation.rmd')
```
to create the standalone HTML presentation.

Or, if you are of the RStudio persuation, simply open `flagg_presentation.rmd` in RStudio and click the "knit" button.

# Building the Paper

To build the paper, open a terminal or command prompt (__not__ an R terminal) in the `writeup` directory and run
```bash
Rscript -e "knitr::knit('flagg_presentation.rnw')"
```
to create the LaTeX file. Next, run
```bash
pdflatex flagg_presentation.tex
```
to create the PDF. You will probably get messages about undefined references. If so, run
```bash
pdflatex flagg_presentation.tex
```
again to cross-reference the figures and bibliography. Finally, run
```bash
Rscript -e "extrafont::embed_fonts('flagg_presentation.pdf')"
```
to embed the correct font.

If you prefer RStudio, open `flagg_writeup.rnw` in RStudio and click the "Compile PDF" button. You still need to embed the fonts, so run
```r
extrafont::embed_fonts('flagg_presentation.pdf')
```
in RStudio's terminal.

# Software Needed to Run the Simulations

All of these are free, but you need to create an account at the VSP site before downloading it.

* [Visual Sample Plan](http://vsp.pnnl.gov)
* [GSLIB (64-bit)](http://www.gslib.com)
* [R](http://www.r-project.org)
* R packages:
    * [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html)
    * [RDCOMClient](http://www.omegahat.net/RDCOMClient/)

# Running the Simulations

The simulations are meant to be run from the root directory of the repository. Open an R terminal and run
```r
source('datasets/easy.r')
source('datasets/medium.r')
source('datasets/hard.r')
```
to generate 3,000 realizations of each site. Each of these commands could take the better part of a day to run, so be patient.

Next, you will want to create some results to analyze. Make sure the GSLIB executables are somewhere in your path, then run
```r
source('datasets/experiment1.r')
source('datasets/experiment2e.r')
source('datasets/experiment2m.r')
source('datasets/experiment2h.r')
```
and go do something else for a week or so while your computer computes many intersections between circles and rectangles.

The only Windows-only part of this project is the illustration of the sampling variability in VSP's detection probability simulation. You can run it with
```r
source('datasets/spacings.r')
```
which should only take a few minutes. Yay.
