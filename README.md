# Visual Sample Plan and Prior Information: What do we Need to Know to Find UXO?
## Kenneth A Flagg, Montana State University

This will soon become the online home of my MS writing project, featuring code and data files. The project is currently hosted on MSU's GitLab. Please email kenneth.flagg@msu.montana.edu to request access.

You can find the presentation [here](presentation/flagg_presentation.md). The paper will be posted once it has been approved.

# Software Needed to Run the Simulations

All of these are free, but you need to create an account at the VSP site before downloading it.

* [Visual Sample Plan](http://vsp.pnnl.gov)
* [GSLIB (64-bit)](http://www.gslib.com)
* [R](http://www.r-project.org)
* R packages:
    * [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html)
    * [RDCOMClient](http://www.omegahat.net/RDCOMClient/)

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
