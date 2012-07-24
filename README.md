# rBiopaxParser

Parses BioPax level 2 files and represents them in R.  
https://github.com/frankkramer/rBiopaxParser  
  
More concretely, `rBiopaxParser`:

 * Can download Biopax data from online resources

 * Parses Biopax owl files

 * Features many functions to select, add, modify and remove from the parsed Biopax models

 * Can produce regulatory graphs from Biopax pathway data and offers functions to merge, diff and transform these graphs in various ways 
 
 * Visualization functions to layout ina a (more or less) beautiful way
 
 * Can write out the (modified) parsed Biopax models 
  
*Prerequisites:*  
This package uses package RCurl to download Biopax files from the web.  
This package uses package XML to parse the Biopax .owl files.  
This package uses package graph and package Rgraphviz to visualize networks.  
To install directly from github you need package devtools installed.
  
Installation or running certain functions MIGHT fail if these are not met. Please read through the following instructions.   
  
*Installing prerequisites for Linux users:*  
XML:   
Make sure your linux has library libxml2 installed. This is almost always the case. Otherwise run in your shell  
<code>
	sudo apt-get install libxml2
</code>  
will fix this issue. You will now be able to install R package XML, this should be automatically done when you install rBiopaxParser, or you can run within R:  
<code>
	install.packages("XML")
</code>  

RCurl:   
RCurl is only needed for a convenience function to download Biopax files directly within R. You can skip this step if you already have the Biopax data downloaded.  
Make sure your linux has library libcurl installed and curl-config in your path. Check out 
<code>
	locate libcurl
	locate curl-config
</code>  
If these are not found (usually the developer version is missing), most Linux users can usually fix this by running   
<code>
	sudo apt-get install libcurl4-openssl-dev
</code>    
You will now be able to install R package RCurl, this should be automatically done when you install rBiopaxParser, or you can run within R:  
<code>
	install.packages("RCurl")
</code>  
If you encounter more problems check out http://www.omegahat.org/RCurl/FAQ.html  
   
  
graph:  
Package graph has moved from CRAN to Bioconductor recently, you might encounter an error saying that package graph is not available for your distribution when calling install.packages("graph").  
Check out http://bioconductor.org/packages/release/bioc/html/graph.html or call
<code>
	source("http://bioconductor.org/biocLite.R")
    biocLite("graph")
</code>  
to install it.  
  
  
Rgraphviz:   
Rgraphviz is used to layout the graphs generated in this package. You can layout and plot these yourself if you want to.  
Make sure your linux has package graphviz installed.  
If this is not the case, many linux users can usually fix this by running   
<code>
	sudo apt-get install graphviz
</code>    
You will now be able to install R package Rgraphviz using:  
<code>
	source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
</code>  
If you encounter more problems check out http://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html  
  
  
devtools:  
Package devtools is available at CRAN. Call
<code>
	install.packages("devtools")
</code>  
to install it.  
  
  
*Installing prerequisites for Windows users:*  
XML & RCurl:   
These packages depend on linux libraries, however Brian Ripley has put some work into this to enable Windows users.  
Check out  http://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/ for these two packages for your R version.  
Download first XML...zip and then RCurl...zip and install them locally on your machine.  
  
  
graph:  
Package graph has moved from CRAN to Bioconductor recently, you might encounter an error saying that package graph is not available for your distribution when calling install.packages("graph").  
Check out http://bioconductor.org/packages/release/bioc/html/graph.html or call
<code>
	source("http://bioconductor.org/biocLite.R")
    biocLite("graph")
</code>  
to install it.  
  
  
Rgraphviz:   
Rgraphviz is used to layout the graphs generated in this package. You can layout and plot these yourself if you want to.  
Make sure your  machine has Graphviz installed, it can be found at: http://www.graphviz.org  
Click on Download -> Windows. 
After installing graphviz you will now be able to install R package Rgraphviz using:  
<code>
	source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
</code>  
If you encounter more problems check out http://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html  
  
  
devtools:  
Package devtools is available at CRAN. For Windows this seems to depend on having Rtools for Windows installed. You can download and install this from:  
http://cran.r-project.org/bin/windows/Rtools/  
To install R package devtools call  
<code>
	install.packages("devtools")
</code>  
  
  
Finally:    
*Installing rBiopaxParser from github:*  
<code>
	install.packages("devtools")
	library(devtools)
	install_github(repo="rBiopaxParser", username="frankkramer")
</code>


 
 