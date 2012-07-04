# rBiopaxParser

Parses BioPax level 2 files and represents them in R.
https://github.com/frankkramer/rBiopaxParser

Installation from github:
<code>
	install.packages("devtools")
	library(devtools)
	install_github(repo="rBiopaxParser", username="frankkramer")
</code>

More concretely, `rBiopaxParser`:

 * Can download Biopax data from online resources

 * Parses Biopax owl files

 * Features many functions to select, add, modify and remove from the parsed Biopax models

 * Can produce regulatory graphs from Biopax pathway data and offers functions to merge, diff and transform these graphs in various ways 
 
 * Visualization functions to layout ina a (more or less) beautiful way
 
 * Can write out the (modified) parsed Biopax models 
 
 