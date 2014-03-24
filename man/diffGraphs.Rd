\name{diffGraphs}
\alias{diffGraphs}
\title{This function returns the different nodes and edges between graph1 and graph2.}
\usage{
diffGraphs(graph1, graph2, colorNodes = TRUE, colors = c("#B3E2CD",
  "#FDCDAC"))
}
\arguments{
  \item{graph1}{graphNEL}

  \item{graph2}{graphNEL}

  \item{colorNodes}{logical}

  \item{colors}{character vector of colors. If
  colorNodes==TRUE these colors are used for graph1 and
  graph2 respectivley.}
}
\value{
Return the diff between the graphs.
}
\description{
This function returns the different nodes and edges between
graph1 and graph2. Layout options of graph1 are kept.
Coloring currently not implemented.
}
\examples{
# load data
 data(biopaxexample)
 pwid1 = "pid_p_100002_wntpathway"
 pwid2 = "pid_p_100146_hespathway"
 mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
 mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
 plotRegulatoryGraph(diffGraphs(mygraph1,mygraph2))
}
\author{
Frank Kramer
}

