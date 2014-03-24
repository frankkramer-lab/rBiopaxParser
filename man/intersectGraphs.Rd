\name{intersectGraphs}
\alias{intersectGraphs}
\title{This function returns a graph computed by the insection of supplied graph1 and graph2.}
\usage{
intersectGraphs(graph1, graph2)
}
\arguments{
  \item{graph1}{graphNEL}

  \item{graph2}{graphNEL}
}
\value{
Returns the intersection of graph1 and graph2.
}
\description{
This function returns a graph computed by the insection of
supplied graph1 and graph2. Layout and weights of graph1
are kept.
}
\examples{
# load data
 data(biopaxexample)
 pwid1 = "pid_p_100002_wntpathway"
 pwid2 = "pid_p_100146_hespathway"
 mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
 mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
 plotRegulatoryGraph(intersectGraphs(mygraph1,mygraph2))
}
\author{
Frank Kramer
}

