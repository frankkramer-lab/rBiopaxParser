\name{diffGraphs}
\alias{diffGraphs}
\title{This function returns the different nodes and edges between graph1 and graph2.}
\usage{
  diffGraphs(graph1, graph2, colorNodes = TRUE,
    colors = c("#B3E2CD", "#FDCDAC"))
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
  This function returns the different nodes and edges
  between graph1 and graph2. Layout options of graph1 are
  kept. Coloring currently not implemented.
}
\author{
  Frank Kramer
}

