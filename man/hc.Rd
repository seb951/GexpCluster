% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render_plots.R
\name{hc}
\alias{hc}
\title{Hierarchical clustering plot}
\usage{
hc(gexp = data[[2]], max_genes = 200)
}
\arguments{
\item{gexp}{gexp df}

\item{max_genes}{max nb of genes to show}
}
\value{
the heatmap
}
\description{
Convert a dataframe of Gexpr to a heatmap
}
\examples{
heatmap = hc();
}
