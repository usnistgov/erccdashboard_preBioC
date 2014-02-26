\name{geneExprTest}
\alias{geneExprTest}
\title{Prepare differential expression testing results for spike-in analysis}
\usage{
geneExprTest(expDat)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}
}
\description{
Prepare differential expression testing results for
spike-in analysis
}
\details{
This function wraps the QuasiSeq differential expression
testing package for datType = "count" or uses the limma
package for differential expression testing if datType =
"array". Alternatively, for count data only, if a correctly
formatted csv file is provided with the necessary DE test
results, then geneExprTest will bypass DE testing (with
reduced runtime). The function will look for a csv file
with the name "filenameRoot.quasiSeq.res.csv" and columns
corresponding to "Feature", "pvals", and "qvals" should be
in the file.
}

