% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Unico.r
\name{association_parametric}
\alias{association_parametric}
\title{Performs parametric statistical testing}
\usage{
association_parametric(
  X,
  Unico.mdl,
  slot_name = "parametric",
  diag_only = FALSE,
  intercept = TRUE,
  X_max_stds = 2,
  Q_max_stds = Inf,
  XQ_max_stds = Inf,
  parallel = TRUE,
  num_cores = NULL,
  log_file = "Unico.log",
  verbose = FALSE,
  debug = FALSE
)
}
\arguments{
\item{X}{An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations. Note that \code{X} must be the same \code{X} used to learn \code{Unico.mdl} (i.e. the original observed 2D mixture used to fit the model).}

\item{Unico.mdl}{The entire set of model parameters estimated by Unico on the 2D mixture matrix (i.e. the list returned by applying function \code{Unico} to \code{X}).}

\item{slot_name}{A string indicating the key for storing the results under \code{Unico.mdl}}

\item{diag_only}{A logical value indicating whether to only use the estimated source-level variances (and thus ignoring the estimate covariance) for controlling the heterogeneity in the observed mixture. if set to FALSE, Unico instead estimates the observation- and feature-specific variance in the mixture by leveraging the entire \code{k} by \code{k} variance-covariance matrix.}

\item{intercept}{A logical value indicating whether to fit the intercept term when performing the statistical testing.}

\item{X_max_stds}{A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the observed mixture value. Only samples whose observed mixture value fall within \code{X_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.}

\item{Q_max_stds}{A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the estimated mixture variance. Only samples whose estimated mixture variance fall within \code{Q_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.}

\item{XQ_max_stds}{A non-negative numeric value indicating, for each feature, the portions of data that are considered as outliers due to the weighted mixture value. Only samples whose weighted mixture value fall within \code{XQ_max_stds} standard deviations from the mean will be used for the statistical testing of a given feature.}

\item{parallel}{A logical value indicating whether to use parallel computing (possible when using a multi-core machine).}

\item{num_cores}{A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.}

\item{log_file}{A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.}

\item{verbose}{A logical value indicating whether to print logs.}

\item{debug}{A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.}
}
\value{
An updated \code{Unico.mdl} object with the the following list of effect size and p-value estimates stored in an additional key specified by \code{slot_name}
\item{gammas_hat}{An \code{m} by \code{k*p1} matrix of the estimated effects of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} covariates on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
\item{betas_hat}{An \code{m} by \code{p2} matrix of the estimated effects of the \code{p2} covariates in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
\item{gammas_hat_pvals}{An \code{m} by \code{k*p1} matrix of p-values for the estimates in \code{gammas_hat} (based on a T-test).}
\item{betas_hat_pvals}{An \code{m} by \code{p2} matrix of p-values for the estimates in \code{betas_hat} (based on a T-test).}
\item{gammas_hat_pvals.joint}{An \code{m} by \code{p1} matrix of p-values for the joint effects (i.e. across all \code{k} sources) of each of the \code{p1} covariates in \code{C1} on each of the \code{m} features in \code{X} (based on a partial F-test). In other words, these are p-values for the combined statistical effects (across all sources) of each one of the \code{p1} covariates on each of the \code{m} features under the Unico model.}
\item{Q}{An \code{m} by \code{n} matrix of weights used for controlling the heterogeneity of each observation at each feature (activated only if \code{debug == TRUE}).}
\item{masks}{An \code{m} by \code{n} matrix of logical values indicating whether observation participated in statistical testing at each feature  (activated only if \code{debug == TRUE}).}
\item{phi_hat}{An \code{m} by \code{k+p1*k+p2} matrix containing the entire estimated effect sizes (including those on source weights) for each feature (activated only if \code{debug == TRUE}).}
\item{phi_se}{An \code{m} by \code{k+p1*k+p2} matrix containing the estimated standard errors associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
\item{phi_hat_pvals}{An \code{m} by \code{k+p1*k+p2} matrix containing the p-values associated with \code{phi_hat} for each feature (activated only if \code{debug == TRUE}).}
}
\description{
Performs parametric statistical testing (T-test) on (1) the marginal effect of each covariate in \code{C1} at source-specific level (2) the joint effect across all sources for each covariate in \code{C1} (3) non-source-specific effect for each covariate in \code{C2}. In the context of bulk genomic data containing a mixture of cell types, these correspond to the marginal effect of each covariate in \code{C1} (potentially including the phenotype of interest) at each cell type, joint tissue-level effect for each covariate in \code{C1}, and tissue-level effect for each covariate in \code{C2}.
}
\details{
If we assume that source-specific values \eqn{Z_{ijh}} are normally distributed, under the Unico model, we have the following:
\deqn{Z_{ij} \sim \mathcal{N}\left(\mu_{j} + (c_i^{(1)})^T \gamma_{jh}, \sigma_{jh}^2 \right)}
\deqn{X_{ij} \sim \mathcal{N}\left(w_{i}^T (\mu_{j} + (c_i^{(1)})^T \gamma_{jh}) + (c_i^{(2)})^T \beta_j,  \text{Sum}\left((w_i w_i^T ) \odot \Sigma_j\right) + \tau_j^2\right)}
For a given feature \eqn{j} under test, the above equation corresponds to a heteroskedastic regression problem with \eqn{X_{ij}} as the dependent variable and \eqn{\{\{w_i\}, \{w_i c_i^{(1)}\}, \{c_i^{(2)}\}\}} as the set of independent variables.
This view allows us to perform parametric statistical testing (T-test for marginal effects and partial F-test for joint effects) by solving a generalized least squares problem with sample \eqn{i} scaled by the inverse of its estimated standard deviation.
}
\examples{
data = simulate_data(n=100, m=2, k=3, p1=1, p2=1, taus_std=0, log_file=NULL)
res = list()
res$params.hat = Unico(data$X, data$W, data$C1, data$C2, parallel=FALSE, log_file=NULL)
res$params.hat = association_parametric(data$X, res$params.hat, parallel=FALSE, log_file=NULL)

}
