# check the difference https://github.com/szcf-weiya/forked_seurat/compare/d9f09de15ddf05fe89a8b16eaa100e3720ee122b...v4.4.0-patch

#' Modified Seurat::DiffTTest by extending the output with statistics in addition to p-values
#'
#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using the Student's t-test
#'
#' @param data.use Data matrix to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param verbose Print a progress bar
#' @importFrom stats t.test
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#'
myDiffTTest <- function(
    data.use,
    cells.1,
    cells.2,
    verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_vals <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])[c("p.value", "statistic")]
      }
    )
  )
  p_val = p_vals[seq(1, length(p_vals), by = 2)]
  tstat = p_vals[seq(2, length(p_vals), by = 2)]

  to.return <- data.frame(p_val, tstat, row.names = rownames(x = data.use))
  return(to.return)
}

#' Modified Seurat::GLMDETest by extending the output with statistics in addition to p-values
#' @param data.use Data to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param min.cells Minimum number of cells threshold
#' @param latent.vars Latent variables to test
#' @param test.use parameterizes the glm
#' @param verbose Print progress bar
#
#' @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#' @importFrom MASS glm.nb
#' @importFrom stats var as.formula glm
myGLMDETest <- function(
    data.use,
    cells.1,
    cells.2,
    min.cells = 3,
    latent.vars = NULL,
    test.use = NULL,
    verbose = TRUE
) {
  group.info <- data.frame(
    group = rep(
      x = c('Group1', 'Group2'),
      times = c(length(x = cells.1), length(x = cells.2))
    )
  )
  rownames(group.info) <- c(cells.1, cells.2)
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars <- if (is.null(x = latent.vars)) {
    group.info
  } else {
    cbind(x = group.info, latent.vars)
  }
  latent.var.names <- colnames(x = latent.vars)
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_vals <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        latent.vars[, "GENE"] <- as.numeric(x = data.use[x, ])
        # check that gene is expressed in specified number of cells in one group
        if (sum(latent.vars$GENE[latent.vars$group == "Group1"] > 0) < min.cells &&
            sum(latent.vars$GENE[latent.vars$group == "Group2"] > 0) < min.cells) {
          warning(paste0(
            "Skipping gene --- ",
            x,
            ". Fewer than ",
            min.cells,
            " cells in both clusters."
          ))
          return(c(2, 0))
        }
        # check that variance between groups is not 0
        if (var(x = latent.vars$GENE) == 0) {
          warning(paste0(
            "Skipping gene -- ",
            x,
            ". No variance in expression between the two clusters."
          ))
          return(c(2, 0))
        }
        fmla <- as.formula(object = paste(
          "GENE ~",
          paste(latent.var.names, collapse = "+")
        ))
        p.estimate <- c(2, 0)
        if (test.use == "negbinom") {
          try(
            expr = p.estimate <- summary(
              object = glm.nb(formula = fmla, data = latent.vars)
            )$coef[2, c(4, 3)],
            silent = TRUE
          )
          return(p.estimate)
        } else if (test.use == "poisson") {
          return(summary(object = glm(
            formula = fmla,
            data = latent.vars,
            family = "poisson"
          ))$coef[2, c(4, 3)])
        }
      }
    )
  )
  p_val = p_vals[seq(1, length(p_vals), by = 2)]
  tstat = p_vals[seq(2, length(p_vals), by = 2)]
  features.keep <- rownames(data.use)
  if (length(x = which(x = p_val == 2)) > 0) {
    features.keep <- features.keep[-which(x = p_val == 2)]
    idx.keep = !p_val == 2
    p_val <- p_val[idx.keep]
    tstat = tstat[idx.keep]
  }
  to.return <- data.frame(p_val, tstat, row.names = features.keep)
  return(to.return)
}

#' Modified Seurat::WilcoxDETest by extending the output with statistics in addition to p-values
#'
#' @param data.use Data matrix to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param verbose Print a progress bar
#' @param ... Extra parameters passed to wilcox.test
#
#' @return Returns a p-value ranked matrix of putative differentially expressed features
#' @importFrom stats wilcox.test qnorm
myWilcoxDETest <- function(
    data.use,
    cells.1,
    cells.2,
    verbose = TRUE,
    ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  overflow.check <- ifelse(
    test = is.na(x = suppressWarnings(length(x = data.use[1, ]) * length(x = data.use[1, ]))),
    yes = FALSE,
    no = TRUE
  )
  limma.check <- TRUE # PackageCheck("limma", error = FALSE)
  if (limma.check[1] && overflow.check) {
    p_vals <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        res = rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])
        return(c(min(2 * min(res[1:2]), 1), mean(res[3:4]) ))
      }
    )
  } else {
    if (getOption('Seurat.limma.wilcox.msg', TRUE) && overflow.check) {
      message(
        "For a more efficient implementation of the Wilcoxon Rank Sum Test,",
        "\n(default method for FindMarkers) please install the limma package",
        "\n--------------------------------------------",
        "\ninstall.packages('BiocManager')",
        "\nBiocManager::install('limma')",
        "\n--------------------------------------------",
        "\nAfter installation of limma, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.limma.wilcox.msg = FALSE)
    }
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    p_vals <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        res = wilcox.test(data.use[x, ] ~ group.info[, "group"], conf.int = T, ...)
        # non-exact test, so go back to pnorm and assume two.sided test
        return(c(res$p.value, abs(qnorm(res$p.value/2)) * sign(res$estimate) ))
      }
    )
  }
  p_val = p_vals[seq(1, length(p_vals), by = 2)]
  tstat = p_vals[seq(2, length(p_vals), by = 2)]
  return(data.frame(p_val, tstat, row.names = rownames(x = data.use)))
}

#' Wilcoxon rank sum test (adapted from `limma::rankSumTestWithCorrelation`)
#' @param index any index vector such that `statistics[index]` contains the values of the statistic for the test group.
#' @param statistics numeric vector giving values of the test statistic.
#' @param correlation numeric scalar, average correlation between cases in the test group.  Cases in the second group are assumed independent of each other and other the first group.
#' @param df degrees of freedom which the correlation has been estimated.
#' @return Numeric vector of length 4 containing the `left.tail` and `right.tail` p-values and the corresponding test statistics.
#' @importFrom stats pt
rankSumTestWithCorrelation <- function(index, statistics, correlation = 0, df = Inf)
  #	Rank sum test as for two-sample Wilcoxon-Mann-Whitney test,
  #	but allowing for correlation between members of test set.
  #	Gordon Smyth and Di Wu
  #	Created 2007.  Last modified 24 Feb 2012.
{
  n <- length(statistics)
  r <- rank(statistics)
  r1 <- r[index]
  n1 <- length(r1)
  n2 <- n-n1
  U <- n1*n2 + n1*(n1+1)/2 - sum(r1)
  mu <- n1*n2/2

  if(correlation==0 || n1==1) {
    sigma2 <- n1*n2*(n+1)/12
  } else {
    sigma2 <- asin(1)*n1*n2 + asin(0.5)*n1*n2*(n2-1) + asin(correlation/2)*n1*(n1-1)*n2*(n2-1) + asin((correlation+1)/2)*n1*(n1-1)*n2
    sigma2 <- sigma2/2/pi
  }

  TIES <- (length(r) != length(unique(r)))
  if(TIES) {
    NTIES <- table(r)
    adjustment <- sum(NTIES*(NTIES+1)*(NTIES-1)) / (n*(n+1)*(n-1))
    sigma2 <- sigma2 * (1 - adjustment)
  }
  zlowertail <- (U+0.5-mu)/sqrt(sigma2)
  zuppertail <- (U-0.5-mu)/sqrt(sigma2)

  #	Lower and upper tails are reversed on output
  #	because R's ranks are the reverse of Mann-Whitney's ranks
  pvalues <- c(less=pt(zuppertail,df=df,lower.tail=FALSE), greater=pt(zlowertail,df=df), zuppertail, zlowertail)
  pvalues
}


#' Gene expression markers of identity classes
#'
#' Wrapper for `Seurat::FindMarkers` that uses modified (extended) tests
#'
#' Finds markers (differentially expressed genes) for identity classes
#'
#' @param object An object
#' @param ... Arguments passed to other methods and to specific DE methods
#' @export
myFindMarkers <- function(object, ...) {
  # environment(Seurat::FindMarkers)$DiffTTest <- myDiffTTest
  # environment(Seurat::FindMarkers)$GLMDETest <- myGLMDETest
  # environment(Seurat::FindMarkers)$WilcoxDETest <- myWilcoxDETest
  utils::assignInNamespace("DiffTTest", myDiffTTest, ns = "Seurat")
  utils::assignInNamespace("GLMDETest", myGLMDETest, ns = "Seurat")
  utils::assignInNamespace("WilcoxDETest", myWilcoxDETest, ns = "Seurat")

  result <- Seurat::FindMarkers(object, ...)
  return(result)
}
