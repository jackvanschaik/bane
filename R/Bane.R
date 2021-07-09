#' Bane class
#' @title Bane Class
#' @docType class
#' @description R6 class to implement the BaNE model
#' @field indep independent variable names
#' @field dep list of variables dependencies
#' @field mu independent variable priors
#' @field lm dependent variable priors
#' @field n Number of independent nodes
#' @field m Number of dependent nodes
#' @field k Number of total nodes
#' @field N Rows of data
#' @field data Data frame containing phenotype data
#' @field all_conds Character vector of all condtion names
#' @field param_names Names of all model parameters
#' @field J Total number of parameters
#' @field M Structure of the parameter matrix
#' @field M_ A helper object
#' @field q_ A helper object
#' @field mu_pr The mean vector for the multivariate normal prior
#' @field Sig_inv The precision matrix for the multivariate normal prior
#' @field LD_Data Data formatted as a list for LaplacesDemon
#' @field ld Output from LaplacesDemon, the posterior sampler
#'
#' @importFrom R6 R6Class
#' @export
Bane <- R6::R6Class("Bane",
                public = list(
                    indep = NULL,
                    dep = NULL,
                    mu = NULL,
                    lm = NULL,
                    n = NULL,
                    m = NULL,
                    k = NULL,
                    N = NULL,
                    data = NULL,
                    all_conds = NULL,
                    param_names = NULL,
                    J = NULL,
                    M = NULL,
                    M_ = NULL,
                    q_ = NULL,
                    mu_pr = NULL,
                    Sig_inv = NULL,
                    LD_Data = NULL,
                    ld = NULL,
                    #' @description Intialize the R6 object
                    #' @return A new `Bane object`
                    intialize = function() {

                    },
                    #' @description Create a new BaNE model
                    #' @param indep A character vector of independent variable names.
                    #' @param dep A named list of dependent variables with their dependencies. List names are the dependent variable names and list elements are the corresponding vector of dependency names.
                    #' @param mu Prior proportions for independent variables.
                    #' @param lm Prior baseline proportions for dependent variables.
                    #' @param data A dataframe whose columns are referenced by dep and indep
                    #' @return Updates the object with side effects
                    create_model = function(indep, dep, mu, lm, data) {
                        # Initialize key model parameters
                        self$indep <- indep
                        self$dep <- dep
                        self$mu <- mu
                        self$lm <- lm
                        self$n <- length(indep)
                        self$m <- length(dep)
                        self$k <- self$m + self$n
                        self$N <- nrow(data)
                        self$all_conds <- c(indep, names(dep))

                        # Validate the model configuration
                        if (!all(unlist(lapply(dep, function(x) all(x %in% self$all_conds))))) {
                            stop("Missing head conditions")
                        }

                        if (!all(self$all_conds %in% names(data))) {
                            stop("Missing conditions in data")
                        }

                        # Add intercept to data
                        self$data <- as.matrix(cbind(
                            data.frame(intercept = rep(1, nrow(data))),
                            data[self$all_conds]
                        ))

                        # Prepare design matrix
                        alpha <- c()
                        beta <- c()
                        A <- matrix(0, nrow = self$m, ncol = self$n)
                        B <- matrix(0, nrow = self$m, ncol = self$m)
                        M <- matrix(0, nrow = self$k, ncol = self$k + 1)
                        colnames(M) <- c("intercept", self$all_conds)
                        rownames(M) <- self$all_conds

                        # Create parameter names / match with design matrix
                        for (i in 1:self$m) {
                            cond <- dep[[i]]
                            for (j in 1:self$n) {
                                if (indep[j] %in% cond) {
                                    alpha <- c(alpha, paste0("a", i, j))
                                    A[i,j] <- 1
                                }
                            }
                            for (j in 1:self$m) {
                                if (names(dep)[j] %in% cond) {
                                    beta <- c(beta, paste0("b", i, j))
                                    B[i,j] <- 1
                                }
                            }
                        }

                        # Finish design matrix
                        self$param_names <- c(paste0("mu", 1:self$n), paste0("lm", 1:self$m), alpha, beta)
                        self$J <- length(self$param_names)
                        M[,1] <- 1
                        M[(self$n+1):self$k,2:(self$n+1)] <- A
                        M[(self$n+1):self$k,(self$n+2):(self$k+1)] <- B
                        self$M <- M

                        # These are some helper objects, to speed up likelihood calculations
                        self$M_ <- as.numeric(M)
                        self$q_ <- which(self$M_ != 0)

                        # Hyperparameters for multivariate normal prior
                        self$mu_pr <- rep(0, self$J)
                        self$mu_pr[1:self$k] <- c(invsig(mu), invsig(lm))
                        self$Sig_inv <- solve(diag(length(self$mu_pr)))

                        # Data object per LaplacesDemon format
                        self$LD_Data <- list(
                            J = self$J,
                            X = self$data,
                            N = nrow(self$data),
                            mon.names = "LP",
                            parm.names = self$param_names
                        )
                    },

                    #' @description Helper function to structure parameters as a matrix
                    #' @param p Vector of model parameters
                    #' @return A matrix of size (k+1) x k
                    params_to_mat = function(p) {
                        Q_ <- rep(0, length(self$M_))
                        Q_[self$q_] <- p
                        matrix(Q_, nrow = self$k)
                    },
                    #' @description Log likelihood
                    #' @param p Vector of model parameters
                    #' @param data A data matrix
                    #' @return The log likelihood value of `p` given `data`
                    LLX = function(p, data) {
                        M_p <- self$params_to_mat(p)
                        s <- data %*% t(M_p)
                        s[] <- vapply(s, sig, numeric(1))
                        s_ <- 1 - s
                        s[] <- vapply(s, log, numeric(1))
                        s_[] <- vapply(s_, log, numeric(1))
                        sum((data[,-1] * s) + (1 - data[,-1])*s_)
                    },
                    #' @description Log prior
                    #' @param p Vector of model parameters
                    #' @return The log prior evaluated at `p`
                    LPR = function(p) {
                        -0.5 * mahalanobis(p, self$mu_pr, self$Sig_inv, inverted = TRUE)
                    },
                    #' @description Log posterior
                    #' @param p Vector of model parameters
                    #' @param data A data matrix
                    #' @return A posterior evalutions in LaplacesDemon format
                    LD_Model = function(p, data) {
                        ll <- self$LLX(p, data$X)
                        lpr <- self$LPR(p)
                        lpo <- ll + lpr

                        list(
                            LP = lpo,
                            Dev = -2 * ll,
                            Monitor = lpo,
                            yhat = 1,
                            parm = p
                        )
                    },
                    #' @description Run the HMC sampler
                    #' @param Iterations passed to `LaplacesDemon`
                    #' @param Status passed to `LaplacesDemon`
                    #' @param Thinning passed to `LaplacesDemon`
                    #' @param eps passed to `LaplacesDemon`
                    #' @param L passed to `LaplacesDemon`
                    #' @return A LaplacesDemon object contain sampler info and posterior draws
                    run_chain = function(Iterations = 2000, Status = 100, Thinning = 5,
                                         eps = 0.1, L = 3) {
                        self$ld <- LaplacesDemon::LaplacesDemon(
                            Model = self$LD_Model,
                            Data = self$LD_Data,
                            Initial.Values = rep(0, self$J),
                            Covar = NULL,
                            Iterations = Iterations,
                            Status = Status,
                            Thinning = Thinning,
                            Algorithm = "HMC",
                            Specs = list(epsilon = rep(eps, self$J), L = L, m = 0.5* diag(self$J))
                        )
                    },
                    #' @description  Mode/ maximum likelihood
                    #' @return Output from the `optim` function maximizng the likelihood function
                    maximum_likelihood = function() {
                        optim(
                            rep(0, self$J),
                            bane_2$LLX,
                            data = bane_2$data,
                            control = list(fnscale = -1, maxit = 1000)
                        )
                    },
                    #' @description Posterior sampled subcohort proportions
                    #' @return A list with subcohort, corresponding proportion draws, and a plot
                    post_subs = function() {
                        v <- length(self$all_conds)
                        comb_mat <- matrix(unlist(lapply(0:((2^v)-1), bitwShiftR, (v-1):0)) %% 2, nrow = v)
                        comb_cnd <- t(rbind(rep(1, 2^v), comb_mat))

                        colnames(comb_cnd) <- c("intercept", self$all_conds)

                        post_subs <- apply(self$ld$Posterior1, 1, function(y) {
                            M <- self$params_to_mat(y)
                            s <- comb_cnd %*% t(M)
                            s[] <- vapply(s, sig, numeric(1))
                            s_ <- 1 - s
                            u <- (s ^ comb_cnd[,-1]) * (s_^ (1 - comb_cnd[,-1]))
                            apply(u, 1, prod)
                        })

                        sub_names <- colnames(comb_cnd)[-1]
                        Z <- do.call(rbind, lapply(1:nrow(comb_cnd), function(i) {
                            sc_data <- tibble::tibble(
                                `Subcohort Name` = paste(paste(sub_names, comb_cnd[i,2:ncol(comb_cnd)], sep="_"), collapse = "_"),
                                `Subcohort Proportion` = post_subs[i,]
                            )
                        }))

                        gg <- ggplot2::ggplot(Z, ggplot2::aes(x = `Subcohort Proportion`, group = `Subcohort Name`, fill = `Subcohort Name`)) +
                            ggplot2::geom_density(alpha = 0.5) +
                            ggplot2::labs(y = "Density", title = "Posterior distributions of subcohort proportions") +
                            ggplot2::theme_minimal()

                        list(subcohorts = comb_cnd, post_subs = post_subs, ggplot = gg)
                    },
                    #' @description Plot the phenotype topology
                    #' @return A ggplot with phenotype topology
                    plot_topology = function() {
                        dep <- self$dep
                        indep <- self$indep
                        L <- lapply(seq(dep), function(i) data.frame(nout = dep[[i]], nin = names(dep)[i]))
                        ig <- igraph::graph_from_edgelist(as.matrix(do.call(rbind, L)))
                        ggn <- ggnetwork::ggnetwork(ig, layout = igraph::layout_as_tree(ig, root = indep, mode = "out"))
                        ggn$dep <- ifelse(ggn$name %in% indep, "independent", "dependent")

                        gg <- ggplot2::ggplot(ggn, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, label = name)) +
                            ggnetwork::geom_edges(color = "black", arrow = ggplot2::arrow(length = ggnetwork::unit(6, "pt"), type = "closed")) +
                            ggnetwork::geom_nodes(ggplot2::aes(color = dep), size = 6) +
                            ggnetwork::geom_nodelabel_repel(ggplot2::aes(color = dep), box.padding = ggnetwork::unit(1, "lines"), alpha = 0.75) +
                            ggnetwork::theme_blank()

                        gg
                    },
                    #' @description Print details about the object
                    #' @return Side effects -- printing
                    print = function() {
                        print("Independent conditions: ")
                        print(self$indep)
                        print("Dependent conditions:" )
                        print(self$dep)
                        print("Prior means for independent conditions: ")
                        print(self$mu)
                        print("Prior means for dependent conditions: ")
                        print(self$lm)
                        print("Total number of conditions: ")
                        print(self$k)
                        print("Sample size: ")
                        print(self$N)
                        print("Data preview: ")
                        print(head(self$data))
                        print("Parameter names: ")
                        print(self$param_names)
                        print("Number of parameters: ")
                        print(self$J)
                        print("Design Matrix: ")
                        print(self$M)
                        print("Prior mean: ")
                        print(self$mu_pr)
                        print("Prior precision matrix: ")
                        print(self$Sig_inv)
                    }
                )
)
