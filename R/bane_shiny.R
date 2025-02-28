#' BaNE Shiny App
#'
#' Runs a GUI for easily implementing the BaNE model
#'
#' @param data A data frame to use the model with
#'
#' @return nothing, usually
#' @export
#' @import shiny
#' @import shinydashboard
#' @import reactable
bane_shiny <- function(data) {
    bane_ <- Bane$new()

    # UI #######################################################################
    ui <- dashboardPage(
        dashboardHeader(title = "BaNE Model Implementation"),
        dashboardSidebar(
            sidebarMenu(
                menuItem("Model Specification", tabName = "model", icon = icon("code-branch")),
                menuItem("Prior Specification", tabName = "prior", icon = icon("sliders-h")),
                menuItem("Model Details", tabName = "details", icon = icon("th")),
                menuItem("Run Posterior Sampler", tabName = "posterior", icon = icon("calculator")),
                menuItem("Posterior Intervals", tabName = "intervals", icon = icon("poll"))
            )
        ),
        dashboardBody(
            tabItems(
                tabItem(tabName = "model",
                        fluidRow(
                            box(width = 12,
                                h3("Condition/Outcome Dependencies"),
                                h4("Allowed Column Names:"),
                                textOutput("col_txt"),
                                textInput("idp_txt", "Independent Conditions", placeholder = "cond_1,cond_2,cond_3"),
                                textAreaInput("dep_txt", "Dependent Conditions", placeholder = "cond_4:cond2\ncond_5:cond_1,cond_4"),
                                actionButton("bane_create", "Create Model")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                h3("Phenotype Topology:"),
                                plotOutput("top_plt")
                            )
                        )

                ),
                tabItem(tabName = "prior",
                        fluidRow(
                            box(width = 12,
                                h3("Prior Parameters"),
                                actionButton("pri_reset", "Reset Page")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                h3("Prior Input"),
                                h4("Prior Mean"),
                                uiOutput("pri_mean"),
                                h4("Prior Covariance Matrix"),
                                uiOutput("pri_cov"),
                                actionButton("pri_use", "Use Parameters")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                h3("Beta to Normal Parameter Estimator"),
                                textInput("beta_a", "Alpha"),
                                textInput("beta_b", "Beta"),
                                actionButton("beta_est", "Estimate"),
                                htmlOutput("est_details"),
                                plotOutput("plt_beta"),
                                plotOutput("plt_norm")
                            )
                        )

                ),
                tabItem(tabName = "details",
                        fluidRow(
                            box(width = 12,
                                h3("Model Details"),
                                actionButton("mod_dtl", "Show Model Details")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                verbatimTextOutput("mod_out")
                            )
                        )

                ),
                tabItem(tabName = "posterior",
                        fluidRow(
                            box(width = 12,
                                h3("HMC Sampler Configuration"),
                                sliderInput("hmc_iter", "Iterations", min = 100, max = 10000, value = 1000),
                                sliderInput("hmc_thin", "Thinning", min = 1, max = 25, value = 5, round = TRUE),
                                sliderInput("hmc_eps", "Epsilon", min = 0.001, max = 0.5, value = 0.1),
                                sliderInput("hmc_leap", "Leapfrog Steps", min = 1, max = 5, value = 3, round = TRUE),
                                actionButton("hmc_run", "Run the sampler")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                h3("Sampler Results"),
                                h4("Chain results"),
                                verbatimTextOutput("hmc_res_1"),
                                h4("Chain consort"),
                                verbatimTextOutput("hmc_res_2")
                            )
                        )

                ),
                tabItem(tabName = "intervals",
                        fluidRow(
                            box(width = 12,
                                h3("Posterior Intervals"),
                                actionButton("int_show", "Show Intervals")
                            )
                        ),
                        fluidRow(
                            box(width = 12,
                                reactableOutput("int_tbl_1"),
                                plotOutput("int_plot")
                            )
                        )
                )
            )
        )
    )

    # Server ###################################################################
    server <- function(input, output) {
        output$col_txt <- renderText(paste(names(data), collapse = ", "))

        observeEvent(input$bane_create, {
            tryCatch({
                indep <- strsplit(input$idp_txt, split=",")[[1]]
                L <- lapply(strsplit(input$dep_txt, split = "\n")[[1]], function(x) strsplit(x, ":")[[1]])
                dep_names <- lapply(L, function(x) x[1])
                dep <- lapply(L, function(x) strsplit(x[2],",")[[1]])
                names(dep) <- dep_names
                mu <- rep(0.5, length(indep))
                lm <- rep(0.5, length(dep))
                bane_$create_model(indep, dep, mu, lm, data)
                output$top_plt <- renderPlot({bane_$plot_topology()})
            }, error = function(e) {
                print("Unable to parse model specification:")
                print(e)
                output$top_plt <- renderPlot({NULL})
            })

        })

        observeEvent(input$pri_reset, {
            tryCatch({
                pr <- bane_$param_names
                def_mean <- matrix(bane_$mu_pr, nrow = 1)
                def_cov <- solve(bane_$Sig_inv)
                colnames(def_mean) <- pr
                rownames(def_cov) <- pr
                colnames(def_cov) <- pr

                output$pri_mean <- renderUI({
                    shinyMatrix::matrixInput(
                        "pri_mat_mean",
                        value = def_mean,
                        cols = list(n = length(p), names = TRUE),
                        rows = list(n = 1),
                        class = "numeric"
                    )
                })

                output$pri_cov <- renderUI({
                    shinyMatrix::matrixInput(
                        "pri_mat_cov",
                        value = def_cov,
                        cols = list(n = length(p)),
                        rows = list(n = length(p)),
                        class = "numeric"
                    )
                })

            }, error = function(e) {
                print("Failed to reset page (has a model been created?)")
                print(e)
            })
        })

        observeEvent(input$pri_use, {
            tryCatch({
                bane_$mu_pr <- input$pri_mat_mean
                bane_$Sig_inv <- solve(input$pri_mat_cov)
            }, error = function(e) {
                print("Failed to assign prior values")
                print(e)
            })
        })

        observeEvent(input$beta_est, {
            tryCatch({
                y <- estimate_normal(input$beta_a, input$beta_b)
                g1 <- ggplot2::ggplot(y$sample, ggplot2::aes(x = beta)) +
                    ggplot2::geom_density() +
                    ggplot2::xlim(0, 1)
                g2 <- ggplot2::ggplot(y$sample, ggplot2::aes(x = norm)) +
                    ggplot2::geom_density()
                output$est_details <- renderText({
                    HTML(sprintf("<p>Beta: Mean %.3f, Variance %.3f</p><p> Normal: Mean %.3f, Variance %.3f</p>",
                            y$mean_b, y$var_b, y$mean_n, y$var_n))
                })
                output$plt_beta <- renderPlot({g1})
                output$plt_norm <- renderPlot({g2})
            }, error = function(e) {
                print("Estimation failed")
                print(e)
            })
        })

        observeEvent(input$mod_dtl, {
            tryCatch({
                output$mod_out <- renderPrint({bane_$print()})
            }, error = function(e) {
                print("Failed to print model details: ")
                print(e)
            })
        })

        observeEvent(input$hmc_run, {
            tryCatch({
                bane_$run_chain(
                    Iterations = input$hmc_iter,
                    Status = ceiling(input$hmc_iter/20),
                    Thinning = input$hmc_thin,
                    eps = input$hmc_eps,
                    L = input$hmc_leap
                )
                output$hmc_res_1 <- renderPrint({print(bane_$ld)})
                output$hmc_res_2 <- renderPrint({LaplacesDemon::Consort(bane_$ld)})
            }, error = function(e) {
                print("Problem with sampler:")
                print(e)
                output$hmc_res <- renderPlot({NULL})
            })
        })

        observeEvent(input$int_show, {
            tryCatch({
                post_subs <- bane_$post_subs()
                output$int_plot <- renderPlot({post_subs$ggplot})
                int_df <- dplyr::mutate_all(as.data.frame(cbind(
                    post_subs$subcohorts,
                    t(apply(post_subs$post_subs, 1, function(x) c(mean=mean(x), median = stats::median(x)))),
                    t(apply(post_subs$post_subs, 1, stats::quantile, c(0.025, 0.975)))
                )), round, 3)
                output$int_tbl_1 <- renderReactable({reactable(int_df)})
            }, error = function(e) {
                print("Error creating confidence intervals, has the sampler been run?")
                print(e)
            })
        })
    }

    shiny::runGadget(ui, server, viewer = shiny::browserViewer())
}

