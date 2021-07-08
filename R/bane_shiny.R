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
    # UI #######################################################################
    ui <- dashboardPage(
        dashboardHeader(title = "BaNE Model Implementation"),
        dashboardSidebar(
            sidebarMenu(
                menuItem("Model Specification", tabName = "model", icon = icon("th")),
                menuItem("Prior Specification", tabName = "prior", icon = icon("th")),
                menuItem("Model Details", tabName = "details", icon = icon("th")),
                menuItem("Posterior Sampling", tabName = "posterior", icon = icon("th"))
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
                            box(width = 12)
                        )

                ),
                tabItem(tabName = "details",
                        fluidRow(
                            box(width = 12)
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
                                plotOutput("hmc_res")
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
                bane_ <- Bane$new(indep, dep, mu, lm, data)
                output$top_plt <- renderPlot({bane_$plot_topology()})
            }, error = function(e) {
                print("Unable to parse model specification:")
                print(e)
                output$top_plt <- renderPlot({NULL})
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
                output$hmc_res <- renderPlot({bane_$post_subs()$ggplot})
            }, error = function(e) {
                print("Problem with sampler:")
                print(e)
                output$hmc_res <- renderPlot({NULL})
            })
        })
    }

    shiny::runGadget(ui, server, viewer = shiny::dialogViewer("BaNE", width = 800))
}
