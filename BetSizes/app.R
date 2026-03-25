library(shiny)
library(ggplot2)

# ------------------------------------------------------------
# Calling-function models
# ------------------------------------------------------------

call_mdf <- function(b) {
    1 / (1 + b)
}

# Linear model determined by endpoints:
# c(0) = alpha
# c(bmax) = c_at_bmax
call_linear <- function(b, alpha, c_at_bmax, bmax) {
    if (bmax <= 0) {
        return(rep(alpha, length(b)))
    }
    slope_val <- (alpha - c_at_bmax) / bmax
    pmin(pmax(alpha - slope_val * b, 0), 1)
}

call_logistic <- function(b, midpoint, steepness) {
    1 / (1 + exp(steepness * (b - midpoint)))
}

call_threshold <- function(b, threshold_b, high_call, low_call) {
    ifelse(b <= threshold_b, high_call, low_call)
}

call_polarized <- function(b, w, shape1_low, shape2_low, shape1_high, shape2_high) {
    thresh <- b / (1 + 2 * b)
    low_tail  <- 1 - pbeta(thresh, shape1 = shape1_low,  shape2 = shape2_low)
    high_tail <- 1 - pbeta(thresh, shape1 = shape1_high, shape2 = shape2_high)
    w * low_tail + (1 - w) * high_tail
}

# Pot-odds-based model:
# Caller continues if hand equity v exceeds threshold t(b)=b/(1+2b)
# and we model v ~ Beta(shape1, shape2).
call_pot_odds <- function(b, shape1, shape2) {
    thresh <- b / (1 + 2 * b)
    1 - pbeta(thresh, shape1 = shape1, shape2 = shape2)
}

# ------------------------------------------------------------
# EV function
# ------------------------------------------------------------

ev_bet <- function(b, e, cfun) {
    cvals <- cfun(b)
    (1 - cvals) * 1 + cvals * (e + (2 * e - 1) * b)
}

make_cfun <- function(model, input, bmax) {
    switch(
        model,

        linear = function(x) {
            call_linear(x, input$alpha, input$c_at_bmax, bmax)
        },

        logistic = function(x) {
            call_logistic(x, input$midpoint, input$steepness)
        },

        threshold = function(x) {
            call_threshold(x, input$threshold_b, input$high_call, input$low_call)
        },

        potodds = function(x) {
            call_pot_odds(x, input$shape1, input$shape2)
        },

        polarized = function(x) {
            call_polarized(
                x,
                input$mix_w,
                input$shape1_low, input$shape2_low,
                input$shape1_high, input$shape2_high
            )
        },

        mdf = function(x) {
            call_mdf(x)
        },

        stop("Unknown calling model")
    )
}

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------

ui <- fluidPage(
    titlePanel("Bet Sizing and Expected Value"),

    sidebarLayout(
        sidebarPanel(
            sliderInput(
                "e",
                "Hero equity if called (e):",
                min = 0, max = 1, value = 0.5, step = 0.01
            ),

            sliderInput(
                "bmax",
                "Maximum bet size shown (in pot units):",
                min = 0, max = 10, value = 2, step = 0.05
            ),

            sliderInput(
                "npts",
                "Plot resolution:",
                min = 100, max = 2000, value = 500, step = 50
            ),

            selectInput(
                "model",
                "Calling function model:",
                choices = c(
                    "Linear" = "linear",
                    "Logistic" = "logistic",
                    "Threshold" = "threshold",
                    "Pot Odds" = "potodds",
                    "Polarized (Beta Mixture)" = "polarized",
                    "Minimum Defense Frequency" = "mdf"
                ),
                selected = "linear"
            ),

            conditionalPanel(
                condition = "input.model == 'linear'",
                sliderInput(
                    "alpha",
                    "Call probability at b = 0:",
                    min = 0, max = 1, value = 1, step = 0.01
                ),
                sliderInput(
                    "c_at_bmax",
                    "Call probability at maximum bet size:",
                    min = 0, max = 1, value = 0, step = 0.01
                )
            ),

            conditionalPanel(
                condition = "input.model == 'logistic'",
                sliderInput(
                    "midpoint",
                    "Logistic midpoint:",
                    min = 0, max = 10, value = 1, step = 0.01
                ),
                sliderInput(
                    "steepness",
                    "Logistic steepness:",
                    min = 0,
                    max = 5,
                    value = 1,
                    step = 0.05
                )
            ),

            conditionalPanel(
                condition = "input.model == 'threshold'",
                sliderInput(
                    "threshold_b",
                    "Threshold bet size:",
                    min = 0, max = 10, value = 1, step = 0.01
                ),
                sliderInput(
                    "high_call",
                    "Call frequency below threshold:",
                    min = 0, max = 1, value = 1, step = 0.01
                ),
                sliderInput(
                    "low_call",
                    "Call frequency above threshold:",
                    min = 0, max = 1, value = 0, step = 0.01
                )
            ),

            conditionalPanel(
                condition = "input.model == 'potodds'",
                sliderInput("shape1", "Beta shape1:", min = 0.5, max = 10, value = 2, step = 0.1),
                sliderInput("shape2", "Beta shape2:", min = 0.5, max = 10, value = 2, step = 0.1)
            ),

            conditionalPanel(
                condition = "input.model == 'polarized'",
                sliderInput("mix_w", "Bluff/weak-hand weight:", min = 0, max = 1, value = 0.5, step = 0.01),

                tags$strong("Weak / bluff component"),
                sliderInput("shape1_low", "shape1 (weak component):", min = 0.5, max = 10, value = 2, step = 0.1),
                sliderInput("shape2_low", "shape2 (weak component):", min = 0.5, max = 10, value = 8, step = 0.1),

                tags$strong("Strong / value component"),
                sliderInput("shape1_high", "shape1 (strong component):", min = 0.5, max = 10, value = 8, step = 0.1),
                sliderInput("shape2_high", "shape2 (strong component):", min = 0.5, max = 10, value = 2, step = 0.1)
            ),

            hr(),

            checkboxInput("show_call", "Show calling function c(b)", value = TRUE),
            checkboxInput("show_ev", "Show EV(b)", value = TRUE),
            checkboxInput("show_opt", "Mark maximizing bet size", value = TRUE)
        ),

        mainPanel(
            plotOutput("mainPlot", height = "450px"),
            plotOutput("optBetByEquityPlot", height = "350px"),

            conditionalPanel(
                condition = "input.model == 'potodds' || input.model == 'polarized'",
                plotOutput("betaPlot", height = "300px")
            ),

            br(),
            tableOutput("summaryTable"),
            br(),
            verbatimTextOutput("formulaText")
        )
    )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {

    grid_data <- reactive({
        if (input$bmax <= 0) {
            b <- 0
        } else {
            b <- seq(0, input$bmax, length.out = input$npts)
        }

        cfun <- make_cfun(input$model, input, input$bmax)

        cvals <- cfun(b)
        evvals <- ev_bet(b, input$e, cfun)

        data.frame(
            b = b,
            call_prob = cvals,
            ev = evvals
        )
    })

    output$mainPlot <- renderPlot({
        df <- grid_data()

        p <- ggplot(df, aes(x = b)) +
            theme_minimal(base_size = 14) +
            labs(
                x = "Bet size b (in pot units)",
                y = NULL,
                title = "Calling Function and Bet EV"
            )

        if (input$show_call) {
            p <- p + geom_line(aes(y = call_prob, linetype = "Call probability"), linewidth = 1)
        }

        if (input$show_ev) {
            p <- p + geom_line(aes(y = ev, linetype = "Expected value"), linewidth = 1)
        }

        if (input$show_opt) {
            i_max <- which.max(df$ev)
            b_opt <- df$b[i_max]
            ev_opt <- df$ev[i_max]

            p <- p +
                annotate("vline", xintercept = b_opt, linetype = "dashed") +
                annotate("point", x = b_opt, y = ev_opt, size = 3)
        }

        p +
            scale_linetype_manual(
                values = c("Call probability" = "dashed", "Expected value" = "solid"),
                name = NULL
            )
    })

    output$optBetByEquityPlot <- renderPlot({
        e_grid <- seq(0, 1, length.out = 101)

        if (input$bmax <= 0) {
            b_grid <- 0
        } else {
            b_grid <- seq(0, input$bmax, length.out = input$npts)
        }

        cfun <- make_cfun(input$model, input, input$bmax)

        opt_b <- sapply(e_grid, function(e_val) {
            ev_vals <- ev_bet(b_grid, e_val, cfun)
            b_grid[which.max(ev_vals)]
        })

        model_label <- switch(
            input$model,
            linear = "Linear",
            logistic = "Logistic",
            threshold = "Threshold",
            potodds = "Pot Odds",
            polarized = "Polarized Mixture",
            mdf = "MDF",
            "Current Model"
        )

        df_opt <- data.frame(
            e = e_grid,
            opt_b = opt_b
        )

        ggplot(df_opt, aes(x = e, y = opt_b)) +
            geom_line(linewidth = 1) +
            theme_minimal(base_size = 14) +
            labs(
                x = "Equity if called (e)",
                y = "Optimal bet size",
                title = paste("Optimal Bet Size as a Function of Equity:", model_label)
            )
    })

    output$betaPlot <- renderPlot({
        req(input$model %in% c("potodds", "polarized"))

        v <- seq(0, 1, length.out = 500)

        b_show <- input$bmax / 2
        threshold <- b_show / (1 + 2 * b_show)

        if (input$model == "potodds") {
            df <- data.frame(
                v = v,
                density = dbeta(v, input$shape1, input$shape2)
            )

            ggplot(df, aes(x = v, y = density)) +
                geom_line(linewidth = 1.2) +
                geom_area(
                    data = subset(df, v >= threshold),
                    alpha = 0.35
                ) +
                annotate("vline", xintercept = threshold, linetype = "dashed") +
                labs(
                    x = "Villain hand equity",
                    y = "Density",
                    title = "Villain Equity Distribution: Beta Model"
                ) +
                theme_minimal()

        } else {
            low_density  <- dbeta(v, input$shape1_low,  input$shape2_low)
            high_density <- dbeta(v, input$shape1_high, input$shape2_high)
            mix_density  <- input$mix_w * low_density + (1 - input$mix_w) * high_density

            df_mix <- data.frame(v = v, density = mix_density)
            df_low <- data.frame(v = v, density = input$mix_w * low_density)
            df_high <- data.frame(v = v, density = (1 - input$mix_w) * high_density)

            ggplot() +
                geom_line(data = df_mix, aes(x = v, y = density), linewidth = 1.2) +
                geom_line(data = df_low, aes(x = v, y = density), linetype = "dotted") +
                geom_line(data = df_high, aes(x = v, y = density), linetype = "dashed") +
                geom_area(
                    data = subset(df_mix, v >= threshold),
                    aes(x = v, y = density),
                    alpha = 0.35
                ) +
                annotate("vline", xintercept = threshold, linetype = "dotdash") +
                labs(
                    x = "Villain hand equity",
                    y = "Density",
                    title = "Villain Equity Distribution: Polarized Mixture"
                ) +
                theme_minimal()
        }
    })

    output$summaryTable <- renderTable({
        df <- grid_data()
        i_max <- which.max(df$ev)

        out <- data.frame(
            Quantity = c(
                "Equity if called",
                "Optimal bet size",
                "Maximum EV",
                "Call frequency at optimal bet"
            ),
            Value = c(
                round(input$e, 4),
                round(df$b[i_max], 4),
                round(df$ev[i_max], 4),
                round(df$call_prob[i_max], 4)
            )
        )

        if (input$model == "linear") {
            slope_val <- if (input$bmax > 0) {
                (input$alpha - input$c_at_bmax) / input$bmax
            } else {
                0
            }

            out <- rbind(
                out,
                data.frame(
                    Quantity = "Implied linear slope",
                    Value = round(slope_val, 4)
                )
            )
        }

        out
    }, striped = TRUE, bordered = TRUE, spacing = "m")

    output$formulaText <- renderText({
        model_text <- switch(
            input$model,

            linear = {
                slope_val <- if (input$bmax > 0) {
                    (input$alpha - input$c_at_bmax) / input$bmax
                } else {
                    0
                }

                paste0(
                    "Current model:\n",
                    "c(0) = ", input$alpha, "\n",
                    "c(b_max) = ", input$c_at_bmax, " with b_max = ", input$bmax, "\n",
                    "Implied slope = ", round(slope_val, 6), "\n",
                    "c(b) = alpha - slope*b, clipped to [0,1]."
                )
            },

            logistic = paste0(
                "Current model:\n",
                "c(b) = 1 / (1 + exp(k*(b - m)))\n",
                "m = ", input$midpoint, "\n",
                "k = ", input$steepness
            ),

            threshold = paste0(
                "Current model:\n",
                "c(b) = ", input$high_call, " for b <= ", input$threshold_b, "\n",
                "c(b) = ", input$low_call, " for b > ", input$threshold_b
            ),

            potodds = paste0(
                "Current model:\n",
                "t(b) = b / (1 + 2b)\n",
                "v ~ Beta(", input$shape1, ", ", input$shape2, ")\n",
                "c(b) = P(v >= t(b)) = 1 - F_Beta(t(b))"
            ),

            polarized = paste0(
                "Current model:\n",
                "v ~ w*Beta(a1,b1) + (1-w)*Beta(a2,b2)\n",
                "w = ", input$mix_w, "\n",
                "Weak component: Beta(", input$shape1_low, ", ", input$shape2_low, ")\n",
                "Strong component: Beta(", input$shape1_high, ", ", input$shape2_high, ")\n",
                "t(b) = b / (1 + 2b)\n",
                "c(b) = w*(1 - F_low(t(b))) + (1-w)*(1 - F_high(t(b)))"
            ),

            mdf = paste0(
                "Current model:\n",
                "c(b) = 1 / (1 + b)\n",
                "This is the minimum defense frequency against a bet of size b into a pot of size 1."
            ),

            "Unknown model"
        )

        paste0(
            model_text,
            "\n\nEV formula:\n",
            "EV(b) = (1 - c(b))*1 + c(b)*(e + (2e - 1)b)"
        )
    })
}

shinyApp(ui = ui, server = server)