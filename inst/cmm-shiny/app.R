library(shiny)
library(COMMultReg)

ui = shinyUI(fluidPage(
	titlePanel("Trinomial CMM vs. Multinomial"),
	sidebarLayout(
		sidebarPanel(
			sliderInput("m", "m:",   min = 1, max = 100, value = 20, step = 1),
			sliderInput("p1", "p1 (unnormalized):", min = 0.01, max = 1, value = 1/3, step = 0.01),
			sliderInput("p2", "p2 (unnormalized):", min = 0.01, max = 1, value = 1/3, step = 0.01),
			sliderInput("p3", "p3 (unnormalized):", min = 0.01, max = 1, value = 1/3, step = 0.01),
			sliderInput("nu", "nu:", min = -3, max = 50, value = 0.8, step = 0.1)
		),
		mainPanel(
			textOutput("NormProb"),
			plotOutput("MultPlot", height = 400),
			plotOutput("CMMPlot", height = 400)
		)
	)
))

server = shinyServer(function(input, output) {
	output$NormProb = renderText({
		p = c(input$p1, input$p2, input$p3)
		p_norm = p / sum(p)
		sprintf("Normalized probabilities: (%0.06f, %0.06f, %0.06f)",
			p_norm[1], p_norm[2], p_norm[3])
	})
	output$MultPlot = renderPlot({
		# Get parameters
		p = c(input$p1, input$p2, input$p3)
		m = input$m
		nu = input$nu

		# Set up grid for plotting
		x_seq = seq(0, m)
		y_seq = seq(0, m)

		f_xy = matrix(NA, length(x_seq), length(y_seq))
		for (i in seq_along(x_seq)) {
			for (j in seq_along(y_seq)) {
				x = x_seq[i]
				y = y_seq[j]
				if (x + y <= m) {
					f_xy[i,j] = dmultinom(c(x, y, m - x - y), size = m, prob = p)
				}
			}
		}

		# Plot
		image(x_seq, y_seq, f_xy, main = "Mult(m, p)", col = rev(heat.colors(10)))
	})
	output$CMMPlot = renderPlot({
		# Get parameters
		p = c(input$p1, input$p2, input$p3)
		m = input$m
		nu = input$nu

		# Set up grid for plotting
		x = seq(0, m)
		y = seq(0, m)

		f_xy = matrix(NA, length(x), length(y))
		for (i in 1:length(x)) {
			for (j in 1:length(y)) {
				if (x[i] + y[j] <= m) {
					xval = t(c(x[i], y[j], m - x[i] - y[j]))
					f_xy[i,j] = d_cmm(xval, m = m, p = p, nu = nu)
				}
			}
		}

		# Plot
		image(x, y, f_xy, main = "CMM(m, p, nu)", col = rev(heat.colors(10)))
	})
})

shinyApp(ui = ui, server = server)
