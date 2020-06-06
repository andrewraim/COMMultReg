library(shiny)
library(COMMultReg)

ui = shinyUI(fluidPage(
	titlePanel("COM-Binomial vs. Binomial"),
	sidebarLayout(
		sidebarPanel(
			sliderInput("m", "m:", min = 1, max = 50, value = 20),
			sliderInput("p", "p:", min = 0.01, max = 0.99, value = 0.5),
			sliderInput("nu", "nu:", min = -3, max = 20, value = 1.0, step = 0.1)
		),
		mainPanel(
			plotOutput("cmbPlot", height = 300),
			plotOutput("binPlot", height = 300)
		)
	)
))

server = shinyServer(function(input, output) {
	output$cmbPlot = renderPlot({
		m = input$m

		# Normalize outside to avoid calculating normalizing constant for each input
		ff = numeric(m+1)
		for (x in 0:m) {
			ff[x+1] = d_cmb(x, m, input$p, input$nu, take_log = FALSE, normalize = FALSE)
		}
		ff = ff / sum(ff)

		# Plot
		barplot(ff, names.arg = 0:m, col = "blue")
		title("CMB(m, p, nu) Density")
		box()
	})
	output$binPlot = renderPlot({
		m = input$m
		ff = dbinom(0:m, prob = input$p, size = m)
		barplot(ff, names.arg = 0:m, col = "green")
		title("Binomial(m, p) Density")
		box()
	})
})

shinyApp(ui = ui, server = server)
