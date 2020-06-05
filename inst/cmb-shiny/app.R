library(shiny)
library(COMStuff)
library(Rcpp)

d_cmb_sample = function(x, m, p, nu, take_log, normalize) {
	n = length(x)
	ff = numeric(n)
	for (i in 1:n) {
		ff[i] = d_cmb(x[i], m, p, nu, take_log, normalize)
	}
	return(ff)
}

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
		ff = d_cmb_sample(x = 0:m, m = m, p = input$p, nu = input$nu, normalize = FALSE, take_log = FALSE)
		# Normalize outside to avoid calculating normalizing constant for each input
		ff = ff / sum(ff)
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
