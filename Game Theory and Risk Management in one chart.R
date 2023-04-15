#### Graphical options ####
background <- TRUE
drawArrows <- TRUE # Arrows are drawn only if background is displayed
drawPoints <- TRUE
nb_cases <- 10
sharpness <- 10000

#### Total wealth ####
x0 <- 5

#### Plot ####
r <- x0/10
p <- seq(from = 0, to = 1, by = 1/sharpness)
n <- length(p)
xNeutral <- r/p # Loss estimates for a risk-neutral entity
xAverse <- x0 * (1 - (1-r/x0)^(1/p)) # Loss estimates for a risk-averse entity

plot(x = p, y = xAverse, type = "l", col = "red",
     xlim = c(0, 1.01), ylim = c(0, 1.1 * x0),
     xlab = "", ylab = "",
     xaxs = "i", yaxs = "i",
     main = "Iso-risk curves for risk-neutral and risk-averse entities")
title(sub = "Risk cost r is defined as the difference between the current wealth \n and the certainty equivalent in a risky environment", cex.sub = 0.8)

for (r in seq(from = x0/nb_cases, to = x0, by = x0/nb_cases)){
  xNeutral <- r/p
  xAverse <- x0 * (1 - (1-r/x0)^(1/p))
  lines(x = p, y = xAverse, col = "red")
  lines(x = p, y = xNeutral, col = "blue")
  text(x = p[n] - 0.03, y = xAverse[n] - x0/125, label = paste("r =", round(r, digits = 2)), cex = 0.8)
}

if (background){
  rect(xleft = 0, ybottom = 0, xright = 0.4, ytop = 0.3*x0,
       col = rgb(red = 30.59/100, green = 89.41/100, blue = 30.59/100, alpha = 0.3), border = "transparent")
  text(x = 0.175, y = 1 * x0/5, label = "Low risk: inexpensive insurance", col = "darkgreen", cex = 0.8)
  
  rect(xleft = 0.4, ybottom = 0, xright = 1.01, ytop = 0.3*x0,
       col = rgb(red = 25.1/100, green = 72.55/100, blue = 100/100, alpha = 0.3), border = "transparent")
  text(x = 0.707, y = 0.22 * x0/5, label = "Frequency losses:", col = "blue", cex = 0.8)
  text(x = 0.8, y = 0.08 * x0/5, label = "self-insured or inexpensive insurance", col = "blue", cex = 0.8)
  
  rect(xleft = 0, ybottom = 0.3*x0, xright = 0.4, ytop = x0,
       col = rgb(red = 25.1/100, green = 72.55/100, blue = 100/100, alpha = 0.3), border = "transparent")
  text(x = 0.09, y = 1.79 * x0/5, label = "Severity losses:", col = "blue", cex = 0.8)
  text(x = 0.174, y = 1.65 * x0/5, label = "relatively inexpensive insurance", col = "blue", cex = 0.8)
  
  rect(xleft = 0.4, ybottom = 0.3*x0, xright = 1.01, ytop = x0,
       col = rgb(red = 100/100, green = 25.1/100, blue = 25.1/100, alpha = 0.3), border = "transparent")
  text(x = 0.76, y = 1.65 * x0/5, label = "High risk: expensive insurance", col = "darkred", cex = 0.8)
  
  title(xlab = "Loss probability", line = 2, cex.lab = 1)
  title(ylab = "Loss estimate", line = 2, cex.lab = 1)
}else{
  title(xlab = "Loss probability", line = -1, cex.lab = 1)
  title(ylab = "Loss estimate", line = -1, cex.lab = 1)
}

if ((background) && (drawArrows)){
  arrows(x0 = 0.55, y0 = 2.5 * x0/5, x1 = 0.3, y1 = 2.5 * x0/5, length = 0.1, col = "gray22")
  text(x = 0.46, y = 2.6 * x0/5, label = "Prevention", col = "gray22", cex = 0.8)
  arrows(x0 = 0.55, y0 = 2.5 * x0/5, x1 = 0.55, y1 = 1 * x0/5, length = 0.1, col = "gray22")
  text(x = 0.61, y = 2 * x0/5, label = "Protection", col = "gray22", cex = 0.8)
}

if(drawPoints){
  points(x = 0.4, y = 0.4 * x0, pch = 16, col = "gray22")
  arrows(x0 = 0.4, y0 = 0, x1 = 0.4, y1 = 0.4 * x0, length = 0, col = "gray22", lty = 3)
  arrows(x0 = 0, y0 = 0.4 * x0, x1 = 0.4, y1 = 0.4 * x0, length = 0, col = "gray22", lty = 3)
  text(x = 0.4, y = 0.425 * x0, label = "Warehouse", cex = 0.8, col = "gray22")
}

text(x = 0.083, y = 1.02 * x0, label = "Current wealth", cex = 0.8)
legend("bottomleft", legend = c("Risk-neutral", "Risk-averse"), col = c("blue", "red"), pch = c("_", "_"))