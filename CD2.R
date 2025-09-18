# Set up for individual plot saving
setwd("C:/Users/subha/Downloads")

# Plot 1: Cumulative event arrivals
pdf("cat_cumulative_arrivals.pdf", width=6, height=4)
plot(arrival_times, 1:n_events, type = "s", lwd = 2,
     xlab = "Years since 1980", ylab = "Cumulative events",
     main = "Cumulative Event Arrivals")
abline(a = 0, b = hpp_rate, col = "red", lty = 2, lwd = 2)
legend("topleft", c("Observed", "HPP expected"), 
       col = c("black", "red"), lty = c(1, 2), lwd = 2)
dev.off()

# Plot 2: Inter-arrival time distribution
pdf("cat_interarrival_dist.pdf", width=6, height=4)
hist(inter_arrivals, breaks = 30, freq = FALSE,
     xlab = "Inter-arrival time (years)", 
     main = "Inter-arrival Distribution",
     col = "lightgray")
curve(dexp(x, rate = hpp_rate), add = TRUE, col = "red", lwd = 2)
dev.off()

# Plot 3: Monthly event counts
pdf("cat_monthly_counts.pdf", width=6, height=4)
hist(monthly_counts, breaks = seq(-0.5, max(monthly_counts) + 0.5, 1),
     xlab = "Events per month", main = "Monthly Event Distribution",
     col = "lightblue")
abline(v = mean_count, col = "red", lwd = 2, lty = 2)
abline(v = median(monthly_counts), col = "blue", lwd = 2, lty = 3)
legend("topright", c("Mean", "Median"), col = c("red", "blue"), lty = c(2, 3))
dev.off()

# Plot 4: Autocorrelation function
pdf("cat_autocorrelation.pdf", width=6, height=4)
acf(monthly_counts, main = "Autocorrelation of Monthly Counts", lag.max = 12)
dev.off()

# Plot 5: Q-Q plot for exponential inter-arrivals
pdf("cat_qq_plot.pdf", width=6, height=4)
qqplot(qexp(ppoints(length(inter_arrivals)), rate = hpp_rate),
       inter_arrivals,
       main = "Q-Q Plot: Exponential Fit",
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles")
abline(0, 1, col = "red", lty = 2, lwd = 2)
dev.off()

# Plot 6: Model comparison
pdf("cat_model_comparison.pdf", width=6, height=4)
comparison_data <- matrix(c(AIC_hpp, AIC_mmpp, BIC_hpp, BIC_mmpp), 
                          nrow = 2, byrow = TRUE)
barplot(comparison_data, beside = TRUE,
        names.arg = c("HPP", "MMPP"),
        main = "Model Comparison",
        ylab = "Information Criterion",
        col = c("lightcoral", "lightblue"),
        legend.text = c("AIC", "BIC"),
        args.legend = list(x = "topright"))
dev.off()