#' @import ggplot2
NULL
#' Draws percentage plot for different type of values
#'
#' Draws plot that graphically shows the percentages of the missing, censored and observed data. It supports
#' generating plots for all major types of censoring including left, right and interval censoring.
#'
#' @param data.indicator matrix including the data type indicators of the original data.
#' @param title title of the generated plot, default is set to "Percentages of different data type".
#'
#' @details The function draws the plot that graphically shows the percentages of the missing, censored and observed
#' data in the dataset. \code{data.indicator} should be a matrix containing the data type indicators as generated in the
#' data preparation step. 0 for missing values, 1 for observed values, and 2 for right censored values, 3 for left censored values,
#' and 4 for interval censored values. \code{title} is the title
#' of the generated plot.
#'
#' @examples
#' data.ind <- simulated.dat[[2]]
#' visual.plot(data.ind)
#'
#' @return The plot that shows the details of the different type of data in the dataset.
#'
#' @export
visual.plot <- function(
  data.indicator,
  title = "Percentages of different data type") {

  c.names <- colnames(data.indicator) # get column names of the data
  n <- nrow(data.indicator)
  p <- ncol(data.indicator)

  data.t <- matrix(NA, nrow = p, ncol = 5)
  colnames(data.t) <- c("Observed", "Missing", "Left_Censored", "Right_Censored", "Interval_Censored")

  for (i in 1:p) {

    summ <- table(data.indicator[, i])
    names.summ <- names(summ)

    if ("0" %in% names.summ) {
      data.t[i, "Missing"] <- summ[names.summ == 0]
    } else {
      data.t[i, "Missing"] <- 0
    }

    if ("1" %in% names.summ) {
      data.t[i, "Observed"] <- summ[names.summ == 1]
    } else {
      data.t[i, "Observed"] <- 0
    }

    if ("2" %in% names.summ) {
      data.t[i, "Right_Censored"] <- summ[names.summ == 2]
    } else {
      data.t[i, "Right_Censored"] <- 0
    }

    if ("3" %in% names.summ) {
      data.t[i, "Left_Censored"] <- summ[names.summ == 3]
    } else {
      data.t[i, "Left_Censored"] <- 0
    }

    if ("4" %in% names.summ) {
      data.t[i, "Interval_Censored"] <- summ[names.summ == 4]
    } else {
      data.t[i, "Interval_Censored"] <- 0
    }

  }

  new.data <- data.frame(id = c.names, data.t)
  new.data$Observed <- new.data$Observed/n * 100
  new.data$Missing <- new.data$Missing/n * 100
  new.data$Right_Censored <- new.data$Right_Censored/n * 100
  new.data$Left_Censored <- new.data$Left_Censored/n * 100
  new.data$Interval_Censored <- new.data$Interval_Censored/n * 100

  percent.mat <- melt(new.data, id.vars = "id")

  ### remove 0 count
  percent.mat <- subset(percent.mat, percent.mat$value != 0.00)
  colnames(percent.mat)[2] <- "Variable"

  ggplot(percent.mat) +
    geom_bar(aes(x = factor(.data$id, levels = c.names), y = .data$value, fill = .data$Variable),
             stat = "identity") +
    coord_flip() +
    labs(x = "Variable", y = "Percentage (%)") +
    scale_fill_manual(breaks = c("Missing", "Observed", "Right_Censored", "Left_Censored", "Interval_Censored"),
                      values = c("#00A9FF", "#00BE67", "#DB72FB", "#FFD64C", "#F8766D")) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}
