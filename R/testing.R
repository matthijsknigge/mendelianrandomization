# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


data("BMI")
data("alzheimer")
data("cholesterol")

int <- Reduce(intersect, list(BMI$SNP, cholesterol$SNP))

outcome <- BMI[BMI$SNP %in% int, ]
exposure <- cholesterol[cholesterol$SNP %in% int, ]

outcome <- outcome[! duplicated(outcome$SNP), ]
exposure <- exposure[! duplicated(exposure$SNP), ]


iv <- mr.wald.ratio(By = outcome$beta, Bx = exposure$beta, By.se = outcome$se, Bx.se = exposure$se)
ivw <- mr.inverse.variance.weighted.method(By = outcome$beta, Bx = exposure$beta, By.se = outcome$se, Bx.se = exposure$se)
egger <- mr.egger.method(By = outcome$beta, Bx = exposure$beta, By.se = outcome$se, Bx.se = exposure$se)


