## Instalamos paquetes
rm(list = ls())

paquetes <- c('matrixcalc', 'wordspace')

instalar <- function(paquete) {
  if (!require(paquete,character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(as.character(paquete), dependecies = TRUE, repos = "http://cran.us.r-project.org")
    library(paquete, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}


lapply(paquetes, instalar)

## Cargamos paquetes necesarios
library("matrixcalc")
library("wordspace")



#source("metadata.R")
source("utils.R")
source("00-load.R")
#source("01-prepare.R")
#source("02-clean.R")

set.seed(231)
n= 10**2

A = matrix(rnorm(n**2), ncol=n)
b = matrix(rnorm(n), ncol=1)
TOL = 10**-8

eliminacion_bloques(A,b,TOL,2)