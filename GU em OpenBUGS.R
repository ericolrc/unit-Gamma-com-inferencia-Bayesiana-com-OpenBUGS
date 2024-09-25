
# Carregar as bibliotecas necessárias
library(R2OpenBUGS)  # Para utilizar o OpenBUGS
library(coda)        # Para análise de cadeias MCMC
library(lattice)     # Para visualização de dados

# Parâmetros simulados: mu = 0.9; phi = 5; n = 30
dados = list(
  x = c(0.8740038, 0.8347652, 0.8273605, 0.9203040, 0.9515342,
        0.9420236, 0.9171107, 0.8562678, 0.8505088, 0.9592614,
        0.9029226, 0.8552615, 0.9200630, 0.9507948, 0.9009879,
        0.9535192, 0.8828526, 0.9527309, 0.9114653, 0.7365916,
        0.8353674, 0.9303087, 0.9133851, 0.9085975, 0.8658272,
        0.9032948, 0.9334139, 0.9110197, 0.9303046, 0.8960092),
  N = 30  # Número de observações
)

# Definindo os parâmetros a serem monitorados
params = c("mu", "phi")

# Função para definir os valores iniciais dos parâmetros
inits <- function() {
  list(mu = 0.2, phi = 2)
}

# Modelo a ser estimado
modelo <- function() {
  for (i in 1:N) {
    dummy[i] <- 0    
    dummy[i] ~ dpois(logLike[i])  # Distribuição Poisson para a variável dummy
    
    # Cálculo da log-verossimilhança
    logLike[i] <- -phi * log(beta) + loggam(phi) - 
      (beta - 1) * log(x[i]) - (phi - 1) * log(-log(x[i])) + C
  }
  
  # Cálculo do beta
  beta <- (pow(mu, 1/phi) / (1 - pow(mu, 1/phi)))
  
  # Priors para mu e phi
  mu ~ dbeta(1, 1)              # Distribuição Beta para mu
  phi ~ dgamma(0.1, 0.1)        # Distribuição Gamma para phi
  C <- 100000                   # Constante
}

# Executar o modelo utilizando o OpenBUGS
saida = bugs(data = dados,
             inits = inits,
             parameters.to.save = params,
             model.file = modelo,
             n.chains = 3,
             n.iter = 2500,
             n.burnin = 1000,
             n.thin = 10,
             codaPkg = TRUE,
             debug = TRUE,
             working.directory = NULL)

# Ler os resultados da análise
line_coda = read.bugs(saida)

# Converter os resultados para um formato MCMC
line_coda2 = as.mcmc(as.matrix(line_coda[, 2:3]))

# Exibir as dimensões dos resultados
dim(line_coda2)

# Plotar as cadeias MCMC
plot(line_coda2, col = "gray")

# Resumo dos parâmetros mu e phi
summary(line_coda2[, 1])  # Resumo para mu
summary(line_coda2[, 2])  # Resumo para phi

# Intervalo de densidade posterior (HPD)
HPDinterval(line_coda2)

# Gráfico da densidade dos parâmetros
densityplot(line_coda2)

# Gráfico da função de autocorrelação
acfplot(line_coda2)
