
# Organizando o banco -----------------------------------------------------


library(dplyr)
library(MVar.pt)
library(Kira)
library(MVN)

dados   <- read.csv('ibge.csv', encoding ='UTF-8')
dados21 <- dados[dados$Ano == 2021, ]
row.names(dados21) <- 1:nrow(dados21)
colunas <- c('vlr_agro', 'vlr_ind', 'vlr_ser', 'vlr_pub',
             'vlr_imp', 'area','tot_pop','pib_capta', 'tot_dom_par','tot_dom_col')
dados2 <- dados21[,c(colunas, 'Nom_reg') ]

result = mvn(data = dados2, subset = "Nom_reg", mvnTest = "hz",
             univariateTest = "AD", univariatePlot = "histogram",
             multivariatePlot = "qq", multivariateOutlierMethod = "adj",
             showOutliers = TRUE, showNewData = TRUE)


co_com <- dados2[result$multivariateOutliers$`Centro-oeste`$Observation,]
nd_com <- dados2[result$multivariateOutliers$Nordeste$Observation,]
nt_com <- dados2[result$multivariateOutliers$Norte$Observation,]
sd_com <- dados2[result$multivariateOutliers$Sudeste$Observation,]
sl_com <- dados2[result$multivariateOutliers$Sul$Observation,]
com_out <- rbind(co_com, nd_com, nt_com, sd_com, sl_com)

co_sem  <- cbind(result$newData$`Centro-oeste`, Nom_reg = 'Centro-oeste')
nd_sem  <- cbind(result$newData$Nordeste, Nom_reg = 'Nordeste')
nt_sem  <- cbind(result$newData$Norte, Nom_reg = 'Norte')
sd_sem  <- cbind(result$newData$Sudeste, Nom_reg = 'Sudeste')
sl_sem  <- cbind(result$newData$Sul, Nom_reg = 'Sul')
sem_out <- rbind(co_sem,nd_sem,nt_sem,sd_sem,sl_sem)



# MANOVA ------------------------------------------------------------------

# Para o banco de outliers
model_com <- manova(as.matrix(com_out[,-ncol(com_out)]) ~ com_out$Nom_reg, data = com_out) 
model_sem <- manova(as.matrix(sem_out[,-ncol(sem_out)]) ~ sem_out$Nom_reg, data = sem_out) 

options(scipen = 999) # tira a notacao cientifica para os numeros muito pequenos

# Pelo menos uma média nas variáveis é diferente pro tratamento
summary(model_com, test = "Wilks", intercept = FALSE)
summary(model_com, test = "Pillai", intercept = FALSE)
summary(model_com, test = "Hotelling-Lawley", intercept = FALSE)
summary(model_com, test = "Roy", intercept = FALSE)

# ------------------

summary(model_sem, test = "Wilks", intercept = FALSE)
summary(model_sem, test = "Pillai", intercept = FALSE)
summary(model_sem, test = "Hotelling-Lawley", intercept = FALSE)
summary(model_sem, test = "Roy", intercept = FALSE)


# Intervalos de confiança -------------------------------------------------

bancos_sem <- list(co_sem, nd_sem, nt_sem, sd_sem, sl_sem)
bancos_com <- list(co_com, nd_com, nt_com, sd_com, sl_com)



IC <- function(bancos,interval = c('B','T'),  sig = 0.05)
{
  intervalos <- matrix(NA,1,2)
  for(i in 1:length(bancos))
  {
    data <- bancos[[i]]
    if(!is.null(data$Nom_reg)){data$Nom_reg <- NULL}
    n   = nrow(data)        # numero de observacoes 
    sigma  = cov(data)      # matriz de covariancia
    mi.amo = colMeans(data) # vetor media amostral
    
    ll = diag(rep(1,ncol(data))) # vetores l's
    p = ncol(data)     # grau de liberdade
    v = nrow(data) - 1 # grau de liberdade
    m = ncol(data)     # numero de variaveis
    
    IC <- matrix(NA,ncol = 2,nrow = ncol(data)) # matriz com os IC's
    colnames(IC) <- c("Lim.Inferior","Lim.Superior")
    rownames(IC) <- colnames(data)
    
    if(interval == "T") { # T^2 Hotelling
      TB = v * p / (v - p + 1) * qf(1 - sig, df1 = p, df2 = v - p + 1, ncp = 0)
      for(j in 1:ncol(data)) {
        Erro.Media <- sqrt(TB) * sqrt((t(ll[,j]) %*% sigma %*% ll[,j]) / n)
        IC[j,1] <- t(ll[,j]) %*% mi.amo - Erro.Media
        IC[j,2] <- t(ll[,j]) %*% mi.amo + Erro.Media
      }
    }
    
    if(interval == "B") { # Bonferroni
      for(j in 1:ncol(data)) {
        Erro.Media <- qt(1 - sig/(2*m), df = v ) * sqrt((t(ll[,j]) %*% sigma %*% ll[,j]) / n)
        IC[j,1] <- t(ll[,j]) %*% mi.amo - Erro.Media
        IC[j,2] <- t(ll[,j]) %*% mi.amo + Erro.Media
      }
    }
    

    
    intervalos <-  round(as.data.frame(na.omit(rbind(intervalos, IC))),2)
  }
  return(intervalos)
}
  
intervalos_B_com <- IC(bancos_com,'B')
intervalos_B_sem <- IC(bancos_sem,'B')
n = nrow(intervalos_B_sem)/5
regioes = c(rep('Centro-Oeste',n),rep('Nordeste',n),rep('Norte',n),rep('Sudeste',n),rep('Sul',n))
intervalos_B_com$regiao = regioes
intervalos_B_sem$regiao = regioes
intervalos_B_com$variavel = colunas
intervalos_B_sem$variavel = colunas
intervalos_B_com$outlier = 'Com'
intervalos_B_sem$outlier = 'Sem'



intervalos <- rbind(intervalos_B_com,intervalos_B_sem)
intervalos$x_pos = paste0(intervalos$regiao,ifelse(intervalos$outlier=='Sem', ' ', ''))
library(ggplot2)

plot.ic <- function(y, ylab, x='x_pos', cor, paleta){
  ax <- ggplot(intervalos[intervalos$variavel == y,], aes(x = x_pos, colour = outlier)) +
    geom_errorbar(aes(ymin = Lim.Inferior, ymax = Lim.Superior), width = 0.2, linewidth = 1) +
    labs(x = 'Região', y = ylab,colour ='Outliers') +
    theme_minimal() +
    coord_flip() +
    scale_color_manual(values = c('#004586', '#FF420E')) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')
  ggsave(paste0('img/', y,'.jpeg'))
  return(ax)
}
plot.ic('vlr_agro', 'Valor Agro (R$1.000)')
plot.ic('vlr_ind', 'Valor Indústria (R$1.000)')
plot.ic('vlr_ser', 'Valor Serviços (R$1.000)')
plot.ic('vlr_pub', 'Valor Serviço Público (R$1.000)')
plot.ic('vlr_imp', 'Valor Impostos (R$1.000)')
plot.ic('area', 'Área (m²)')
plot.ic('tot_pop', 'Total da população')
plot.ic('pib_capta', 'PIB per capita (R$1.000)')
plot.ic('tot_dom_par', 'Total domicílios')
plot.ic('tot_dom_col', 'Total de domicílios')

# geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
  

# Pra separar por banco 
# n          <- nrow(instervalos_B_com)/length(bancos_com)
# IC_co <- intervalo[1:n          ,]; rownames(IC_co) <- colunas
# IC_ne <- intervalo[(n+1):(2*n)  ,]; rownames(IC_ne) <- colunas
# IC_n  <- intervalo[(2*n+1):(3*n),]; rownames(IC_n ) <- colunas
# IC_se <- intervalo[(3*n+1):(4*n),]; rownames(IC_se) <- colunas
# IC_s  <- intervalo[(4*n+1):(5*n),]; rownames(IC_s ) <- colunas

