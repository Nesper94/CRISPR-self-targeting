library("tidyverse")
#library("dplyr")

#Database importing
spacers = read.csv2("self-target-proteins.tsv", header = T, dec = ".", sep = "\t" )
#Object to counting (Major to minor) number of registers (species)
N_sp = spacers %>% count(Species)
N_sp = N_sp [order(-N_sp$n),]

#Min families reguster
N_fl = spacers %>% count(spacers$Family)
N_fl = N_fl [order(-N_fl$n),]

hist(N_fl$n)
mean(N_fl$n)

#Where do we cut? 
Comm =  N_fl%>% count(N_fl$n)
Comm = Comm [order(-Comm$n),]

Comm = Comm %>% rename( "Registros_con" = "N_fl$n",
                "Familias" = "n")
      #Relative frequencies of observations in the db
Fr_rel = (Comm$Familias/sum(Comm$Familias))*100 
Comm$Fr_rel = Fr_rel

              #SP= Sabiendo que el 13% de las familias sólo tiene un registrio de SS, es necesario
              #filtrarlas

#Filter of low observations families

min_obs = 3

#New SS filtered table 
spa_cut = spacers %>%
  group_by(Family) %>%
  filter(n() >= min_obs) %>%
  ungroup()


#Random sampling of family registers
spa_filt = spa_cut %>% group_split(Family) %>% lapply(function(group) sample_n(group, 3)) %>%  bind_rows()

#------------#  QUESTION #---------------#
#DEBERÍAMOS ITERAR VARIAS VECES EL PROCESO? PROS = MÁS REPRESENTATIVIDAD DE GRUPOS GRANDES
#CONTRAS= SESGO DE DATOS ÚNICOS EN GRUPOS PEQUEÑOS


#Corroboration of all Family registers = 3

corroboration_s_f =  spa_filt%>% count(spa_filt$Family)
boxplot(corroboration_s_f$n)

#Export normaliced data
write.table(spa_filt, "filtered_s-t-p.tsv", sep = "\t", row.names = FALSE)
