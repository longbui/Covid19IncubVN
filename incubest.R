
library(ggplot2)
library(dplyr)
library(lubridate)


vietnam <- readxl::read_excel("data/COVID19_VN_1311.xlsx") %>% janitor::clean_names()
vietnam$type <- ifelse(vietnam$been_to_wuhan_hubei == "Yes" | vietnam$been_to_covid_19_countries == "Yes", "imported" ,"local")
vietnam$gender[vietnam$gender=="f"] <- "Female"
vietnam$gender[vietnam$gender=="m"] <- "Male"
vietnam$gender <- Hmisc::capitalize(vietnam$gender)

i.cases <- vietnam %>% filter(!is.na(date_of_onset)) %>% filter(type=="local")  %>% filter(nationality =="vn") 
table(i.cases$gender)
prop.table(table(i.cases$gender))

incub.dat <- i.cases%>% dplyr::select(id, first_contact, first_case_of_original, last_contact,date_of_onset, date_of_confirm, date_of_public,date_of_admission,gender, age)

#preparing data input for estimation
incub <- incub.dat %>% filter(!is.na(date_of_onset)) %>% 
  dplyr::select(id,first_contact, last_contact, date_of_onset, age,gender, date_of_admission) %>%
  rename(ER = last_contact,
         EL = first_contact,
         SL = date_of_onset) %>% mutate(SR = SL) %>% dplyr::select(id, ER, EL, SR, SL, gender, age, date_of_admission) 

incub <- incub %>% mutate(E_int = difftime(ER, EL, units = "days") %>% as.numeric(),
                          S_int = difftime(SR,SL, units= "days") %>% as.numeric(),
                          type=as.numeric(S_int==0) + as.numeric(E_int==0)) %>%
                          dplyr::select(id, EL, ER, SL, SR, type, gender, age, date_of_admission) %>%
                          filter(SL > ER) %>%
                          as.data.frame()

incub <- incub %>% mutate(EL_new = difftime(EL, ER, units = "days") %>% as.integer(),
                          ER_new = difftime(ER, ER, units = "days") %>% as.integer(),
                          SL_new = difftime(SL, ER, units = "days") %>% as.integer(),
                          SR_new = difftime(SR, ER, units = "days") %>% as.integer(),
                          PL = difftime(date_of_admission, ER, units = "days") %>% as.integer()  
                          ) %>% dplyr::select (id, EL_new, ER_new, SL_new, SR_new, type, gender, age, PL) %>%
                         rename(EL = EL_new, ER = ER_new, SL = SL_new, SR = SR_new) 
incub<- incub %>% tidyr::drop_na(EL,ER,SL,SR)
             
#start rstan in modelstan.R

#Visulizing results
## mean of infection moment
Emed = tE$`50%`
incub.new <- incub %>% bind_cols(as.data.frame(Emed))

#Figure 2
over.plot <- ggplot(incub.new, aes(y=factor(id))) + 
  geom_segment(aes(x=EL, xend=ER, yend=factor(id)), 
              show.legend = TRUE) +
  geom_segment(aes(x=ER, xend=SR, yend=factor(id)), 
               show.legend = TRUE) +
  geom_point(aes(x=EL, y=factor(id), shape="EL"), size=1.5) +
  geom_point(aes(x=Emed, y=factor(id), shape="Emean"), size=1.5) +
  geom_point(aes(x=ER, y=factor(id), shape="ER"), size=1.5) +
  geom_point(aes(x=SR, y=factor(id), shape="SR"), size=1.5) +
  #geom_point(aes(x=PL, y=factor(id), colour="PL"), size=0.5) +
  geom_segment(aes(x=Emed, xend=SR, yend=factor(id)), 
               show.legend = TRUE) +
  geom_vline(xintercept = 0, color = "grey") +
  scale_x_continuous("Days from most recent possible exposure") +
  scale_y_discrete("Case") +
  scale_shape_discrete(name="",  labels =c("Earliest exposure", "Possible moment of infection",
                                            "Most recent possible exposure", "Symtomp onset")) +
  guides(shape = guide_legend(ncol=2)) +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal", legend.margin = margin(t = -0.2, unit='cm'),
        axis.ticks.y= element_blank(),
        axis.text.x=element_text(color="black")) 

#plot cumulative distribution
  estDat <- sapply(seq(0.001,0.999, by = 0.001), function(p) quantile(qweibull(p = p, shape = alpha, scale = sigma), probs = c(0.025, 0.5, 0.975)))
  colnames(estDat) <- as.character(seq(0.001,0.999, by = 0.001))
  estDat  <- as.data.frame(t(estDat))
  estDat<- cbind(rownames(estDat), data.frame(estDat, row.names=NULL))
  colnames(estDat) <- c("est","low", "median", "high")
  estDat$est = as.character(estDat$est)
  estDat$low = as.character(estDat$low)
  estDat$median = as.character(estDat$median)
  estDat$high = as.character(estDat$high)
  
  pDat <- estDat  %>% filter(est=='0.025' | est=='0.5' | est=='0.975')

  
  #estDat$est = as.numeric(estDat$est) 
 dist.Plot <- ggplot(estDat) + 
    geom_line(aes(x=as.numeric(est), y=as.numeric(median), group = 1), color = "blue", size =.5) +
    geom_line(aes(x=as.numeric(est), y=as.numeric(low), group = 1), linetype =3, color = "blue", size = .5) + 
    geom_line(aes(x=as.numeric(est), y=as.numeric(high), group = 1),  linetype =3, color = "blue", size = .5) +
    geom_ribbon(aes(ymin=as.numeric(low), ymax=as.numeric(high), x=as.numeric(est)), alpha = 0.25) +
    #geom_point(data = pDat, aes(x=as.numeric(est), y=as.numeric(median))) +
    geom_segment(dat = pDat, aes(x=as.numeric(est), y = as.numeric(low), xend = as.numeric(est), yend = as.numeric(high)), color = "red")+
    geom_label(aes(label = c("abc"), x=(as.numeric(low) + as.numeric(est))/2), size=3) +
    #coord_cartesian(xlim = c(0,1), ylim =c(0,20)) +
    scale_x_continuous("Cumulative probability distribution",breaks = seq(0, 1, by = .1)) +
    scale_y_continuous("Incubation period (in days)",breaks = seq(0, 20, by = 2)) +
    coord_flip(xlim = c(0,1), ylim =c(1,20)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) 

 #distribution uith percentiles
 ggplot(estDat) + 
   geom_line(aes(x=as.numeric(est), y=as.numeric(median), group = 1), color = "blue", size =.5) +
   geom_line(aes(x=as.numeric(est), y=as.numeric(low), group = 1), linetype =3, color = "blue", size = .5) + 
   geom_line(aes(x=as.numeric(est), y=as.numeric(high), group = 1),  linetype =3, color = "blue", size = .5) +
   geom_ribbon(aes(ymin=as.numeric(low), ymax=as.numeric(high), x=as.numeric(est)), alpha = 0.25) +
   #geom_point(data = pDat, aes(x=as.numeric(est), y=as.numeric(median))) +
   geom_segment(dat = pDat, aes(x=as.numeric(est), y = as.numeric(low), xend = as.numeric(est), yend = as.numeric(high)), color = "red")+
   geom_text(aes(x=0.975, y=7, label = "97.5^th~percentile"), parse = TRUE) +
   geom_text(aes(x=0.5, y=10, label = "50^th~percentile"), parse = TRUE) +
   geom_text(aes(x=0.025, y=5, label = "2.5^th~percentile"), parse = TRUE) +
   #coord_cartesian(xlim = c(0,1), ylim =c(0,20)) +
   scale_x_continuous("Cumulative probability distribution",breaks = seq(0, 1, by = .1)) +
   scale_y_continuous("Incubation period (in days)",breaks = seq(0, 20, by = 2)) +
   coord_flip(xlim = c(0,1), ylim =c(1,20)) +
   theme_bw() +
   theme(panel.grid.minor = element_blank())  
 
 



