rm(list=ls())

# setup
need <- c('tidyverse','modelsummary','haven','fixest',"car") # list packages needed
have <- need %in% rownames(installed.packages()) # checks packages you have
if(any(!have)) install.packages(need[!have]) # install missing packages
invisible(lapply(need, library, character.only=T))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data <- read_dta("PopulismTerrorism.dta")
data <- as_factor(data)

#############################twfe########################################
load("cleaned/data.RData")
#total
outcome <- "PopulistTotal~"
homeattacks <- "HomeAttacksY1+HomeAttacksY2+HomeAttacksY3+HomeAttacksY4"
fe <- "|Country+Year"
reg_total_attack_none <- feols(as.formula(paste0(outcome,homeattacks,fe)), weights=data$weight,data)
controls <- "+GrowthY1+GrowthY2+GrowthY3+UnemploymentY1+UnemploymentY2+UnemploymentY3+
Export_penetrationY1+Export_penetrationY2+Export_penetrationY3+
Import_penetrationY1+Import_penetrationY2+Import_penetrationY3+
Immigration_percentY1+Immigration_percentY2+Immigration_percentY3+
Emigration_percentY1+Emigration_percentY2+Emigration_percentY3"
reg_total_attack_control <- feols(as.formula(paste0(outcome,homeattacks,controls,fe)), weights=data$weight,data)
geoattacks <- "+geoAttacksY1+geoAttacksY2+geoAttacksY3+geoAttacksY4"
reg_total_attack_geo <- feols(as.formula(paste0(outcome,homeattacks,geoattacks,controls,fe)), weights=data$weight,data)
culattacks <- "+culAttacksY1+culAttacksY2+culAttacksY3+culAttacksY4"
reg_total_attack_cul <- feols(as.formula(paste0(outcome,homeattacks,culattacks,controls,fe)), weights=data$weight,data)

homecasualty <- "HomencasY1+HomencasY2+HomencasY3+HomencasY4"
geocasualty <- "+geoncasY1+geoncasY2+geoncasY3+geoncasY4"
culcasualty <- "+culncasY1+culncasY2+culncasY3+culncasY4"
reg_total_casualty_none <- feols(as.formula(paste0(outcome,homecasualty,fe)), weights=data$weight,data)
reg_total_casualty_control <- feols(as.formula(paste0(outcome,homecasualty,controls,fe)), weights=data$weight,data)
reg_total_casualty_geo <- feols(as.formula(paste0(outcome,homecasualty,geocasualty,controls,fe)), weights=data$weight,data)
reg_total_casualty_cul <- feols(as.formula(paste0(outcome,homecasualty,culcasualty,controls,fe)), weights=data$weight,data)

coef_map <- c("HomeAttacksY1"="HomeTerrorism$_{t-1}$",
              "HomeAttacksY2"="HomeTerrorism$_{t-2}$",
              "HomeAttacksY3"="HomeTerrorism$_{t-3}$",
              "HomeAttacksY4"="HomeTerrorism$_{t-4}$",
              "HomencasY1"="HomeTerrorism$_{t-1}$",
              "HomencasY2"="HomeTerrorism$_{t-2}$",
              "HomencasY3"="HomeTerrorism$_{t-3}$",
              "HomencasY4"="HomeTerrorism$_{t-4}$",
              "geoAttacksY1"="ForeignTerrorism$_{t-1}$",
              "geoAttacksY2"="ForeignTerrorism$_{t-2}$",
              "geoAttacksY3"="ForeignTerrorism$_{t-3}$",
              "geoAttacksY4"="ForeignTerrorism$_{t-4}$",
              "geoncasY1"="ForeignTerrorism$_{t-1}$",
              "geoncasY2"="ForeignTerrorism$_{t-2}$",
              "geoncasY3"="ForeignTerrorism$_{t-3}$",
              "geoncasY4"="ForeignTerrorism$_{t-4}$",
              "culAttacksY1"="ForeignTerrorism$_{t-1}$",
              "culAttacksY2"="ForeignTerrorism$_{t-2}$",
              "culAttacksY3"="ForeignTerrorism$_{t-3}$",
              "culAttacksY4"="ForeignTerrorism$_{t-4}$",
              "culncasY1"="ForeignTerrorism$_{t-1}$",
              "culncasY2"="ForeignTerrorism$_{t-2}$",
              "culncasY3"="ForeignTerrorism$_{t-3}$",
              "culncasY4"="ForeignTerrorism$_{t-4}$")
twfe_total <- list(reg_total_attack_none, reg_total_attack_control, 
                   reg_total_attack_geo, reg_total_attack_cul,
                   reg_total_casualty_none, reg_total_casualty_control,
                   reg_total_casualty_geo, reg_total_casualty_cul)
for (i in 1:length(twfe_total)){
  twfe_total[[i]]$df.residual <- parameters::degrees_of_freedom(twfe_total[[i]])
  home <- names(coef(twfe_total[[i]]))[grepl("Home",names(coef(twfe_total[[i]])))]
  twfe_total[[i]]$homeF <- lht(twfe_total[[i]],test='F',
                               c(paste0(home,"=0")))
  foreign <- character()
  foreign <- names(coef(twfe_total[[i]]))[grepl("geo|cul",names(coef(twfe_total[[i]])))]
  if (length(foreign)>0){
    twfe_total[[i]]$foreignF <- lht(twfe_total[[i]],test='F',
                                    c(paste0(foreign,"=0"))
    )
  }
}

gen_stars <- function(p){
  if (p<0.001){
    return("***")
  } else if (p<0.01){
    return("**")
  } else if (p<0.05){
    return("*")
  } else if (p<0.1){
    return("+")
  } else {
    return("")
  }
}

glance_custom.fixest <- function(x, ...) {
  data.frame(
    "controls" = ifelse("GrowthY1" %in% names(coef(x)),
                        "Yes","No"),
    'vcov.type'="Clustered",
    'weighted'= "Yes",
    'CountryFE' = "Yes",
    'YearFE' = "Yes",
    'HomeF' = paste0(round(x$homeF[2,3],2),gen_stars(x$homeF[2,4]),
                   " [",x$homeF[2,2],",",x$homeF[2,1],"]"),
    "ForeignF" = if(!"foreignF" %in% names(x)){
      ""
    } else {
      paste0(round(x$foreignF[2,3],2),gen_stars(x$foreignF[2,4]),
             " [",x$foreignF[2,2],",",x$foreignF[2,1],"]")
    }
  )
}

gof_map <- list(list("raw"="controls", clean="Controls", fmt=NULL),
                list("raw"="weighted", clean="Weighted", fmt=NULL),
                list("raw"="HomeF", "clean" = "$F_{home}$", "fmt"=NULL),
                list("raw"="ForeignF", "clean" = "$F_{foreign}$", "fmt"=NULL),
                list("raw" = "nobs", "clean" = "$N$", "fmt" = 0),
                list("raw" = "adj.r.squared", "clean" = "$R^2$ Adj.", "fmt" = 3),
                list("raw" = "r2.within.adjusted", "clean" = "$R^2$ Within Adj.", "fmt" = 3),
                list('raw' = 'vcov.type', 'clean' = 'Std.Errors', 'fmt' = NULL),
                list("raw"="CountryFE", "clean" = "Country FE", "fmt" = NULL),
                list("raw"="YearFE", "clean" = "Year FE", "fmt" = NULL))

#table_twfe_total <- 
  msummary(twfe_total, gof_map = gof_map,
         coef_map=coef_map, stars=T) |>
  add_header_above(c(" "=1,"$c=0$"=2,"$c=c^{geo}$"=1,"$c=c^{cul}$",
                     "$c=0$"=2,"$c=c^{geo}$"=1,"$c=c^{cul}$")) |>
  add_header_above(c(" "=1,"Total Attacks"=4,"Total Casualties"=4))

save(twfe_total, gof_map, coef_map, glance_custom.fixest, 
     gen_stars, file="twfe_total.RData")

#right
outcome <- "PopulistRight~"
homeattacks <- "HomeleftproxyY1+HomeleftproxyY2+HomeleftproxyY3+HomeleftproxyY4+
HomerightproxyY1+HomerightproxyY2+HomerightproxyY3+HomerightproxyY4+
HomeIslClaimY1+HomeIslClaimY2+HomeIslClaimY3+HomeIslClaimY4"
fe <- "|Country+Year"
reg_right_attack_none <- feols(as.formula(paste0(outcome,homeattacks,fe)),data)
controls <- "+GrowthY1+GrowthY2+GrowthY3+UnemploymentY1+UnemploymentY2+UnemploymentY3+
Export_penetrationY1+Export_penetrationY2+Export_penetrationY3+
Import_penetrationY1+Import_penetrationY2+Import_penetrationY3+
Immigration_percentY1+Immigration_percentY2+Immigration_percentY3+
Emigration_percentY1+Emigration_percentY2+Emigration_percentY3"
reg_right_attack_control <- feols(as.formula(paste0(outcome,homeattacks,controls,fe)),data)
geoattacks <- "+geoleftproxyY1+geoleftproxyY2+geoleftproxyY3+geoleftproxyY4+
georightproxyY1+georightproxyY2+georightproxyY3+georightproxyY4+
geoIslClaimY1+geoIslClaimY2+geoIslClaimY3+geoIslClaimY4"
reg_right_attack_geo <- feols(as.formula(paste0(outcome,homeattacks,geoattacks,controls,fe)),data)
culattacks <- "+culleftproxyY1+culleftproxyY2+culleftproxyY3+culleftproxyY4+
culrightproxyY1+culrightproxyY2+culrightproxyY3+culrightproxyY4+
culIslClaimY1+culIslClaimY2+culIslClaimY3+culIslClaimY4"
reg_right_attack_cul <- feols(as.formula(paste0(outcome,homeattacks,culattacks,controls,fe)),data)

homecasualty <- "Homencas_leftY1+Homencas_leftY2+Homencas_leftY3+Homencas_leftY4+
Homencas_rightY1+Homencas_rightY2+Homencas_rightY3+Homencas_rightY4+
Homencas_IslClY1+Homencas_IslClY2+Homencas_IslClY3+Homencas_IslClY4"
geocasualty <- "+geoncas_leftY1+geoncas_leftY2+geoncas_leftY3+geoncas_leftY4+
geoncas_rightY1+geoncas_rightY2+geoncas_rightY3+geoncas_rightY4+
geoncas_IslClY1+geoncas_IslClY2+geoncas_IslClY3+geoncas_IslClY4"
culcasualty <- "+culncas_leftY1+culncas_leftY2+culncas_leftY3+culncas_leftY4+
culncas_rightY1+culncas_rightY2+culncas_rightY3+culncas_rightY4+
culncas_IslClY1+culncas_IslClY2+culncas_IslClY3+culncas_IslClY4"
reg_right_casualty_none <- feols(as.formula(paste0(outcome,homecasualty,fe)),data)
reg_right_casualty_control <- feols(as.formula(paste0(outcome,homecasualty,controls,fe)),data)
reg_right_casualty_geo <- feols(as.formula(paste0(outcome,homecasualty,geoattacks,controls,fe)),data)
reg_right_casualty_cul <- feols(as.formula(paste0(outcome,homecasualty,culattacks,controls,fe)),data)
twfe_right <- list(reg_right_attack_control, 
                   reg_right_attack_geo, reg_right_attack_cul,
                   reg_right_casualty_control,
                   reg_right_casualty_geo, reg_right_casualty_cul)
for (i in 1:length(twfe_right)){
  twfe_right[[i]]$df.residual <- parameters::degrees_of_freedom(twfe_right[[i]])
  home <- names(coef(twfe_right[[i]]))[grepl("Home",names(coef(twfe_right[[i]])))]
  twfe_right[[i]]$homeF <- lht(twfe_right[[i]],test='F',
                               c(paste0(home,"=0")))
  foreign <- character()
  foreign <- names(coef(twfe_right[[i]]))[grepl("geo|cul",names(coef(twfe_right[[i]])))]
  if (length(foreign)>0){
    twfe_right[[i]]$foreignF <- lht(twfe_right[[i]],test='F',
                                    c(paste0(foreign,"=0"))
    )
  }
}
coef_map <- c("HomeleftproxyY1"="HomeLeft$_{t-1}$",
              "HomeleftproxyY2"="HomeLeft$_{t-2}$",
              "HomeleftproxyY3"="HomeLeft$_{t-3}$",
              "HomeleftproxyY4"="HomeLeft$_{t-4}$",
              "Homencas_leftY1"="HomeLeft$_{t-1}$",
              "Homencas_leftY2"="HomeLeft$_{t-2}$",
              "Homencas_leftY3"="HomeLeft$_{t-3}$",
              "Homencas_leftY4"="HomeLeft$_{t-4}$",
              "HomerightproxyY1"="HomeRight$_{t-1}$",
              "HomerightproxyY2"="HomeRight$_{t-2}$",
              "HomerightproxyY3"="HomeRight$_{t-3}$",
              "HomerightproxyY4"="HomeRight$_{t-4}$",
              "Homencas_rightY1"="HomeRight$_{t-1}$",
              "Homencas_rightY2"="HomeRight$_{t-2}$",
              "Homencas_rightY3"="HomeRight$_{t-3}$",
              "Homencas_rightY4"="HomeRight$_{t-4}$",
              "HomeIslClaimY1"="HomeIslam$_{t-1}$",
              "HomeIslClaimY2"="HomeIslam$_{t-2}$",
              "HomeIslClaimY3"="HomeIslam$_{t-3}$",
              "HomeIslClaimY4"="HomeIslam$_{t-4}$",
              "Homencas_IslClY1"="HomeIslam$_{t-1}$",
              "Homencas_IslClY2"="HomeIslam$_{t-2}$",
              "Homencas_IslClY3"="HomeIslam$_{t-3}$",
              "Homencas_IslClY4"="HomeIslam$_{t-4}$",
              "geoleftproxyY1"="ForeignLeft$_{t-1}$",
              "geoleftproxyY2"="ForeignLeft$_{t-2}$",
              "geoleftproxyY3"="ForeignLeft$_{t-3}$",
              "geoleftproxyY4"="ForeignLeft$_{t-4}$",
              "geoncas_leftY1"="ForeignLeft$_{t-1}$",
              "geoncas_leftY2"="ForeignLeft$_{t-2}$",
              "geoncas_leftY3"="ForeignLeft$_{t-3}$",
              "geoncas_leftY4"="ForeignLeft$_{t-4}$",
              "georightproxyY1"="ForeignRight$_{t-1}$",
              "georightproxyY2"="ForeignRight$_{t-2}$",
              "georightproxyY3"="ForeignRight$_{t-3}$",
              "georightproxyY4"="ForeignRight$_{t-4}$",
              "geoncas_rightY1"="ForeignRight$_{t-1}$",
              "geoncas_rightY2"="ForeignRight$_{t-2}$",
              "geoncas_rightY3"="ForeignRight$_{t-3}$",
              "geoncas_rightY4"="ForeignRight$_{t-4}$",
              "geoIslClaimY1"="ForeignIslam$_{t-1}$",
              "geoIslClaimY2"="ForeignIslam$_{t-2}$",
              "geoIslClaimY3"="ForeignIslam$_{t-3}$",
              "geoIslClaimY4"="ForeignIslam$_{t-4}$",
              "geoncas_IslClY1"="ForeignIslam$_{t-1}$",
              "geoncas_IslClY2"="ForeignIslam$_{t-2}$",
              "geoncas_IslClY3"="ForeignIslam$_{t-3}$",
              "geoncas_IslClY4"="ForeignIslam$_{t-4}$",
              "culleftproxyY1"="ForeignLeft$_{t-1}$",
              "culleftproxyY2"="ForeignLeft$_{t-2}$",
              "culleftproxyY3"="ForeignLeft$_{t-3}$",
              "culleftproxyY4"="ForeignLeft$_{t-4}$",
              "culncas_leftY1"="ForeignLeft$_{t-1}$",
              "culncas_leftY2"="ForeignLeft$_{t-2}$",
              "culncas_leftY3"="ForeignLeft$_{t-3}$",
              "culncas_leftY4"="ForeignLeft$_{t-4}$",
              "culrightproxyY1"="ForeignRight$_{t-1}$",
              "culrightproxyY2"="ForeignRight$_{t-2}$",
              "culrightproxyY3"="ForeignRight$_{t-3}$",
              "culrightproxyY4"="ForeignRight$_{t-4}$",
              "culncas_rightY1"="ForeignRight$_{t-1}$",
              "culncas_rightY2"="ForeignRight$_{t-2}$",
              "culncas_rightY3"="ForeignRight$_{t-3}$",
              "culncas_rightY4"="ForeignRight$_{t-4}$",
              "culIslClaimY1"="ForeignIslam$_{t-1}$",
              "culIslClaimY2"="ForeignIslam$_{t-2}$",
              "culIslClaimY3"="ForeignIslam$_{t-3}$",
              "culIslClaimY4"="ForeignIslam$_{t-4}$",
              "culncas_IslClY1"="ForeignIslam$_{t-1}$",
              "culncas_IslClY2"="ForeignIslam$_{t-2}$",
              "culncas_IslClY3"="ForeignIslam$_{t-3}$",
              "culncas_IslClY4"="ForeignIslam$_{t-4}$")
#table_twfe_total <- 
msummary(twfe_right, gof_map = gof_map,
         coef_map=coef_map, stars=T) |>
  add_header_above(c(" "=1,"$c=0$"=1,"$c=c^{geo}$"=1,"$c=c^{cul}$"=1,
                     "$c=0$"=1,"$c=c^{geo}$"=1,"$c=c^{cul}$"=1)) |>
  add_header_above(c(" "=1,"Total Attacks"=3,"Total Casualties"=3))

save(twfe_right, gof_map, coef_map, glance_custom.fixest, 
     gen_stars, file="twfe_right.RData")

#right simple
data <- data |>
  mutate(
    Homeleftproxy = HomeleftproxyY1+HomeleftproxyY2+HomeleftproxyY3+HomeleftproxyY4,
    Homencas_left = Homencas_leftY1+Homencas_leftY2+Homencas_leftY3+Homencas_leftY4,
    Homerightproxy = HomerightproxyY1+HomerightproxyY2+HomerightproxyY3+HomerightproxyY4,
    Homencas_right = Homencas_rightY1+Homencas_rightY2+Homencas_rightY3+Homencas_rightY4,
    HomeIslClaim = HomeIslClaimY1+HomeIslClaimY2+HomeIslClaimY3+HomeIslClaimY4,
    Homencas_IslCl = Homencas_IslClY1+Homencas_IslClY2+Homencas_IslClY3+Homencas_IslClY4,
    geoleftproxy = geoleftproxyY1+geoleftproxyY2+geoleftproxyY3+geoleftproxyY4,
    geoncas_left = geoncas_leftY1+geoncas_leftY2+geoncas_leftY3+geoncas_leftY4,
    georightproxy = georightproxyY1+georightproxyY2+georightproxyY3+georightproxyY4,
    geoncas_right = geoncas_rightY1+geoncas_rightY2+geoncas_rightY3+geoncas_rightY4,
    geoIslClaim = geoIslClaimY1+geoIslClaimY2+geoIslClaimY3+geoIslClaimY4,
    geoncas_IslCl = geoncas_IslClY1+geoncas_IslClY2+geoncas_IslClY3+geoncas_IslClY4,
    culleftproxy = culleftproxyY1+culleftproxyY2+culleftproxyY3+culleftproxyY4,
    culncas_left = culncas_leftY1+culncas_leftY2+culncas_leftY3+culncas_leftY4,
    culrightproxy = culrightproxyY1+culrightproxyY2+culrightproxyY3+culrightproxyY4,
    culncas_right = culncas_rightY1+culncas_rightY2+culncas_rightY3+culncas_rightY4,
    culIslClaim = culIslClaimY1+culIslClaimY2+culIslClaimY3+culIslClaimY4,
    culncas_IslCl = culncas_IslClY1+culncas_IslClY2+culncas_IslClY3+culncas_IslClY4)
outcome <- "PopulistRight~"
homeattacks <- "Homeleftproxy+Homerightproxy+HomeIslClaim"
fe <- "|Country+Year"
reg_right_attack_none <- feols(as.formula(paste0(outcome,homeattacks,fe)),data)
controls <- "+GrowthY1+GrowthY2+GrowthY3+UnemploymentY1+UnemploymentY2+UnemploymentY3+
Export_penetrationY1+Export_penetrationY2+Export_penetrationY3+
Import_penetrationY1+Import_penetrationY2+Import_penetrationY3+
Immigration_percentY1+Immigration_percentY2+Immigration_percentY3+
Emigration_percentY1+Emigration_percentY2+Emigration_percentY3"
reg_right_attack_control <- feols(as.formula(paste0(outcome,homeattacks,controls,fe)),data)
geoattacks <- "+geoleftproxy+georightproxy+geoIslClaim"
reg_right_attack_geo <- feols(as.formula(paste0(outcome,homeattacks,geoattacks,controls,fe)),data)
culattacks <- "+culleftproxy+culrightproxy+culIslClaim"
reg_right_attack_cul <- feols(as.formula(paste0(outcome,homeattacks,culattacks,controls,fe)),data)

homecasualty <- "Homencas_left+Homencas_right+Homencas_IslCl"
geocasualty <- "+geoncas_left+geoncas_right+geoncas_IslCl"
culcasualty <- "+culncas_left+culncas_right+culncas_IslCl"
reg_right_casualty_none <- feols(as.formula(paste0(outcome,homecasualty,fe)),data)
reg_right_casualty_control <- feols(as.formula(paste0(outcome,homecasualty,controls,fe)),data)
reg_right_casualty_geo <- feols(as.formula(paste0(outcome,homecasualty,geocasualty,controls,fe)),data)
reg_right_casualty_cul <- feols(as.formula(paste0(outcome,homecasualty,culcasualty,controls,fe)),data)
twfe_right <- list(reg_right_attack_control, 
                   reg_right_attack_geo, reg_right_attack_cul,
                   reg_right_casualty_control,
                   reg_right_casualty_geo, reg_right_casualty_cul)
for (i in 1:length(twfe_right)){
  twfe_right[[i]]$df.residual <- parameters::degrees_of_freedom(twfe_right[[i]])
  home <- names(coef(twfe_right[[i]]))[grepl("Home",names(coef(twfe_right[[i]])))]
  twfe_right[[i]]$homeF <- lht(twfe_right[[i]],test='F',
                               c(paste0(home,"=0")))
  foreign <- character()
  foreign <- names(coef(twfe_right[[i]]))[grepl("geo|cul",names(coef(twfe_right[[i]])))]
  if (length(foreign)>0){
    twfe_right[[i]]$foreignF <- lht(twfe_right[[i]],test='F',
                                    c(paste0(foreign,"=0"))
    )
  }
}
coef_map <- c("Homeleftproxy"="HomeLeft",
              "Homencas_left"="HomeLeft",
              "Homerightproxy"="HomeRight",
              "Homencas_right"="HomeRight",
              "HomeIslClaim"="HomeIslam",
              "Homencas_IslCl"="HomeIslam",
              "geoleftproxy"="ForeignLeft",
              "geoncas_left"="ForeignLeft",
              "georightproxy"="ForeignRight",
              "geoncas_right"="ForeignRight",
              "geoIslClaim"="ForeignIslam",
              "geoncas_IslCl"="ForeignIslam",
              "culleftproxy"="ForeignLeft",
              "culncas_left"="ForeignLeft",
              "culrightproxy"="ForeignRight",
              "culncas_right"="ForeignRight",
              "culIslClaim"="ForeignIslam",
              "culncas_IslCl"="ForeignIslam")
msummary(twfe_right, gof_map = gof_map,
         coef_map=coef_map, stars=T) |>
  add_header_above(c(" "=1,"$c=0$"=1,"$c=c^{geo}$"=1,"$c=c^{cul}$"=1,
                     "$c=0$"=1,"$c=c^{geo}$"=1,"$c=c^{cul}$"=1)) |>
  add_header_above(c(" "=1,"Total Attacks"=3,"Total Casualties"=3))

save(twfe_right, gof_map, coef_map, glance_custom.fixest, 
     gen_stars, file="twfe_right.RData")
