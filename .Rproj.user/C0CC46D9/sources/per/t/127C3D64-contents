library(tidyverse)
library(deSolve)

setup_initial_conditions <- function(N_total, num_dose_1, num_dose_2, VE1, VE2, initial_I, initial_E) {
  n_dose2_only <- num_dose_2                  
  n_dose1_only <- num_dose_1 - num_dose_2     
  n_unvax <- N_total - num_dose_1             
  VE1 <- as.numeric(VE1); VE2 <- as.numeric(VE2)
  
  S_init <- n_unvax - initial_E - initial_I
  V1_init <- n_dose1_only * VE1
  S1_init <- n_dose1_only * (1 - VE1)
  V2_init <- n_dose2_only * VE2
  S2_init <- n_dose2_only * (1 - VE2)
  
  init_state <- c(
    S = S_init, E = initial_E, I = initial_I, R = 0,
    V1 = V1_init, S1 = S1_init, E1 = 0, I1 = 0, R1 = 0,
    V2 = V2_init, S2 = S2_init, E2 = 0, I2 = 0, R2 = 0,
    V3 = 0, S3 = 0, E3 = 0, I3 = 0, R3 = 0,
    Cum_Incidence = initial_I 
  )
  return(init_state)
}

seir_epidemic_model <- function(time, state, parameters, switch_cl, switch_kcl) {
  with(as.list(c(state, parameters)), {
    
    N <- S+E+I+R + V1+S1+E1+I1+R1 + V2+S2+E2+I2+R2 + V3+S3+E3+I3+R3
    lambda <- beta * (I + I1 + I2 + I3) / N 

    is_campaign_active <- (time >= day_camp_start & time <= day_camp_end)
    
    rate_cl <- ifelse(switch_cl == TRUE & is_campaign_active, v_cl, 0)
    rate_kcl <- ifelse(switch_kcl == TRUE & is_campaign_active, v_kcl, 0)

    dS <- -lambda*S - v_tcmr*S - rate_cl*S - rate_kcl*S
    dE <- lambda*S - sigma*E
    dI <- sigma*E - gamma*I
    dR <- gamma*I

    dV1 <- (v_tcmr + rate_cl + rate_kcl)*S*VE1 - rate_cl*V1 - rate_kcl*V1
    dS1 <- (v_tcmr + rate_cl + rate_kcl)*S*(1-VE1) - lambda*S1 - rate_cl*S1 - rate_kcl*S1
    dE1 <- lambda*S1 - sigma*E1
    dI1 <- sigma*E1 - gamma*I1
    dR1 <- gamma*I1

    dV2 <- (v_tcmr + rate_cl + rate_kcl)*S1*VE2 + (rate_cl + rate_kcl)*V1*VE2 - rate_kcl*V2
    dS2 <- (v_tcmr + rate_cl + rate_kcl)*S1*(1-VE2) + (rate_cl + rate_kcl)*V1*(1-VE2) - lambda*S2 - rate_kcl*S2
    dE2 <- lambda*S2 - sigma*E2
    dI2 <- sigma*E2 - gamma*I2
    dR2 <- gamma*I2

    dV3 <- rate_kcl*V2*VE3 + rate_kcl*S2*VE3
    dS3 <- rate_kcl*V2*(1-VE3) + rate_kcl*S2*(1-VE3) - lambda*S3
    dE3 <- lambda*S3 - sigma*E3
    dI3 <- sigma*E3 - gamma*I3
    dR3 <- gamma*I3

    dCum_Incidence <- sigma*E + sigma*E1 + sigma*E2 + sigma*E3
    
    return(list(c(dS, dE, dI, dR, 
                  dV1, dS1, dE1, dI1, dR1, 
                  dV2, dS2, dE2, dI2, dR2, 
                  dV3, dS3, dE3, dI3, dR3, 
                  dCum_Incidence)))
  })
}

df <- readRDS("cases_clean.rds")
df1 <- df |> filter(episode == 1) |> select(dates, birth, age_month, tinh_trang, episode)

df1_final <- df1 |>
  mutate(dates = as.Date(dates)) |>
  group_by(dates) |>
  summarise(cases = n(), .groups = "drop") |>
  complete(dates = seq.Date(min(dates), max(dates), by = "day"), fill = list(cases = 0)) |>
  mutate(
    cum_cases = cumsum(cases),
    time = row_number() - 1 
  )

thoi_gian_dich_2018 <- df1_final$time

parameters_2018 <- c(
  beta = 1.0,   
  sigma = 1/14, 
  gamma = 1/10, 
  
  VE1 = 0.95, 
  VE2 = 0.98, 
  VE3 = 0.98,
  
  v_tcmr = 0.001,      
  v_cl = 0.01,         
  v_kcl = 0.015,    

  day_camp_start = 100, 
  day_camp_end = 160    
)

init_2018 <- setup_initial_conditions(
  N_total = 1200000, 
  num_dose_1 = 1080000, 
  num_dose_2 = 960000,  
  VE1 = parameters_2018["VE1"], 
  VE2 = parameters_2018["VE2"], 
  initial_I = df1_final$cases[1], 
  initial_E = df1_final$cases[1] * 3 
)

cost_function_total <- function(beta_guess) {
  temp_params <- parameters_2018
  temp_params["beta"] <- beta_guess

  out_fit <- ode(y = init_2018, times = thoi_gian_dich_2018, func = seir_epidemic_model, 
                 parms = temp_params, switch_cl = FALSE, switch_kcl = TRUE)
  
  total_mo_hinh <- tail(as.data.frame(out_fit)$Cum_Incidence, 1)
  total_thuc_te <- tail(df1_final$cum_cases, 1) 
  
  return((total_thuc_te - total_mo_hinh)^2)
}

fit_result <- optim(par = 0.5, fn = cost_function_total, method = "Brent", lower = 0.01, upper = 2.0)
beta_chuan_2018 <- fit_result$par

parameters_2018["beta"] <- beta_chuan_2018

run_sim_cases <- function(sw_cl, sw_kcl) {
  out <- ode(y = init_2018, times = thoi_gian_dich_2018, func = seir_epidemic_model, 
             parms = parameters_2018, switch_cl = sw_cl, switch_kcl = sw_kcl)
  return(tail(as.data.frame(out)$Cum_Incidence, 1))
}

ca_thuc_te_kcl <- run_sim_cases(sw_cl = FALSE, sw_kcl = TRUE)
ca_gia_dinh_base <- run_sim_cases(sw_cl = FALSE, sw_kcl = FALSE)
ca_gia_dinh_cl <- run_sim_cases(sw_cl = TRUE, sw_kcl = FALSE)
cat(sprintf("=> Số ca ngăn ngừa được:   %.0f ca\n", ca_gia_dinh_base - ca_thuc_te_kcl))
