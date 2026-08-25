state_policies <- read.csv("zzz/state_policies_2026-08-24.csv")
summary(state_policies)
head(state_policies)

## Fix error
fix_obs <- with(
    state_policies,
    which(
        variable == "w_gayrights_civilunions_marriage" &
        policy_short_description == "Conversion Therapy Ban"
    ))
state_policies$variable[fix_obs] <- "gayrights_conversion_therapy"

sort(with(state_policies, tapply(policy_real, year, \(x) sum(!is.na(x)))))

state_policies_2010_2012 <-
    subset(state_policies, year %in% 2010:2012 & state_abb != "DC") |>
    subset(variable != "w_gayrights_ban_conversiontherapy")

state_policies_2010_2012$policy_variable <- state_policies_2010_2012$variable
state_policies_2010_2012$value <- state_policies_2010_2012$policy
state_policies_2010_2012$value_real <- state_policies_2010_2012$policy_real

state_policies_2010_2012 <- subset(
    state_policies_2010_2012,
    select = c(year, state_abb, state_name, policy_variable,
               policy_short_description, policy_longer_description,
               value, value_real)
)

summary(state_policies_2010_2012)

state_policies_2010_2012[is.na(state_policies_2010_2012$value), ]

state_policies_2010_2012 |>
    subset(policy_variable == "w_gayrights_civilunions_marriage") |>
    (\(d) table(d$year, d$state_abb))()

saveRDS(state_policies_2010_2012, file = "data/state_policies_2010_2012.rda")


str(state_policies_2010_2012)
summary(state_policies_2010_2012)
table(state_policies_2010_2012$policy_variable)
