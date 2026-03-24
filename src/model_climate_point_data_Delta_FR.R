#Title: model_climate_point_data_Delta_FR.R
#Author: FR
#Date: 01.25.2026

library(dplyr)
library(nnet)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(gtsummary)
library(ggplot2)
library(forcats)
library(trend)

#home=setwd("C:\\Users\\frusere\\OneDrive - University of Idaho\\Desktop\\UOI\\DATA")

#load survey data
F<-read.csv("./data/survey/FW2L_DC.csv") 

#load climate data
daily_data_BV<-read.csv("./data/nclimgrid_point_data/daily_data_BV.csv")
daily_data<-read.csv("./data/nclimgrid_point_data/tmax_nclimgrid.csv")


# Convert all variables to factors except Q2 and Q8
F <- F %>%
  mutate(across(
    .cols = -c(ID2, ID8),      # exclude Q2 and Q8
    .fns = ~ as.factor(.)
  ))

F$ID2<-as.numeric(F$ID2)
F$ID8<-as.numeric(F$ID8)

# Create a new categorical variable based on ID2
F$ID2_Dcat <- cut(
  F$ID2,
  breaks = c(-Inf, 9, 20, 30, 40, Inf),  # define the intervals
  labels = c("0-9", "10-20", "21-30", "31-40", ">40"),  # assign labels
  right = TRUE  # intervals include the right endpoint
)

F <- F %>%
  mutate(
    ID31_3level = case_when(
      ID31 == "Advanced college/Graduate degre" ~ "Advanced/Post graduate level",
      ID31 %in% c("College undergraduate degree",
                  "Some College",
                  "Technical or vocational school") ~ "College graduate level",
      ID31 %in% c("Some high school",
                  "High school graduate") ~ "High school level",
      TRUE ~ NA_character_
    )
  )


F <- F %>%
  mutate(
    ID31_3level = factor(
      ID31_3level,
      levels = c("High school level", "College graduate level", "Advanced/Post graduate level"),
      ordered = TRUE
    )
  )

F <- F %>%
  mutate(
    ID32_3level = case_when(ID32 == "Unemployed" ~ "Unemployed",
                            ID32 == "Farming/Forestry" ~ "Farming/Forestry",
                            ID32 == "Retired" ~ "Retired",
                            TRUE ~ "Employed in various sectors"   # all remaining categories
    )
  )

F$ID32_3level<-as.factor(F$ID32_3level)
F$ID32_3level <- factor(F$ID32_3level,
                        levels = c("Unemployed",
                                   "Farming/Forestry",
                                   "Employed in various sectors",
                                   "Retired"))

levels(F$ID32_3level)

F <- F %>%
  mutate(IDMID8 = (ID8 - 32) * 5/9)

F6_Bovill<-F %>%
  filter(Nearest_Town == "Bovill")



# Bovill senslope from 1951-2025

bovill_data <- daily_data_BV 

colnames(bovill_data)
str(bovill_data)
bovill_data <- daily_data_BV

yearly_avg <- bovill_data %>%
  group_by(year) %>%
  summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")


sen <- sens.slope(
  yearly_avg$mean_temp,
  yearly_avg$year
)

sen_slope <- as.numeric(sen$estimates)  # °C per year
slope_label <- paste0("Sen slope = ", round(sen_slope, 3), " °C/year")

yearly_avg <- yearly_avg %>%
  mutate(
    sen_line = median(mean_temp) +
      sen_slope * (year - median(year))
  )


p <- ggplot(yearly_avg, aes(x = year, y = mean_temp)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkblue") +
  geom_line(aes(y = sen_line), color = "red", linewidth = 1) +
  annotate(
    "text",
    x = min(yearly_avg$year) + 2,
    y = max(yearly_avg$mean_temp) - 0.5,
    label = slope_label,
    hjust = 0,
    size = 5,
    color = "red"
  ) +
  labs(
    title = "Trend in Annual Average Temperature in Bovill",
    x = "Year",
    y = "Mean Annual Temperature (°C)"
  ) +
  theme_classic()

p

########################################################################################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    Slope = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      yearly_avg <- bovill_data %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")
      
      if (nrow(yearly_avg) >= 2) {
        as.numeric(
          sens.slope(
            yearly_avg$mean_temp,
            yearly_avg$year
          )$estimates
        )
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()


t1<-ggplot(F6_Bovill, aes(x = ID6, y = Slope)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    title = "Relationship between perceived change and temperature Trend (Slope)",
    x = "que_6",
    y = "Slope of Temperature Trend (°C/year)"
  ) +
  theme_classic()

t1

#boxplot

t1_box <- ggplot(F6_Bovill, aes(x = factor(ID6), y = Slope)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(
    title = "Temperature Trend (Slope) by perceived change (que_6)",
    x = "que_6 (ID6)",
    y = "Slope of Temperature Trend (°C/year)"
  ) +
  theme_classic()

t1_box

# Number of days with Temp => 25, =>30,=> 35 Bovill

yearly_hot_days <- daily_data_BV %>%
  group_by(year) %>%
  summarize(
    days_25 = sum(value >= 25, na.rm = TRUE),
    days_30 = sum(value >= 30, na.rm = TRUE),
    days_35 = sum(value >= 35, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate 90th and 95th percentile for each calendar day across all years
daily_percentiles <- daily_data_BV %>%
  group_by(month, day) %>%
  summarize(
    p90 = quantile(value, 0.9, na.rm = TRUE),
    p95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

daily_data_BV <- daily_data_BV %>%
  left_join(daily_percentiles, by = c("month", "day"))

daily_data_BV <- daily_data_BV %>%
  mutate(
    exceed_90 = value > p90,
    exceed_95 = value > p95
  )

yearly_extremes <- daily_data_BV %>%
  group_by(year) %>%
  summarize(
    days_above_90 = sum(exceed_90, na.rm = TRUE),
    days_above_95 = sum(exceed_95, na.rm = TRUE)
  )

sen_90 <- sens.slope(yearly_extremes$days_above_90, yearly_extremes$year)$estimates
sen_95 <- sens.slope(yearly_extremes$days_above_95, yearly_extremes$year)$estimates

#######################################################################################################


# Step 1: Compute Sen's slope for each threshold
thresholds <- c("days_25", "days_30", "days_35")

# Function to get Sen slope line
get_sen_line <- function(years, values) {
  if(length(values) >= 2) {
    sen <- sens.slope(values, years)
    slope <- as.numeric(sen$estimates)
    line <- median(values) + slope * (years - median(years))
    return(line)
  } else {
    return(rep(NA, length(values)))
  }
}

# Step 2: Add Sen slope lines to the dataset
yearly_hot_days <- yearly_hot_days %>%
  mutate(
    line_25 = get_sen_line(year, days_25),
    line_30 = get_sen_line(year, days_30),
    line_35 = get_sen_line(year, days_35)
  )

# Step 3: Plot
p <- ggplot(yearly_hot_days, aes(x = year)) +
  geom_line(aes(y = days_25, color = "≥25°C"), linewidth = 1) +
  geom_line(aes(y = line_25, color = "≥25°C"), linetype = "dashed") +
  geom_line(aes(y = days_30, color = "≥30°C"), linewidth = 1) +
  geom_line(aes(y = line_30, color = "≥30°C"), linetype = "dashed") +
  geom_line(aes(y = days_35, color = "≥35°C"), linewidth = 1) +
  geom_line(aes(y = line_35, color = "≥35°C"), linetype = "dashed") +
  scale_color_manual(
    name = "Threshold",
    values = c("≥25°C" = "blue", "≥30°C" = "orange", "≥35°C" = "red")
  ) +
  labs(
    title = "Number of Hot Days per Year in Bovill",
    x = "Year",
    y = "Number of Hot Days"
  ) +
  theme_classic() +
  theme(legend.position = "top")



p

##############################################################################################################################
library(dplyr)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_25_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 25, na.rm = TRUE)) %>%
        pull()
    },
    
    days_30_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 30, na.rm = TRUE)) %>%
        pull()
    },
    
    days_35_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 35, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()


View(F6_Bovill)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_above_90_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(exceed_90, na.rm = TRUE)) %>%
        pull()
    },
    
    days_above_95_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(exceed_95, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()


View(F6_Bovill)
#############################################################################################################33

yearly_extremes <- daily_data_BV %>%
  group_by(year) %>%
  summarize(
    days_25 = sum(value >= 25, na.rm = TRUE),
    days_30 = sum(value >= 30, na.rm = TRUE),
    days_35 = sum(value >= 35, na.rm = TRUE),
    days_above_90 = sum(exceed_90, na.rm = TRUE),
    days_above_95 = sum(exceed_95, na.rm = TRUE),
    .groups = "drop"
  )


library(trend)
library(dplyr)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_days_30 = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      ts_data <- yearly_extremes %>%
        filter(year >= start_year, year <= end_year) %>%
        select(year, days_30)
      
      if (nrow(ts_data) >= 2) {
        as.numeric(
          sens.slope(ts_data$days_30, ts_data$year)$estimates
        )
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()
###########################################################################################################################################

View(F6_Bovill)

########################################################################################################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    
    slope_days_25 = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      d <- yearly_extremes %>%
        filter(year >= start_year, year <= end_year)
      
      if (nrow(d) >= 2)
        as.numeric(sens.slope(d$days_25, d$year)$estimates)
      else NA_real_
    },
    
    slope_days_30 = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      d <- yearly_extremes %>%
        filter(year >= start_year, year <= end_year)
      
      if (nrow(d) >= 2)
        as.numeric(sens.slope(d$days_30, d$year)$estimates)
      else NA_real_
    },
    
    slope_days_35 = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      d <- yearly_extremes %>%
        filter(year >= start_year, year <= end_year)
      
      if (nrow(d) >= 2)
        as.numeric(sens.slope(d$days_35, d$year)$estimates)
      else NA_real_
    },
    
    slope_days_95 = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      d <- yearly_extremes %>%
        filter(year >= start_year, year <= end_year)
      
      if (nrow(d) >= 2)
        as.numeric(sens.slope(d$days_above_95, d$year)$estimates)
      else NA_real_
    }
    
  ) %>%
  ungroup()

View(F6_Bovill)

#################################################################################################


summary(F6_Bovill)

summary(F6_Bovill)

##########################################################################################################################################################
library(dplyr)

# Compute participant-specific number of hot days
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    hot_days = {
      # Participant exposure period
      start_year <- max(1951, 2025 - ID2)  # duration of residence
      end_year   <- 2025
      
      # Use participant-specific hot day threshold from IDMID8
      threshold <- IDMID8 
      
      # Filter daily data to participant's exposure period
      participant_days <- bovill_data %>%
        filter(year >= start_year, year <= end_year)
      
      # Count number of days exceeding threshold
      sum(participant_days$value >= threshold, na.rm = TRUE)
    }
  ) %>%
  ungroup()


library(dplyr)
library(trend)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_hot_days = {
      
      # Participant exposure period
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      # Participant-specific threshold
      threshold <- IDMID8
      
      # Build annual time series of hot days
      yearly_hot <- bovill_data %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(
          hot_days = sum(value >= threshold, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Compute Sen's slope if enough data
      if (nrow(yearly_hot) >= 2) {
        as.numeric(
          sens.slope(
            yearly_hot$hot_days,
            yearly_hot$year
          )$estimates
        )
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()


t1<-ggplot(F6_Bovill, aes(x = ID6, y = slope_hot_days)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    title = "Relationship between perceived change and temperature Trend (Slope)",
    x = "que_6",
    y = "Slope of Temperature Trend (°C/year)"
  ) +
  theme_classic()
t1

#boxplot

t1_box <- ggplot(F6_Bovill, aes(x = factor(ID6), y = slope_hot_days)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, color = "darkgreen") +
  labs(
    title = "Hot-days trend (slope) by perceived change (que_6)",
    x = "que_6 (ID6)",
    y = "Slope of hot days trend (units/year)"
  ) +
  theme_classic()

t1_box

############################################################################################################################################################################################



season_filter <- function(data, start_month, end_month) {
  data %>% filter(month >= start_month, month <= end_month)
}


bovill_AO <- bovill_data %>%
  season_filter(4, 10)

yearly_avg_AO <- bovill_AO %>%
  group_by(year) %>%
  summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")

sen_AO <- sens.slope(yearly_avg_AO$mean_temp, yearly_avg_AO$year)
sen_AO$estimates

# May–September
bovill_MS <- bovill_data %>% season_filter(5, 9)

# June–August
bovill_JJA <- bovill_data %>% season_filter(6, 8)
###########################################################################################################################
library(trend)
library(dplyr)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_temp_AO = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      yearly_avg <- bovill_data %>%
        season_filter(4, 10) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")
      
      if (nrow(yearly_avg) >= 2)
        as.numeric(sens.slope(yearly_avg$mean_temp, yearly_avg$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()
###########################################################################################################################################################
library(trend)
library(dplyr)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_temp_JJA = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      yearly_avg <- bovill_data %>%
        season_filter(6, 8) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")
      
      if (nrow(yearly_avg) >= 2)
        as.numeric(sens.slope(yearly_avg$mean_temp, yearly_avg$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()
####################################################################################################################
library(trend)
library(dplyr)

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_temp_MS = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      yearly_avg <- bovill_data %>%
        season_filter(5, 9) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(mean_temp = mean(value, na.rm = TRUE), .groups = "drop")
      
      if (nrow(yearly_avg) >= 2)
        as.numeric(sens.slope(yearly_avg$mean_temp, yearly_avg$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()
###########################################################################################################
#Seasonal hot-day counts (absolute thresholds)
#APRIL–OCTOBER example

yearly_hot_AO <- daily_data_BV %>%
  season_filter(4, 10) %>%
  group_by(year) %>%
  summarize(
    days_25 = sum(value >= 25, na.rm = TRUE),
    days_30 = sum(value >= 30, na.rm = TRUE),
    days_35 = sum(value >= 35, na.rm = TRUE),
    .groups = "drop"
  )

#########################################################################################################################
#Seasonal hot-day counts (absolute thresholds)
#MAY–SEPT example
yearly_hot_MS <- daily_data_BV %>%
  season_filter(5, 9) %>%
  group_by(year) %>%
  summarize(
    days_25 = sum(value >= 25, na.rm = TRUE),
    days_30 = sum(value >= 30, na.rm = TRUE),
    days_35 = sum(value >= 35, na.rm = TRUE),
    .groups = "drop"
  )
#####################################################################################################################################
#Seasonal hot-day counts (absolute thresholds)
#JUN–AUG example



yearly_hot_JJA <- daily_data_BV %>%
  season_filter(6, 8) %>%
  group_by(year) %>%
  summarize(
    days_25 = sum(value >= 25, na.rm = TRUE),
    days_30 = sum(value >= 30, na.rm = TRUE),
    days_35 = sum(value >= 35, na.rm = TRUE),
    .groups = "drop"
  )


#####################################################################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_30_AO_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(4, 10) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 30, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

#########################################################################################################################################

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_30_JJA_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(6, 8) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 30, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

m <- lm(days_30_JJA_exp ~ ID2, data = F6_Bovill)
F6_Bovill$resid_days30JJA <- resid(m)

ggplot(F6_Bovill, aes(x = factor(ID6), y = resid_days30JJA)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(
    title = "Residual JJA ≥30°C exposure (after removing years lived) by perception (ID6)",
    x = "Perceived change (ID6)",
    y = "Residual exposure (days)"
  ) +
  theme_classic()



#################################################################################################################################

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_30_MS_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(5, 9) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 30, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()


############################################################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_35_AO_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(4, 10) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 35, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

##########################################################################################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_35_JJA_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(6, 8) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 35, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

#############################################################################
F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_35_MS_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(5, 9) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 35, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

###############################################################################################3

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_25_AO_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(4, 10) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 25, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()



F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_25_JJA_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(6, 8) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 25, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()
###########################################################################################################



F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    days_25_MS_exp = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      
      daily_data_BV %>%
        season_filter(5, 9) %>%
        filter(year >= start_year, year <= end_year) %>%
        summarize(sum(value >= 25, na.rm = TRUE)) %>%
        pull()
    }
  ) %>%
  ungroup()

######################################################################################################################
# Participant-specific seasonal Sen’s slope of hot days APRIL-OCTOBER


F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_hot_AO = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      threshold  <- IDMID8
      
      yearly_hot <- bovill_data %>%
        season_filter(4, 10) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(
          hot_days = sum(value >= threshold, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(yearly_hot) >= 2)
        as.numeric(sens.slope(yearly_hot$hot_days, yearly_hot$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()


###########################################################################################################################
#Participant-specific seasonal Sen’s slope of hot days MAY-SEPTEMBER

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_hot_MS = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      threshold  <- IDMID8
      
      yearly_hot <- bovill_data %>%
        season_filter(5, 9) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(
          hot_days = sum(value >= threshold, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(yearly_hot) >= 2)
        as.numeric(sens.slope(yearly_hot$hot_days, yearly_hot$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()

########################################################################################################################

#Participant-specific seasonal Sen’s slope of hot days  JUN-AUG

F6_Bovill <- F6_Bovill %>%
  rowwise() %>%
  mutate(
    slope_hot_JJA = {
      start_year <- max(1951, 2025 - ID2)
      end_year   <- 2025
      threshold  <- IDMID8
      
      yearly_hot <- bovill_data %>%
        season_filter(6, 8) %>%
        filter(year >= start_year, year <= end_year) %>%
        group_by(year) %>%
        summarize(
          hot_days = sum(value >= threshold, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(yearly_hot) >= 2)
        as.numeric(sens.slope(yearly_hot$hot_days, yearly_hot$year)$estimates)
      else NA_real_
    }
  ) %>%
  ungroup()
########################################################################################################################################################
