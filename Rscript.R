#   SARS-CoV2 Research Project in Ghana using data from GISAID, and owid
#   In collaboration with Nidia (USA), Nicholas and Ann (Ghana) and Elijah et al (ghana)
#   script written by Olayinka Sunday OKOH 
#   22nd April, 2021
#   The data were downloaded on 19th April, 2021 from gisaid and owid

library(readr) # To read in the tsv file
library(dplyr) # For data manipulation and wrangling
library(ggplot2) # For Visualization
library(readxl) # For reading in excel files
library(tmap) # For Maps
library(tmaptools) #  For Maps
library(sf) # For Maps


#   Read in  the file as raw_data
raw_data<- read_tsv("raw_data/Ghana_gisaid_hcov-19_2021_04_19_15.tsv", col_names = TRUE,
                    col_types = cols(strain ="f", virus = "c", gisaid_epi_isl = "-",
                                     genbank_accession = "-", date = "D", region = "f",
                                     country = "c",  division = "c", location = "c",
                                     region_exposure = "c", country_exposure = "c", division_exposure = "c",
                                     segment = "c", host = "f", age = "i", sex = "f",
                                     Nextstrain_clade = "f", pangolin_lineage = "f", GISAID_clade = "f",
                                     originating_lab = "c", submitting_lab = "c", authors = "-", url = "-", title = "-",
                                     paper_url = "-",date_submitted = "D" 
                    ))


#Change file to ghana_gisaid and select some columns to use for the analysis
ghana_gisaid <- raw_data %>% select(-c(`Virus name`, `Accession ID`, Comment, `Comment type`, `Additional host information`, Passage, Host)) %>% as.data.frame()


# Exploratory analysis
dim(ghana_gisaid) # 360 sequence submitted

# Get the number of sequences submitted from each state
Ghana_sequences <- ghana_gisaid %>% select(Location) %>% group_by(Location) %>% summarize (sequences = n()) %>% mutate(Percentage = round((sequences/sum(sequences)*100),1))

#   Export the tables as Sequences submitted to gisaid from each ghana region 
write.csv(Ghana_sequences, "outputs/Submitted_Sequences.csv")

#   Of the 360 sequences submitted from ghana, 261 could not be traced to any region 


#lineages circulating in ghana are 36
length(unique(ghana_gisaid$Lineage)) # They are 36

# Analyze the the SARS-CoV2 lineages in ghana
Lineages <- ghana_gisaid %>% group_by(Lineage) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1)) %>% arrange(desc(N))
write.csv(Lineages, "outputs/Lineages.csv")

# Analyze the GISAID clades in ghana
GISAID_clade <- ghana_gisaid %>% group_by(Clade) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1)) %>% arrange(desc(N))
write.csv(GISAID_clade, "outputs/Clades.csv")

# Analyze the Sequencing Technology used in ghana
Sequencing_Technology <- ghana_gisaid %>% group_by(`Sequencing technology`) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1)) %>% arrange(desc(N))
write.csv(Sequencing_Technology, "outputs/Sequencing_Tech.csv")

# Analyze the Specimens used in Ghana
Specimen <- ghana_gisaid %>% group_by(Specimen) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1)) %>% arrange(desc(N))
write.csv(Specimen, "outputs/Specimen.csv")

# Assebly method used in Ghana
Assembly_method <- ghana_gisaid %>% group_by(`Assembly method`) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1)) %>% arrange(desc(N))
write.csv(Assembly_method, "outputs/Assembly_method.csv")





#read in data from OWID - Our World in Data
owid_data <- read_csv("raw_data/owid-covid-data.csv")

ghana_owid <- as.data.frame(ghana_owid)

ghana_owid <- owid_data %>% filter(location == "Ghana") %>% select(c(location, date, total_cases, total_cases_per_million,
                                                                     total_deaths, total_deaths_per_million, total_tests, total_tests_per_thousand, positive_rate,
                                                             tests_per_case, total_vaccinations_per_hundred)) %>% as.data.frame()

# Effect of Lockdown on COVID-19 incidence
ggplot(ghana_owid, aes(x = date, y = total_cases)) + 
  geom_line(color = 'red') + 
  geom_vline(xintercept = as.Date(c("2020-03-30", "2020-07-01", "2020-10-01"))) +
  labs(title = "Incidence of COVID-19 in Ghana", y = "Number of Cases", x = "Date")
  ggsave("outputs/ghana_covid19_incidence.jpg", height =  120, width = 200, units = "mm")
  ggsave("outputs/ghana_covid19_incidence.pdf", height =  120, width = 200, units = "mm")

  
  # Deaths from COVID-19
  ggplot(ghana_owid, aes(x = date, y = total_deaths)) + 
    geom_line(color = 'red') + 
    geom_vline(xintercept = as.Date(c("2020-03-30", "2020-07-01", "2020-10-01"))) +
    labs(title = "Deaths from COVID-19 in Ghana", y = "Number of deaths", x = "Date")
  ggsave("outputs/deaths_from_covid19.jpg", height =  120, width = 200, units = "mm")
  ggsave("outputs/deaths_from_covid19.pdf", height =  120, width = 200, units = "mm")
 
  
  # COVID-19 Tests in Ghana
  ggplot(ghana_owid, aes(x = date, y = total_tests)) + 
    geom_line(color = 'red') + 
    geom_vline(xintercept = as.Date(c("2020-03-30", "2020-07-01", "2020-10-01"))) +
    labs(title = "COVID-19 tests in Ghana", y = "Number of tests", x = "Date")
  ggsave("outputs/covid19_tests.jpg", height =  120, width = 200, units = "mm")
  ggsave("outputs/covid19_tests.pdf", height =  120, width = 200, units = "mm")
  
# Positivity rate
  ggplot(ghana_owid, aes(x = date, y = positive_rate)) + 
    geom_line(color = 'red') + 
    #geom_vline(xintercept = as.Date(c("2020-03-30", "2020-07-01", "2020-10-01"))) +
    labs(title = "COVID-19 Positivity rate in Ghana", y = "", x = "Date")
  ggsave("outputs/positive_rate.jpg", height =  120, width = 200, units = "mm")
  ggsave("outputs/positive_rate.pdf", height =  120, width = 200, units = "mm")
  
  
  

#    ====================== MAPS ================================
#   No data for this yet
# Read in the ghana shape file which was downloaded from gadm.org
ghana_sh <- st_read("raw_data/gadm36_GHA_shp/gadm36_GHA_1.shp")

#Let COVID-19 data in ghana be ghana_data and rename State as NAME_1 since that is what is used in Tmap
#ghana_data <- ghana_cases_and_sequences %>% rename(NAME_1 = State)

#Join the map data and covid-19 data and name it as ghana
#ghana <- left_join(ghana_sh, ghana_data, by = "NAME_1")

##Map of COVID-19 Cases as reported by NCDC
ghana_cases_map <- tm_shape(ghana_sh) +
  tm_fill(col = "Confirmed", style = "cont", palette = "Reds", title = "COVID-19 cases") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 cases \nacross ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_cases_map, filename = "outputs/ghana_Covid_19_cases_map.pdf")
tmap_save(ghana_cases_map, filename = "outputs/ghana_Covid_19_cases_map.jpg")

##Map  of ghana Recoveries
ghana_Recoveries_map <- tm_shape(ghana) +
  tm_fill(col = "Recoveries", style = "cont", palette = "Blues", title = "Recoveries") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Recoveries \nacross ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Recoveries_map, filename = "outputs/ghana_Covid_19_Recoveries_map.pdf")
tmap_save(ghana_Recoveries_map, filename = "outputs/ghana_Covid_19_Recoveries_map.jpg")


##Map  of ghana Deaths
ghana_Deaths_map <- tm_shape(ghana) +
  tm_fill(col = "Deaths", style = "cont", palette = "Reds", title = "Deaths") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Deaths \nacross ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Deaths_map, filename = "outputs/ghana_Covid_19_Deaths_map.pdf")
tmap_save(ghana_Deaths_map, filename = "outputs/ghana_Covid_19_Deaths_map.jpg")


##Map  of ghana Testing
ghana_Test_map <- tm_shape(ghana) +
  tm_fill(col = "Tests", style = "cont", palette = "Blues", title = "COVID-19 tests") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Testing \nacross ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Test_map, filename = "outputs/ghana_Covid_19_Testing_map.pdf")
tmap_save(ghana_Test_map, filename = "outputs/ghana_Covid_19_Testing_map.jpg")



##Map  of ghana sequences
ghana_sequences_map <- tm_shape(ghana) +
  tm_fill(col = "sequences", style = "cont", palette = "Blues", title = "Sequences in GISAID") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 sequences \nacross ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_sequences_map, filename = "outputs/ghana_Covid_19_sequences_map.pdf")
tmap_save(ghana_sequences_map, filename = "outputs/ghana_Covid_19_sequences_map.jpg")



##Map  of ghana Cases_1M
ghana_Cases_1M_map <- tm_shape(ghana) +
  tm_fill(col = "Cases_1M", style = "cont", palette = "Reds", title = "Cases per million") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 cases per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Cases_1M_map, filename = "outputs/ghana_Covid_19_Cases_1M_map.pdf")
tmap_save(ghana_Cases_1M_map, filename = "outputs/ghana_Covid_19_Cases_1M_map.jpg")



##Map  of ghana Tests_1M
ghana_Tests_1M_map <- tm_shape(ghana) +
  tm_fill(col = "Test_1M", style = "cont", palette = "Blues", title = "Tests per million") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Tests per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Tests_1M_map, filename = "outputs/ghana_Covid_19_Tests_1M_map.pdf")
tmap_save(ghana_Tests_1M_map, filename = "outputs/ghana_Covid_19_Tests_1M_map.jpg")



##Map  of ghana Positive_1KTests
ghana_Positive_1KTests_map <- tm_shape(ghana) +
  tm_fill(col = "Positive_1KTests", style = "cont", palette = "Blues", title = "Positive per 1,000 tests") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Tests per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Positive_1KTests_map, filename = "outputs/ghana_Covid_19_Positive_1KTests_map.pdf")
tmap_save(ghana_Positive_1KTests_map, filename = "outputs/ghana_Covid_19_Positive_1KTests_map.jpg")


##Map  of ghana Sequences_1HCases
ghana_Sequences_1HCases_map <- tm_shape(ghana) +
  tm_fill(col = "Sequences_1HCases", style = "cont", palette = "Blues", title = "Sequences per 100 cases") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Tests per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Sequences_1HCases_map, filename = "outputs/ghana_Covid_19_Sequences_1HCases_map.pdf")
tmap_save(ghana_Sequences_1HCases_map, filename = "outputs/ghana_Covid_19_Sequences_1HCases_map.jpg")


####Map  of ghana Recovery_1HCases
ghana_Recovery_1HCases_map <- tm_shape(ghana) +
  tm_fill(col = "Recovery_1HCases", style = "cont", palette = "Blues", title = "Recoveries per 100 cases") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Tests per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Recovery_1HCases_map, filename = "outputs/ghana_Covid_19_Recovery_1HCases_map.pdf")
tmap_save(ghana_Recovery_1HCases_map, filename = "outputs/ghana_Covid_19_Recovery_1HCases_map.jpg")


##Map  of ghana Deaths_1HCases
ghana_Deaths_1HCases_map <- tm_shape(ghana) +
  tm_fill(col = "Deaths_1HCases", style = "cont", palette = "Reds", title = "Deaths per 100 cases") +
  tm_borders(col = "black", lwd = 1) +
  tm_text("NAME_1") 
#tm_credits("COVID-19 Tests per 1M \npopulation across ghana", size = 0.8, position = c(0.6, 0.01))
#tm_layout(frame = FALSE, legend.width = 0.7, legend.position = c(0.85, 0.8), legend.title.size = 0.7) 
tmap_save(ghana_Deaths_1HCases_map, filename = "outputs/ghana_Covid_19_Deaths_1HCases_map.pdf")
tmap_save(ghana_Deaths_1HCases_map, filename = "outputs/ghana_Covid_19_Deaths_1HCases_map.jpg")

# Projected COVID-19 Cases in ghana based on no of positive tests per 1000 tests

# ghana_projected_cases_map <- tm_shape(ghana) + 
#  tm_fill(col = "Projected_cases", style = "cont", palette = "Reds", title = "Projected COVID-19 Cases in ghana") +
#  tm_borders(col = "black", lwd = 1) +
 # tm_text("NAME_1") 
#tmap_save(ghana_projected_cases_map, filename = "outputs/ghana_projected_cases_map.pdf")
#tmap_save(ghana_projected_cases_map, filename = "outputs/ghana_projected_cases_map.jpg")


#=============== END of MAPS ========================================



#================== Bar charts and other plots
# No data for this yet
# Regions that submitted sequences
ggplot(ghana_sequences, aes(y = reorder(State, sequences), x = sequences)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "SARS-CoV-2 sequences submitted in GISAID from ghana", y = "States", x = "Sequences")
ggsave("outputs/ghana_seq_submitted_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("outputs/ghana_seq_submitted_bar.pdf", height =  120, width = 200, units = "mm")


##Sequences submitted per 100 cases
##States that submitted sequences
ggplot(na.omit(ghana_cases_and_sequences), aes(y = reorder(State, Sequences_1HCases), x = Sequences_1HCases)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(Sequences_1HCases,2), "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "Sequences submitted per 100 COVID-19 cases from ghana", y = "States", x = "%")
ggsave("outputs/ghana_seq_submitted_per1H_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("outputs/ghana_seq_submitted_per1H_bar.pdf", height =  120, width = 200, units = "mm")


##Lineages circulating in ghana 
ggplot(Lineages, aes(y = reorder(Lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "SARS-CoV-2 lineages circulating in ghana", y = "Lineages", x = "Number of sequences")
ggsave("outputs/ghana_pango_lineages_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("outputs/ghana_pango_lineages_bar.pdf", height =  120, width = 200, units = "mm")



##  GISAID Clade diversity
ggplot(GISAID_clade, aes(y = reorder(Clade, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "GISAID clade diversity in ghana", y = "GISAID clades", x = "Number of sequences")
ggsave("outputs/ghana_gisaid_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("outputs/ghana_gisaid_bar.pdf", height =  120, width = 200, units = "mm")

# =========================== END OF SCRIPT=================================================
