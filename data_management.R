# Load the readxl package
library(readxl)
library(openxlsx)
library(dplyr)

##BIOSPECIMEN DATA
# Read the Excel file
biodata <- read_excel("biospecimenData_20240325.xlsx")
# Select the desired columns
biospecimen_data <- biodata[c("Participant ID", "External Id", "Biospecimens Id", "Tissue Type (Source Text)")]
# Rename the columns if needed
colnames(biospecimen_data) <- c("Participant_ID", "External_Id", "Biospecimens_Id", "Tissue_Type")
# Print
head(biospecimen_data)

##CLINICAL DATA
# Read the Excel file and select the 'Diagnoses' sheet
clindata <- read_excel("clinicalData_20240325.xlsx", sheet = "Diagnoses")
# Select the desired columns
clinical_data <- clindata[c("Participant ID", "External Id", "Diagnosis (Source Text)", "Age at Diagnosis (Days)", "Tumor Location")]
# Rename the columns if needed
colnames(clinical_data) <- c("Participant_ID", "External_Id", "Diagnosis", "Age_at_Diagnosis_Days", "Tumor_Location")
# Print
head(clinical_data)

##GENOMIC DATA
# Read the TSV file
gendata <- read.delim("participant-data-type.tsv", header = TRUE, stringsAsFactors = FALSE)
# Select the desired columns
genomic_data <- gendata[c("File.ID", "Participants.ID", "Family.Id", "Data.Type", "Data.Category", "File.Format")]
# Rename the columns if needed
colnames(genomic_data) <- c("Genome_File_ID", "Participant_ID", "Family_ID", "Data_Type", "Data_Category", "File_Format")
# Print
head(genomic_data)

##MERGE DATASETS
# Merge Biospecimen and Clinical datasets
merged_data <- merge(biospecimen_data, clinical_data, by = "Participant_ID", all = TRUE)
# Merge with Genomic dataset
merged_data <- merge(merged_data, genomic_data, by = "Participant_ID", all = TRUE)
# Selecting only the specified columns
dataset <- merged_data[c("Participant_ID", "External_Id.x", "Family_ID", "Diagnosis", "Tumor_Location", "Age_at_Diagnosis_Days", "Biospecimens_Id", "Tissue_Type", "Genome_File_ID", "Data_Type", "Data_Category", "File_Format")]
# Renaming columns for clarity
colnames(dataset) <- c("Participant_ID", "External_ID", "Family_ID", "Diagnosis", "Tumor_Location", "Age_at_Diagnosis_Days", "Biospecimens_ID", "Tissue_Type", "Genome_File_ID", "Data_Type", "Data_Category", "File_Format")
# Print
head(dataset)

##CASES WITH GENOMIC FILES
# Filter out rows with missing Genome_File_ID information
gen_dataset <- dataset[complete.cases(dataset$Genome_File_ID), ]
# Print the first few rows of the filtered dataset
head(gen_dataset)
# Count the total number of unique Participant_ID
total_participants <- length(unique(gen_dataset$Participant_ID))
# Print
print(total_participants) #406
#write out
write.xlsx(gen_dataset, "genomic_dataset.xlsx", rowNames = FALSE)

# Group by Participant_ID and summarize the data
summary_table <- gen_dataset %>%
  group_by(Participant_ID) %>%
  summarize(
    Total_Genome_File_ID = n_distinct(Genome_File_ID),
    Data_Types = paste(unique(Data_Type), collapse = ", ")
  )

# Print the summary table
print(summary_table)
write.xlsx(summary_table, "summary_table_genomic_dataset.xlsx", rowNames = FALSE)

# Create a table to summarize the counts of unique Data_Type
data_type_summary <- table(gen_dataset$Data_Type)
# Convert the table to a data frame for better readability
data_type_summary_df <- as.data.frame(data_type_summary)
# Rename the columns for clarity
colnames(data_type_summary_df) <- c("Data_Type", "Count")
# Print the summary table
print(data_type_summary_df)
write.xlsx(data_type_summary_df, "data_type_summary_genomict.xlsx", rowNames = FALSE)