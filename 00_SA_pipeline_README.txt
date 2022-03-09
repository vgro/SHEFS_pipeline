# READ ME for South Africa SDM pipeline
# author: Vivienne Groner
# date: 09.03.2022
# source: my_path


DATA PRE-PROCESSING, SDM RUN ON CLUSTER, POST-PROCESSING

1. prepare GBIF occurence records (download, clean, rarify, create background and pseudoabsences, extract environmental data)
	01_SA_prep_occurrence_2021.R (loop over list of species in folders, /04_occurrence_records/GBIF_raw)

2. create sh and R jobs for cluster
	02_job_maker.R
	SA_dummy.R
	SA_dummy.sh
	03_make_batch_files.R

3. run SDM on cluster (present (1979-2013) and RCP45, RCP85 for 5 GCMs)
	species specific R and sh file
	batch files

	SA_functions_CBER.R
	02_environmental_data/CHELSA/
	02_environmental_data/CHELSA_future/
	04_occurrence_records/GBIF_ready/environmental/
	04_occurrence_records/GBIF_ready/points/

4. copy tar.zip files from myriad to RDS 
	(login myriad -> login xxxxxxx@transfer02)

5. untar folders
	04_untar.R
	cmd < cmd_call.txt

6. reorder SDM output in new structure (-> I copied core folders to /data/)


LAND USE AND INTENSITY

1. resample land cover data to 1 km resolution (reclassifed to PREDICTS and high level land cover classes in ArcMap)
	05_SA_NLC_2020_resample.R

2. create land use and intensity suitability map
	calculate effect sizes for each group (arthropods, birds, amphibia, reptiles, mammals, gastropods)
		06_PREDICTS_subset.R 
		
	create map for 2020 and 2080 ssp245 and ssp585
		07_create_LCC_scenarios_CropNatUrbWat.R


DATA ANALYSIS

SA_functions_CBER_analysis.R

read SDM output and create dataframe with presence locations
	08_SDM_df_climate_current.R
	08_SDM_df_climate_future.R

and calculate LUI effect:
	08_SDM_df_climate_lui_current.R
	08_SDM_df_climate_lui_future.R


SENSITIVITY TESTS

distribution of AUCs for all species
	code: 09_sensitivity_analysis_AUC.R

sensitivity test of buffer size and density for pseudoabsences
	code: 10_sensitivity_analysis_pseudoabsence.R
        output: /04_occurrence_records/buffer_test/

sensitivity test LUI factors
	code: 11_sensitivity_analysis_LUI_factors.R
	output: /02_environmental_data/LUI_sensitivity/




