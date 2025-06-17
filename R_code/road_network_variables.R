#processing transportation files
folderpath<-"datafiles_new/Burghardt_BUA/CBSA_stats_cumulative_by_decade/"
# Note the relevant data can be downloaded here: https://figshare.com/articles/dataset/Historical_road_network_statistics_for_core-based_statistical_areas_in_the_U_S_1900_-_2010_/19584088
# the original data is not provided in this github project
files  <- list.files(folderpath)
files<-paste0(folderpath,files)
tables <- lapply(files, read.csv, header = TRUE)
road_network_vars <- do.call(rbind , tables)
road_network_vars<-subset(road_network_vars, year%in%c(2010, 2015))
road_network_vars$year_original<-road_network_vars$year
road_network_vars$year<-road_network_vars$year_original+1
  #ifelse(road_network_vars$year_original==2010, 2011, road_network_vars$year_original)
road_network_vars<-road_network_vars[c("msaid","year","patch_bupl", "patch_bua", "all_bupl", "all_bua","num_nodes","distance", "mean_local_gridness")]
rm(tables)
