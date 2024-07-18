# Temporal-spatial-associations-of-large-groups-of-laying-hens-in-a-quasi-commercial-barn
Code for the manuscript "Temporal-spatial associations of large groups of laying hens in a quasi-commercial barn"

•	CombineGMMresultsOfAllWeeksOnePen.R
In this script, we combine the results the Gaussian Mixture Model (GMM) for each week, and we execute some analysis of interest, like extracting social networks or strong durable ties. This script must be run for each pen separately and the pen name must be specified at the beginning of the script. Note that this script is specifically for the analyses regarding the entire pen (interior part + wintergarden). To function the script needs several files, all are included in the data provided (folder data).
•	CombineGMMresultsOfAllWeeksOnePen_WG.R
This script is different from the previous one because it contains analysis specifically for the wintergarden. It must be run for each pen separately and the pen name must be specified at the beginning of the script. 
•	CombineGMMresultsOfAllWeeksOnePen_interior.R
It contains the analysis specifically for the interior part of the pen. It must be run for each pen separately and the pen name must be specified at the beginning of the script. Moreover, it compares the association matrices within the same week of the interior part and the wintergarden (which means that the result of the wintergarden analysis should already be available). 
•	GaussianMixtureModel_pen_week_singlepen_singleweek.R
This is the script that was used to calculate the association matrices for each week. We used a Gaussian Mixture Model from the package asnipe() – see main article for citations. The code is supplied but the data for running the code is not supplied huge amount of data. If interested, please contact the corresponding author. 
•	impact weather on wintergarden use.R
This scripts models the impact of weather – more specifically temperature and cloud coverage – on the proportion hens use the wintergarden. The script does this for all the pens together.  

DATA (supplied with the main manuscript)
There is one folder for each pen. Inside each pen’s folder, there are the folder ‘interior’, where the association matrices of the interior part of the pen are saved; ‘WG, instead, where the matrices of the wintergarden are saved. The folder ‘other’ is for other files that are necessary for the analysis of that specific pen. Data is mainly saved as R environment (*.RData). Each file contains data to build the social network for each week. The name of the file reports the pen, the week and which analysis is (if not specified, it is the entire pen, ‘WG’ for wintergarden and ‘interior’ for the interior part). The code provided puts these files together depending on which social networks you want to build. Just take care of setting the right working directory.
Inside the folder data there are also some files that are cumulative for all the pens.
