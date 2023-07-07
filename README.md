# Uncovering Care Disparities Using Process Mining
This repository provides the code to run a disparity analysis using process mining techniques. It is applied to patients within MIMIC-IV that have sepsis but can be extended to any other event log. To extract the event log of sepsis patients from MIMIC please
1. Apply the sql statements in mimic_extraction/ to MIMIC-IV and store the resulting tables.
2. Apply the Full_preprocessing_notebook.ipynb to retrieve the event log.

After that, you can run the analysis on the event log. You can do so by:
1. Running the main.py script to perform the hypothesis test on the distribution of variants across patient attributes and to generate the aggregated variants with timing and SOFA information.
2. Applying the statistics.py script to generate descriptive statistics and testing for their distribution across patient groups.
3. Applying the visualization.py script to generate the visualizations of variants, along with the aggregated timing and SOFA information for all patient groups.
