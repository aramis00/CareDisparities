import pandas as pd
import numpy as np
from scipy.stats import iqr
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sns
res = {}
variants = {}
controls = {}
frequencies = {}

sens_char = "gender"

stats = pd.read_csv("case_stats.csv")
print(stats.columns)
print(stats)
for c in ["anchor_age"]:#,"charlson_comorbidity_index"]:
    sns.histplot(data=stats, x=c, hue=sens_char)
    plt.show()
    print(c)
    grouped =stats.groupby(sens_char)[c].agg(["median",iqr])
    print(grouped)
    native = stats[stats[sens_char] == "F"][c]#'Native english speaker'][c]
    nonnative = stats[stats[sens_char] =="M"][c]# 'Non-native speaker'][c]
    #Mann-Whitney
    stat, p_value = mannwhitneyu(native, nonnative)
    print('Mann-Whitney U Test Statistics=%.3f, p=%.3f' % (stat, p_value))


for c in ["ethnicity","language"]:#"language"
    grouped = stats.groupby(sens_char)[c].value_counts()
    print(grouped)
    grouped = stats.groupby(sens_char)[c].value_counts(normalize=True)
    print(grouped)
    # Create a contingency table
    contingency_table = pd.crosstab(stats[sens_char], stats[c])
    # Perform the chi-square test
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    print("Chi-square statistic =", chi2)
    print("p-value =", p)

