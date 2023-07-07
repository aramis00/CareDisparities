# This script can be used to take the event log and perform the test whether variant distributions are significantly
# different and to perform the computations for the aggregation of timing and sofa information in variants
# across classes
import itertools
import pickle
import process_mining_utils as pm
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
from scipy.stats import bootstrap
import numpy as np
def filter_icd(row,diags):
    for d in diags:
        if row[d] != 0:
            return True

    return False

def avg(nums,axis):
    return sum(nums)/len(nums)

log_raw = pd.read_csv("files/event_log.csv")
log_raw["starttime"] = pd.to_datetime(log_raw["starttime"])
log_raw["endtime"] = pd.to_datetime(log_raw["endtime"])

#We determined the most detrimental comorbidities and perform the analysis for each of those
for comorbs in [["congestive_heart_failure","myocardial_infarct","chronic_pulmonary_disease"]]:#[["age_score"]]:#,["renal_disease"],["congestive_heart_failure","myocardial_infarct","chronic_pulmonary_disease"]]:#[["renal_disease"]]:#,["mild_liver_disease","severe_liver_disease"],["congestive_heart_failure","myocardial_infarct","chronic_pulmonary_disease"],["malignant_cancer","malignant_cancer"]]:
    print("Calculating results for comorbidities:")
    print(comorbs)
    #The log can be filtered for low frequency variants, i.e., outliers. This can be done before or after applying the comorbidity filer
    log = pm.clean_low_frequency_variants(log_raw, 50)
    #The log is filtered for comorbidities here
    log["filter"] = log.apply(lambda x: filter_icd(x,comorbs),axis=1)
    #We could apply further filters
    #log = log[log["filter"]]
    #log["filter"] = log["charlson_comorbidity_index_r"].apply(lambda x: x >= 0 and x <= 4)
    log = log[log["filter"]]
    del log["filter"]

    print("Filtering Done")
    #print(log)

    #We calculate the process mining artifacts
    variant_dict = pm.log_to_variants_dict(log)
    class_frequency_dict = {}
    frequencies = pm.variants_dict_to_frequency_dist(variant_dict)
    #print(frequencies)

    # We generate a table with the variable for which we would like to generate statistics later
    stats_variables = ["anchor_age", "charlson_comorbidity_index", "gender", "ethnicity", "language"]
    stats_log = log_raw.groupby("stay_id")[stats_variables].agg("first")
    # print(stats_log)
    stats_log.to_csv("case_stats.csv", index=False)




    #The attributes for which we would like to split patients and investiage for disparities
    cols = ["anchor_year_group", "gender", "ethnicity", "language"]
    logs = {}

    #We initialize a dictionary for each attribute
    for c in cols:
        logs[c] = {}
        #control_dict[c] = {}
    for c in cols:
        vals = log[c].unique()
        logs[c] = {elem : pd.DataFrame() for elem in vals}
        for key in logs[c].keys():
            logs[c][key] = log[:][log[c] == key]
            logs[c][key]
    #print(logs)

    results_dict = {}

    #We conduct the comparison for each attribute. We determine the possible values of this attribute and compute the
    #statistics for each variant and each attribute value
    for c in cols:
        results_dict[c] = {}
        overall_var_dict = {}
        class_frequency_dict[c] = {}
        plt.rcParams["figure.figsize"] = (20, 3)
        dist = {}
        vals = list(logs[c].keys())
        #We determine the variants that are present for each attribute value
        for k in vals:
            variant_dict_c = pm.log_to_variants_dict(logs[c][k])
            overall_var_dict[k] = variant_dict_c
            dist[k] = pm.variants_dict_to_frequency_dist(variant_dict_c)


        #print(dist)
        #We determine the set of all variants contained in the event log by computing the union of variants across
        #attribute values
        all_variants = set()
        for k in vals:
            vars = list(dist[k].keys())
            all_variants = all_variants.union(set(vars))
        all_variants = list(all_variants)



        # Column for the dataframe to contain the frequency of variants across attribute values
        column_list = ["Variant"] + [v for v in vals]

        # rows for the frequency dataset
        lists = {i: [str(i)] for i in range(0, len(all_variants))}

        #constructing the frequency dataset
        for i in range(0,len(all_variants)):
            var = all_variants[i]
            for k in vals:
                if var in dist[k].keys():
                    lists[i].append(dist[k][var])
                else:
                    lists[i].append(0)
        draw_df = pd.DataFrame([lists[k] for k in lists.keys()],
                               columns=column_list)
        #Draw the distribution oover variants
        draw_df.plot(x="Variant",
                     kind='bar',
                     stacked=False,
                     title='Variant Distribution over ' + c)
        plt.savefig(c +comorbs[0]+'.png', dpi = 400)

        # We construct the frequency distribution to feed into the statistical test
        freq_list = {}
        for k in vals:
            freq = []
            for v in all_variants:
                if v in overall_var_dict[k].keys():
                    freq.append(len(overall_var_dict[k][v]))
                else:
                    freq.append(1)
            freq_list[k] = freq

        #for each combination of attribute value pairs, we test whether the distribution over variants is
        #significantly different
        for (v1, v2) in itertools.combinations(vals, r=2):
            cont_tab = [freq_list[v1],freq_list[v2]]
            #print(cont_tab)
            g, p, dof, expctd = chi2_contingency(cont_tab)
            #print(v1)
            #print(v2)
            print("DISTRIBUTION OVER VARIANTS: RESULTS for "+v1+" and "+v2)
            print(g)
            print("p-value: ")
            print(p)


        #We now aggregate the timing and sofa information in the variants
        for i in range(0, len(all_variants)):
            var = all_variants[i]
            results_dict[c][var] = {}
            class_frequency_dict[c][var] = {}
            (case, interval, variant_array, variant_string) = variant_dict[var][0]

            #We construct a schema for the differences we would like to time
            diffs = [ (i,i+1) for i in range(0,len(variant_array[0])-1)]
            diffs += [(interval[i].index('s'),interval[i].index('e')) for i in interval.keys()]
            #print(variant_array)

            #We collect the times and sofa scores for each of the differences and collect across all cases, i.e., patients
            diff_map_k = {}
            sofa_map_k = {}
            diff_map_k_dist = {}
            skip = False
            for k in vals:
                diff_map = {}
                sofa_map = {}
                if var not in overall_var_dict[k].keys():
                    skip = True
                    break

                results_dict[c][var][k] = {}

                #We retrieve all a cases of the variant for this attribute value
                case_list = overall_var_dict[k][var]
                class_frequency_dict[c][var][k] = len(case_list)

                #We iterate over each case and add the individual timing statistics to our aggregation
                for (case, interval, variant_array, variant_string) in case_list:
                    for (s,e) in diffs:
                        #find timestamps
                        s_time = None
                        e_time = None
                        for int_k in interval.keys():
                            if interval[int_k][s] != '':
                                if interval[int_k][s] == 's':
                                    s_time = log.loc[int_k]["starttime"]
                                else:
                                    s_time = log.loc[int_k]["endtime"]
                            if interval[int_k][e] != '':
                                if interval[int_k][e] == 's':
                                    e_time = log.loc[int_k]["starttime"]
                                else:
                                    e_time = log.loc[int_k]["endtime"]
                        if (s,e) not in diff_map.keys():
                            diff_map[(s,e)] = []
                        diff_map[(s,e)].append((e_time- s_time).total_seconds())

                    for e_id in interval.keys():
                        s_time = log.loc[e_id]["starttime"]
                        e_time = log.loc[e_id]["endtime"]
                        start_p = interval[e_id].index('s')
                        end_p = interval[e_id].index('e')
                        if (start_p, end_p) not in diff_map.keys():
                            diff_map[(start_p, end_p)] = []
                        if (start_p, end_p) not in sofa_map.keys():
                            sofa_map[(start_p, end_p)] = []
                        diff_map[(start_p, end_p)].append((e_time - s_time).total_seconds())
                        sofa_map[(start_p, end_p)].append(log.loc[e_id]["SOFA"])

                #We calculate the average over the collected metrics
                diff_map_k[k] = {(s,e): (-1,sum(diff_map[(s,e)])/len(diff_map[(s,e)])) for (s,e) in diffs}
                sofa_map_k[k] = {(s,e): (-1,sum(sofa_map[(s,e)])/len(sofa_map[(s,e)])) for (s,e) in sofa_map.keys()}

                diff_map_k_dist[k] = {(s,e) : diff_map[(s,e)] for (s,e) in diffs}
                for (s,e) in diffs:
                    sofa_av = (-1,-1)
                    # We bootstrap this calculation to get confidence intervals
                    if sum(diff_map[(s, e)]) > 0.00001 and len(diff_map[(s, e)])>10:
                        res = bootstrap((np.asarray(diff_map[(s, e)], dtype=np.float32),), np.mean, confidence_level=0.9, random_state=33)
                        #print(res.confidence_interval.low)
                        #print(res.confidence_interval.high)
                        diff_map_k[k][(s,e)] = (res.confidence_interval.low,res.confidence_interval.high)
                    elif len(diff_map[(s, e)])<0.00001:
                        diff_map_k[k][(s, e)] = (0, 0)

                    if (s,e) in sofa_map_k[k].keys():
                        sofa_av = sofa_map_k[k][(s,e)]
                        if len(sofa_map[(s, e)]) > 10:
                            res = bootstrap((np.asarray(sofa_map[(s, e)], dtype=np.float32),), np.mean,
                                            confidence_level=0.9, random_state=33)
                            #print(res.confidence_interval.low)
                            #print(res.confidence_interval.high)
                            sofa_av = (res.confidence_interval.low, res.confidence_interval.high)
                    results_dict[c][var][k][(s, e)] = (diff_map_k[k][(s,e)],sofa_av)
            if skip:
                if var in results_dict[c]:
                    del results_dict[c][var]
                continue
            #print(pd.DataFrame(diff_map_k))
            #analyze
            # This code could be used to perform hypothesis test on the difference between aggregated timing information
            # reg_df_l = []
            # for (s,e) in diffs:
            #
            #     for k in vals:
            #         for v in diff_map_k_dist[k][(s,e)]:
            #             age = log.loc[s]["admission_age"]
            #
            #             reg_df_l.append({"class":k,"time":v,"age":age})
            #     reg_df = pd.DataFrame(reg_df_l)
            #     one_hot = pd.get_dummies(reg_df["class"])
            #     one_hot_original_names = list(one_hot.columns.values)
            #     one_hot.columns = one_hot.columns.str.replace(" ", "_")
            #     one_hot.columns = one_hot.columns.str.replace("-", "_")
            #     one_hot.columns = one_hot.columns.str.replace("_", "")
            #     #one_hot.columns = [ "Q('"+col+"')" for col in one_hot.columns]
            #     # Drop column B as it is now encoded
            #     reg_df = reg_df.drop("class", axis=1)
            #     # Join the encoded df
            #     reg_df = reg_df.join(one_hot)
            #     #print(reg_df)
            #     #reg_df.columns = reg_df.columns.str.replace(" ", "_")
            #     #reg = LinearRegression().fit(reg_df["class"], reg_df["time"])
            #     result = sm.ols(formula="time ~ " + " + ".join([ "Q('"+col+"')" for col in one_hot.columns])+ " + age", data=reg_df).fit()
            #     print("(" + str(s) + "," + str(e) + ")")
            #     print(result.summary())
            #     for k in vals:
            #         results_dict[c][var][k][(s, e)] = (results_dict[c][var][k][(s, e)], result.pvalues.loc["Q('"+list(one_hot.columns.values)[one_hot_original_names.index(k)]+"')"])
            #     #result = sm.ols(x = one_hot, y= reg_df["time"]).fit()

    #print(results_dict)

    #We write the results of our aggregation to different files

    with open('result'+comorbs[0]+'.pkl', 'wb') as fp:
        pickle.dump(results_dict, fp)

    with open('variants'+comorbs[0]+'.pkl', 'wb') as fp:
        pickle.dump(variant_dict, fp)

    with open('frequency'+comorbs[0]+'.pkl', 'wb') as fp:
        pickle.dump(class_frequency_dict, fp)
