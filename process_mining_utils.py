######################
# We provide process minning utility functions in this file, allowing users to construct the different process
# mining artifacts from cases to variants.
######################
import numpy as np



def case_to_interval(case):
    arr = case.to_numpy()
    indexes = list(case.index)
    arr = np.c_[np.c_[indexes],arr]
    start_tuples = [(arr[i][0],"s",arr[i][3]) for i in range(0,len(arr))]
    end_tuples = [(arr[i][0],"e",arr[i][4])for i in range(0,len(arr))]
    all_tuples = [(e,s,t) for (e,s,t) in start_tuples + end_tuples if case.loc[e]["activity"]!="ADMIT" and case.loc[e]["activity"]!="DISCHARGE"]
    all_tuples.sort(key = lambda x:x[2])
    all_tuples+=[(e,s,t) for (e,s,t) in start_tuples + end_tuples if case.loc[e]["activity"]=="DISCHARGE"]
    all_tuples = [(e, s, t) for (e, s, t) in start_tuples + end_tuples if case.loc[e]["activity"] == "ADMIT"] + all_tuples
    interval_case = {arr[x][0]:["" for y in range(0,len(arr)*2)] for x in range(0,len(arr))}
    for counter in range(0,len(all_tuples)):
        (i, ty, ti) = all_tuples[counter]
        interval_case[i][counter] = ty
    return interval_case

def interval_to_variant(interval, case):
    variant_array = []
    keys = interval.keys()
    acts = [case.loc[i][1] for i in keys]
    keys, acts = zip(*sorted(zip(keys, acts)))
    for k in keys:
        lane_list = []
        for elem in interval[k]:
            if elem != '':
                lane_list.append(case.loc[k][1])
            else:
                lane_list.append("None")
        variant_array.append(lane_list)
    return variant_array


def log_to_variants_dict(log):
    cases = log.groupby("stay_id")
    cases = [x for _, x in cases]
    variant_dict = {}
    for case in cases:
        interval = case_to_interval(case)
        variant_array = interval_to_variant(interval, case)
        variant_string = ''.join(["".join(lane) for lane in variant_array])
        if variant_string not in variant_dict:
            variant_dict[variant_string] = []
        variant_dict[variant_string].append((case, interval, variant_array, variant_string))
    #print("N of vars: "+str(len(variant_dict.keys())))
    return variant_dict

def variants_dict_to_frequency_dist(variant_dict):
    frequencies = {key: len(variant_dict[key]) for key in variant_dict.keys()}
    total_freq = sum([v for k,v in frequencies.items() ])
    rel_freq = {key: frequencies[key]/total_freq for key in frequencies.keys()}
    return rel_freq

def clean_low_frequency_variants(log,thresh):
    variant_dict = log_to_variants_dict(log)
    to_delete = []
    for v in variant_dict.keys():
        if len(variant_dict[v]) < thresh:
            for (case, interval, variant_array, variant_string) in variant_dict[v]:
                to_delete += list(case.index)
    log = log.drop(to_delete)
    return log