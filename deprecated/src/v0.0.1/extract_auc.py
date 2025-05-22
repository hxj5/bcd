# extract_auc.py, converted from extract_auc.R

import pandas as pd
import os
import glob

def get_table1(cnv_scale, cnv_type, metric, fn):
    dat = pd.read_csv(fn, sep='\t')
    dat.columns.values[-1] = f"{cnv_type}_{dat.columns[-1]}"
    return dat

def get_table2(cnv_scale, metric, dir_path):
    metric_dir_id = "s5_roc" if metric == "roc" else "s6_prc"
    
    all_dat = None
    i = 1
    for cnv_type in ["copy_gain", "copy_loss", "loh"]:
        cnv_type_dir = os.path.join(dir_path, cnv_type, "result", metric_dir_id)
        if not os.path.exists(cnv_type_dir):
            cnv_type_dir = os.path.join(dir_path, cnv_type, metric_dir_id)
            if not os.path.exists(cnv_type_dir):
                raise FileNotFoundError(f"dir '{cnv_type_dir}' does not exist!")
        
        fn_list = glob.glob(os.path.join(cnv_type_dir, "*.auc.df.tsv"))
        fn = fn_list[0]
        
        dat = get_table1(cnv_scale=cnv_scale, cnv_type=cnv_type, metric=metric, fn=fn)
        
        if cnv_type == "loh":
            casper_rows = dat[dat['method'] == "casper"]
            if len(casper_rows) > 1:
                min_val = casper_rows.iloc[:, -1].min()
                for idx in range(len(dat)):
                    if dat.iloc[idx]['method'] == "casper" and dat.iloc[idx, -1] == min_val:
                        break
                dat = dat.drop(idx)
        
        if i == 1:
            all_dat = dat
        else:
            dat = dat.iloc[:, [0, -1]]
            all_dat = all_dat.merge(dat, on="method", how="left")
        i += 1
    
    return all_dat

def get_table3(cnv_scale, dir_path):
    roc_dat = get_table2(cnv_scale=cnv_scale, metric="roc", dir_path=dir_path)
    prc_dat = get_table2(cnv_scale=cnv_scale, metric="prc", dir_path=dir_path)
    all_dat = pd.concat([roc_dat, prc_dat.iloc[:, -2:]], axis=1)
    return all_dat