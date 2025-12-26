# predict.py - predict subclonal structure.


import os
from logging import info
from ..utils.base import assert_e



def run_predict(
    tool_list, out_dir, 
    k,
    verbose = True
):
    # check args.
    assert len(tool_list) > 0
    os.makedirs(out_dir, exist_ok = True)
    
    out_fn_list = []
    for tool in tool_list:
        tid = tool.tid.lower()
        info("predict subclonal structure for '%s' ..." % tid)

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)
        out_fn = os.path.join(res_dir, "%s_predictions.tsv" % tid)
        
        if tid == "calicost":
            tool.predict(
                out_fn = out_fn,
                verbose = verbose
            )

        elif tid == "copykat":
            tool.predict(
                out_fn = out_fn,
                k = k,
                tmp_dir = res_dir,
                verbose = verbose
            )
        
        elif tid == "infercnv":
            out_fn = tool.predict(
                res_dir,
                k = k,
                dist = 'euclidean',
                hclust = 'ward.D2',
                verbose = verbose
            )

        elif tid == "numbat":
            tool.predict(
                out_fn,
                verbose = verbose     
            )
                
        elif tid == "xclone":
            pass
        
        else:
            raise ValueError(f"Error: unknown tool id '{tid}'.")

        out_fn_list.append(out_fn)


    res = dict(
        # out_fns : list of str
        #   Output subclonal structure prediction files.
        #   Each value is a prediction TSV file, in the same order
        #   with `tool_list`.
        out_fns = out_fn_list
    )
    
    return(res)
