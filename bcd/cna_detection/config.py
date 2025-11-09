# config.py


import sys



class Config:
    def __init__(self):
        self.sid = None
        self.tool_list = None
        self.out_dir = None
        self.truth_fn = None
        self.cell_anno_fn = None
        self.gene_anno_fn = None
        self.cna_type_list = None
        self.overlap_how = "isec-cells"
        self.max_n_cutoff = 1000
        self.fig_width = 4.25,
        self.fig_height = 3.25,
        self.fig_dpi = 300,
        self.fig_dec = 3,
        self.fig_legend_xmin = 0.5,
        self.fig_legend_ymin = 0.12,
        self.verbose = True


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
            
        s =  "%s\n" % prefix
        s += "%ssid = %s\n" % (prefix, self.sid)
        s += "%slen(tool_list) = %d\n" % (prefix, len(self.tool_list))
        s += "%stid list = '%s'\n" % (prefix, ", ".join([tool.tid for tool in self.tool_list]))
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%struth_fn = %s\n" % (prefix, self.truth_fn)
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sgene_anno_fn = %s\n" % (prefix, self.gene_anno_fn)
        s += "%scna_type_list = %s\n" % (prefix, str(self.cna_type_list))
        s += "%soverlap_how = %s\n" % (prefix, self.overlap_how)
        s += "%smax_n_cutoff = %s\n" % (prefix, str(self.max_n_cutoff))
        s += "%sfig_width = %f\n" % (prefix, self.fig_width)
        s += "%sfig_height = %f\n" % (prefix, self.fig_height)
        s += "%sfig_dpi = %d\n" % (prefix, self.fig_dpi)
        s += "%sfig_dec = %d\n" % (prefix, self.fig_dec)
        s += "%sfig_legend_xmin = %f\n" % (prefix, self.fig_legend_xmin)
        s += "%sfig_legend_ymin = %f\n" % (prefix, self.fig_legend_ymin)
        s += "%sverbose = %s\n" % (prefix, str(self.verbose))
        s += "%s\n" % prefix

        fp.write(s)
