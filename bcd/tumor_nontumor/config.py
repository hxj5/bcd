# config.py


import sys



class Config:
    def __init__(self):
        self.sid = None
        self.tool_list = None
        self.out_dir = None
        self.truth_fn = None
        self.tumor_labels = None
        self.overlap_how = "isec"
        self.fig_dpi = 300
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
        s += "%stumor_labels = %s\n" % (prefix, str(self.tumor_labels))
        s += "%soverlap_how = %s\n" % (prefix, self.overlap_how)
        s += "%sfig_dpi = %d\n" % (prefix, self.fig_dpi)
        s += "%sverbose = %s\n" % (prefix, str(self.verbose))
        s += "%s\n" % prefix

        fp.write(s)
