# args.py



class ToolArgs:
    def __init__(
        self, 
        tid = None, 
        has_gain = True, has_loss = True, has_loh = True
    ):
        self.tid = tid
        self.has_gain = has_gain
        self.has_loss = has_loss
        self.has_loh = has_loh
        
    def has_cna_type(self, cna_type):
        if cna_type == "gain":
            return(self.has_gain)
        elif cna_type == "loss":
            return(self.has_loss)
        elif cna_type == "loh":
            return(self.has_loh)
        else:
            return(None)


        
class InferCNVArgs(ToolArgs):
    def __init__(self, obj_fn):
        super().__init__(tid = "inferCNV", has_loh = False)
        self.obj_fn = obj_fn
        


class NumbatArgs(ToolArgs):
    def __init__(self, obj_fn):
        super().__init__(tid = "Numbat")
        self.obj_fn = obj_fn
