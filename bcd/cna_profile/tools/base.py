# base.py



class Tool:
    def __init__(
        self, 
        tid = None,
        has_gain = True, has_loss = True, has_loh = True
    ):
        self.tid = tid
        self.__has_cna_type = {
            'gain': has_gain,
            'loss': has_loss,
            'loh': has_loh
        }


    def has_cna_type(self, cna_type):
        if cna_type in self.__has_cna_type:
            return(self.__has_cna_type[cna_type])
        else:
            raise ValueError("invalid CNA type '%s'" % cna_type)
