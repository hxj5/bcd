# base.py


import os
from ..utils.base import assert_e


class Tool:
    def __init__(
        self, 
        tid,
        obj_path,
        out_dir
    ):
        self.tid = tid
        self.obj_path = obj_path
        self.out_dir = out_dir
        
        self.__pp()
        
    def __pp(self):
        assert_e(self.obj_path)
        os.makedirs(self.out_dir, exist_ok = True)
