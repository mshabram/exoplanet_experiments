# -*- coding: utf-8 -*-
"""
Created on Fri Apr 01 16:07:42 2016

@author: jcat
"""

def memoize(function):
    cache = {}
    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function
    
    
# show example usage
if __name__ == "__main__":
  import cProfile


def load_model_table_slow(loaded):
    import numpy as np
    if(not loaded):
        for j in range(100):
            model_table = np.genfromtxt('cdpp_slope_to_detection_efficiency_model_table.csv',delimiter=",")
        return model_table

@memoize
# load_model_table
def load_model_table_fast(loaded):
    import numpy as np
    if(not loaded):
        for j in range(100):
            model_table = np.genfromtxt('cdpp_slope_to_detection_efficiency_model_table.csv',delimiter=",")
        return model_table

print "Calling load_model_table_slow ..."
cProfile.run("load_model_table_slow(False)")

print "Calling load_model_table_fast ..."
cProfile.run("load_model_table_fast(True)")