# -*- coding: utf-8 -*-
"""
Created on Fri Apr 01 16:07:42 2016

@author: jcat
"""
class detection_efficiency_model( object ):
    """
    Given inputs of MES and CDPP slope, and the function load_model_table,
    compute estimated detection efficiency,
    based on quadratic mapping between CDPPP slope and detection efficiency,
    as a function of MES

    If MES is outside range covered in table, extrapolate.
    Otherwise, interpolate between MES bins to get model detection 
    efficiency at desired MES.
    
    !!!!! NOTE: It is up to the user to ensure that the input CDPP slope
    is in the valid range of -0.6 to + 0.4
    
    Clean up the model by setting negative values to zero and values 
    greater than one to one.
    """

    def __init__(self,MES,CDPPslope,load_model_table):
        self.model_table = load_model_table(True)
        self.MES = MES
        self.CDPPslope = CDPPslope
    
    def detection_efficiency(self):
        MESlist = self.model_table[:,0]
        MES = self.MES
        CDPPslope = self.CDPPslope
        
        # Extrapolate when the MES value is outside
        #   the range covered by the table
        
        # If MES < minimum MES, assume detection efficiency
        #   is 0
        if MES < MESlist[0]:
            detection_efficiency = 0
            
        # If MES > maximum MES, assume detection efficiency 
        #   has the same value as at maximum MES
        elif MES > MESlist[-1]:
            detection_efficiency = self.model_table[-1,1] + self.model_table[-1,2]*CDPPslope + self.model_table[-1,3]*CDPPslope**2
 
        # Interpolate when the MES value is within the range
        #   coverered by the table:
        else:
            from numpy import interp as interp
            detection_efficiency_all = self.model_table[:,1] + self.model_table[:,2]*CDPPslope + self.model_table[:,3]*CDPPslope**2
            MES_all = self.model_table[:,0]   
            detection_efficiency = interp(MES,MES_all,detection_efficiency_all)

        # Clean up negative values and values that exceed 1
        if detection_efficiency < 0:
            detection_efficiency = 0
        elif detection_efficiency > 1:
            detection_efficiency = 1
            
        # Return computed detection efficiency
        return detection_efficiency

def memoize(function):
# from http://programmingzen.com/2009/05/18/memoization-in-ruby-and-python/
    cache = {}
    def memoized_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return memoized_function
    
    
# show example usage
if __name__ == "__main__":
  import cProfile
  
# Memoize, so that first call to load_model_table (inside detection_efficiency_model.py) has to load the table
# But on subsequent calls the table is cached!
@memoize
def load_model_table(True):
    import numpy as np
    model_table = np.genfromtxt('cdpp_slope_to_detection_efficiency_model_table.csv',delimiter=",")
    return model_table

# Timings show that first call has to load_model_table has to load the table
# But on subsequent calls the table is cached!
cProfile.run("de1 = detection_efficiency_model( 6.125,-0.1,load_model_table )")
cProfile.run("de2 = detection_efficiency_model( 9.111,-0.1,load_model_table )")
cProfile.run("de3 = detection_efficiency_model( 10.98,-0.1,load_model_table )")
cProfile.run("de4 = detection_efficiency_model( 2,-0.1,load_model_table )")
cProfile.run("de5 = detection_efficiency_model( 30,-0.1,load_model_table )")

# Results of test runs
print 'detection efficiency for MES = 6.125, CDPP slope = - 0.1 is ' + str(de1.detection_efficiency())
print 'detection efficiency for MES = 9.211, CDPP slope = - 0.1 is ' + str(de2.detection_efficiency())
print 'detection efficiency for MES = 10.981, CDPP slope = - 0.1 is ' + str(de3.detection_efficiency())
print 'detection efficiency for MES = 2, CDPP slope = - 0.1 is ' + str(de4.detection_efficiency())
print 'detection efficiency for MES = 30, CDPP slope = - 0.1 is ' + str(de5.detection_efficiency())

# Second time the table is in cache
# print "Calling load_model_table ..."
# cProfile.run("model_table2 = load_model_table(True)")
#cProfile.run("de2 = detection_efficiency_model( 6.125,-0.1,load_model_table(True) )")

