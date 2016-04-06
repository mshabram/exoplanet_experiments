# -*- coding: utf-8 -*-
class detection_efficiency_model( object ):
    """
    Given inputs of MES and CDPP slope, compute estimated detection efficiency,
    based on quadratic mapping between CDPPP slope and detection efficiency,
    as a function of MES

    If MES is outside range covered in table, extrapolate.
    Otherwise, interpolate between MES bins to get model detection 
    efficiency at desired MES.
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
    
