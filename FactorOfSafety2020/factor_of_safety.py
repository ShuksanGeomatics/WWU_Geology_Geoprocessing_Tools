# -*- coding: utf-8 -*-
'''
Calculates Factor of Safety using Infinite Slope (a slope stability model outlined by Hammond et al. (1992))
'''
__author__ = "Gerry Gabrisch GISP"
__date__ = "December 3, 2020"
__copyright__ = "MIT"
__credits__ = "Robert J  Mitchell PhD, Western Washington University, Bellingham, WA"
__version__ = "Compatible with both Python 3.6/ArcPro v2.6 geoprocessing tool or Python 2.7/ArcGIS 10.7."

import traceback
import sys

try:
       
        import arcpy
        from arcpy.sa import *
        arcpy.AddMessage('Running factor of safety using Infinite Slope per Hammond et al. (1992)')
        arcpy.AddMessage('Western Washington University Department of Geology, 12/03/2020')
        arcpy.AddMessage('Checking out required extensions.')
        arcpy.CheckOutExtension("spatial")
        arcpy.CheckOutExtension("3D")
        
        #######################                          UNCOMMENT FOR TESTING                   #######################
        ###User Defined Environment Factors 
        ## To prevent overwriting outputs change overwriteOutput option to False.
        #arcpy.env.overwriteOutput = True
        ##raster output cell size in coordinate reference system units of measure.  For example, 10 meters if inputs are in UTM.
        #arcpy.env.cellSize = "10"
        ##Set the geoprocessing extent.  User can enter values or use the maximum extent of the input rasters.
        ##arcpy.env.extent = "549036.905972784 5395667.22147943 555546.905972784 5400767.22147943"
        #arcpy.env.extent = "MAXOF"
        #arcpy.env.workspace = r"E:\InfiniteSlope\Lab8\Lab8\InfiniteSlope2014\SmithCreekFOS.gdb"        
       
        ##User Defined Parameters...
        #Dw_D_ratio = 0.9058680904958534
        #tree_root_strength_as_cohesion = 2		# Cs (kPa)
        #vegetation_surcharge =	0.5			# qo (kPa)
        #soil_cohesion =	15.0				# Cs(kPa)
        #effective_internal_angle_of_friction = 25.0	# phi (degrees o)
        #dry_soil_unit_weight =	14.5			# gamma dry(kN/m3)
        #moist_soil_unit_weight = 16.0			# gamma moist(kN/m3)
        #saturated_soil_unit_weight =	18.0		# gamma sat(kN/m3)
        #water_unit_weight = 9.81			# gamma water(kN/m3)
        
        
        ##The name of the input slope raster (with slope in degrees) from the workspace geodatabase...
        #surface_slope_raster_in_degrees = r"slope_degree"
        ##The name of the input total soil thickness raster from the workspace geodatabase...
        #total_soil_thickness_raster = r"soildepth_high"
        
        ##The user defined name of the output factor of safety raster. This will be the resulting name in the workspace geodatabase...
        #factor_of_safety = r"factor_of_safety_Gerry20201203"
        #########################################################################################################
        
        
        #inputs
        arcpy.env.workspace = arcpy.GetParameterAsText(0)
        Dw_D_ratio = arcpy.GetParameter(1)
        tree_root_strength_as_cohesion = arcpy.GetParameter(2)
        vegetation_surcharge =	arcpy.GetParameter(3)	
        soil_cohesion =	arcpy.GetParameter(4)		
        effective_internal_angle_of_friction = arcpy.GetParameter(5)
        dry_soil_unit_weight =	arcpy.GetParameter(6)
        moist_soil_unit_weight = arcpy.GetParameter(7)
        saturated_soil_unit_weight = arcpy.GetParameter(8)
        water_unit_weight = arcpy.GetParameter(9)
        surface_slope_raster_in_degrees = arcpy.GetParameterAsText(10)
        total_soil_thickness_raster =arcpy.GetParameterAsText(11)
        #output
        factor_of_safety = arcpy.GetParameterAsText(12)   
        
        
        # To prevent overwriting outputs change overwriteOutput option to False.
        arcpy.env.overwriteOutput = arcpy.GetParameterAsText(13)
        #raster output cell size in coordinate reference system units of measure.  For example, 10 meters if inputs are in UTM.
        
        
        arcpy.env.cellSize = "10"
        #Set the geoprocessing extent.  User can enter values or use the maximum extent of the input rasters.
        #arcpy.env.extent = "549036.905972784 5395667.22147943 555546.905972784 5400767.22147943"
        arcpy.env.extent = "MAXOF"        
        
        def FactorSafety(surface_slope_raster_in_degrees, total_soil_thickness_raster, Dw_D_ratio, moist_soil_unit_weight, dry_soil_unit_weight, saturated_soil_unit_weight, water_unit_weight, vegetation_surcharge, tree_root_strength_as_cohesion, soil_cohesion, effective_internal_angle_of_friction, factor_of_safety):  
                
                #Constant for the conversion from degrees slope to slope in radians - pi/180
                constant_value = 0.01745329252
                
                #convert all inputs to floats to avoid integer division rounding issues with Python
                tree_root_strength_as_cohesion = float(tree_root_strength_as_cohesion)
                vegetation_surcharge =	float(vegetation_surcharge)
                soil_cohesion =	float(soil_cohesion)
                effective_internal_angle_of_friction = float(effective_internal_angle_of_friction)
                dry_soil_unit_weight =	float(dry_soil_unit_weight)
                moist_soil_unit_weight = float(moist_soil_unit_weight)
                saturated_soil_unit_weight =	float(saturated_soil_unit_weight)
                water_unit_weight = float(water_unit_weight)
               
                
                #Convert slope in degrees raster to slope in radians...
                slope_radians = arcpy.sa.Times(surface_slope_raster_in_degrees, constant_value)
                slope_radians.save("slope_radians")
                
                #Calculate the Dw  = soil thickness raster * Dw/D ratio..
                Dw = 'Dw'
                Dw = arcpy.sa.Times(total_soil_thickness_raster, Dw_D_ratio)        
                Dw.save('Dw')
        
        
                #Calculat the resisting force (numerator of the FOS equation)...
                arcpy.AddMessage('Calculating the resisting force raster...')
                resisiting_forces = "resisiting_forces"
                
                resisiting_forces = tree_root_strength_as_cohesion + soil_cohesion + Square(Cos(slope_radians)) * (vegetation_surcharge + moist_soil_unit_weight * (total_soil_thickness_raster - Dw) + (saturated_soil_unit_weight - water_unit_weight) *Dw) * Tan(effective_internal_angle_of_friction * constant_value)
                resisiting_forces.save('resisiting_forces')
        
                #Calculate the driving force (denominator of the FOS equation)...
                arcpy.AddMessage('Calculating the driving force raster...')
                driving_forces = "driving_forces"
                driving_forces = Sin(slope_radians) * Cos(slope_radians) * (vegetation_surcharge + moist_soil_unit_weight * (total_soil_thickness_raster - Dw) + saturated_soil_unit_weight * Dw)
                driving_forces.save('driving_forces')
        
                #Calculate the factor of safety (resisting force/driving force)...
                arcpy.AddMessage('Calculating the factor of safety')
                factor_of_safetytemp = 'factor_of_safetytemp'
                factor_of_safetytemp = arcpy.sa.Divide(resisiting_forces, driving_forces)
                factor_of_safetytemp.save(factor_of_safety)
                
                
        FactorSafety(surface_slope_raster_in_degrees, total_soil_thickness_raster, Dw_D_ratio, moist_soil_unit_weight, dry_soil_unit_weight, saturated_soil_unit_weight, water_unit_weight, vegetation_surcharge, tree_root_strength_as_cohesion, soil_cohesion, effective_internal_angle_of_friction, factor_of_safety)
        arcpy.AddMessage('Finished without errors!')

except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        # Return python error messages for use in script tool or Python window
        arcpy.AddError(pymsg)
        arcpy.AddError(msgs)
        # Print Python error messages for use in Python / Python window
        print(pymsg)
        print(msgs)
