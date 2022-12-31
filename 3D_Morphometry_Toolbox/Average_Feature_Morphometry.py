## IMPORT MODULES ##

import arcpy
import pandas as pd
import numpy as np
import math
import os

print("Modules imported")
arcpy.AddMessage("Modules imported")

# ------------------------------------------------------------
## DEFINE VARIABLES

env_workspace = arcpy.GetParameterAsText(3)
transect_data = arcpy.GetParameterAsText(0)
feature_data = arcpy.GetParameterAsText(1)
outline_data = arcpy.GetParameterAsText(2)
output_folder = env_workspace

print("Variables defined")
arcpy.AddMessage("Variables defined")

# --------------------------------------------------------------
## CREATE FILE GEODATABASE (fdgb)

FGDB = arcpy.management.CreateFileGDB(output_folder, "fgdb", "CURRENT")

print("File geodatabase created")
arcpy.AddMessage("File geodatabase created")

# -----------------------------------------------------------
## CONVERT TRANSECT DATA SHAPEFILE TO FEATURE CLASS

Transect_Morphometry_FC = arcpy.FeatureClassToFeatureClass_conversion(transect_data, FGDB, "Transect_Morphometry_FC")

# ------------------------------------------------------------
## EXPORT Transect_Data ATTRIBUTE DATA TO EXCEL FILE

# Define local variables
in_table = Transect_Morphometry_FC
out_xls = output_folder+"\\"+"Transect_Morphometry.xls"

# Run conversion
arcpy.TableToExcel_conversion(in_table, out_xls)

print("Transect data exported to excel file")
arcpy.AddMessage("Transect data exported to excel")

# -----------------------------------------------------------
## READ / WRITE FILE ##

input_file = pd.read_excel(out_xls)
df = input_file

# ------------------------------------------------------------
## CREATE NEW COLUMN 'AV_HEIGHT'

arcpy.AddField_management(transect_data, "AV_HEIGHT", "LONG")

print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## AVERAGE TRANSECT 'HEIGHT'

def average_feature(group):
    g = group['HEIGHT'].mean()
    group['AV_HEIGHT'] = g

    return group

df = df.groupby('Feature_ID').apply(average_feature)
df.head()

print("Feature height averaged")
arcpy.AddMessage("Feature height averaged")

# -------------------------------------------------------------
## CREATE NEW COLUMN 'AV_ASYMMETRY'

arcpy.AddField_management(transect_data, "AV_ASYMM", "LONG")

print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## AVERAGE TRANSECT 'ASYMMETRY'

def average_feature(group):
    g = group['ASYMMETRY'].mean()
    group['AV_ASYMM'] = g

    return group

df = df.groupby('Feature_ID').apply(average_feature)
df.head()

print("Feature asymmetry averaged")
arcpy.AddMessage("Feature asymmetry averaged")

# -------------------------------------------------------------
## CREATE NEW COLUMN 'AV_SLOPE'

arcpy.AddField_management(transect_data, "AV_SLOPE", "LONG")

print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## AVERAGE TRANSECT 'SLOPE'

def average_feature(group):
    g = group['Avg_Slope'].mean()
    group['AV_SLOPE'] = g

    return group

df = df.groupby('Feature_ID').apply(average_feature)
df.head()

print("Feature slope averaged")
arcpy.AddMessage("Feature slope averaged")

# -------------------------------------------------------------
## CREATE NEW COLUMN 'AV_WIDTH'

arcpy.AddField_management(transect_data, "AV_WIDTH", "LONG")

print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## AVERAGE TRANSECT 'WIDTH'

def average_feature(group):
    g = group['WIDTH'].mean()
    group['AV_WIDTH'] = g

    return group

df = df.groupby('Feature_ID').apply(average_feature)
df.head()

print("Feature width averaged")
arcpy.AddMessage("Feature width averaged")

# -------------------------------------------------------------
## CREATE NEW COLUMN 'TOTAL_VOLUME'

arcpy.AddField_management(transect_data, "TOTAL_VOL", "LONG")

print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## SUM TRANSECT 'CROSS-SECTIONAL VOLUME'

def average_feature(group):
    g = group['CS_VOL'].sum()
    group['TOTAL_VOL'] = g

    return group

df = df.groupby('Feature_ID').apply(average_feature)
df.head()

print("Feature volume calculated")
arcpy.AddMessage("Feature volume calculated")

# -------------------------------------------------------------
## REMOVE DUPLICATE ROWS

df = df.drop_duplicates('Feature_ID', keep='first')

# -------------------------------------------------------------
## SAVE DF AS NEW EXCEL(.XLS) file

df.to_csv(output_folder+"\\"+"Feature_Data.csv")

print("Exported to new excel file")
arcpy.AddMessage("Exported to new excel file")

# -------------------------------------------------------------
## JOIN AVERAGED TRANSECT DATA TO FEATURE DATA

# set local variables
inFeatures = feature_data
joinTable = output_folder+"\\"+"Feature_Data.csv"
joinField = "Feature_ID"

Feature_MorphometryCopy = arcpy.CopyRows_management(joinTable, output_folder+"\\"+"Feature_MorphometryCopy")

Feature_Morphometry_FL = arcpy.MakeFeatureLayer_management(inFeatures, "Feature_Morphometry_FL")
arcpy.JoinField_management(Feature_Morphometry_FL, joinField, Feature_MorphometryCopy, joinField)
Feature_MorphometryJOIN = arcpy.CopyFeatures_management(Feature_Morphometry_FL, output_folder+"\\"+"Av_Feature_Morphometry1")

print("Averaged transect data joined to shapefile")
arcpy.AddMessage("Averaged transect data joined to shapefile")

# -------------------------------------------------------------
## JOIN AVERAGED FEATURE DATA TO OUTLINE SHAPEFILE

# set local variables
inFeatures = outline_data
joinTable = output_folder+"\\"+"Av_Feature_Morphometry1.shp"
joinField = "Feature_ID"

OutlineCopy = arcpy.CopyRows_management(joinTable, output_folder+"\\"+"OutlineCopy")

Outline_FL = arcpy.MakeFeatureLayer_management(inFeatures, "Outline_FL")
arcpy.JoinField_management(Outline_FL, joinField, joinTable, joinField, ["SINUOSITY", "LENGTH", "CL_Z", "AV_BASE_TE", "AV_HEIGHT", "AV_ASYMM", "AV_SLOPE", "AV_WIDTH", "TOTAL_VOL"])
Feature_MorphometryJOIN = arcpy.CopyFeatures_management(Outline_FL, output_folder+"\\"+"Av_Feature_Morphometry")

print("Averaged Feature Morphometry shapefile created")
arcpy.AddMessage("Averaged Feature Morphometry shapefile created")


# ----------------------------------------------------------------------------------------
## DELETE UNWANTED COLUMNS - AV_FEATURE_MORPHOMETRY ##

# Set local variables
inFeatures = Feature_MorphometryJOIN
Av_Feature_Morphometry = output_folder+"\\"+"Av_Feature_Morphometry.shp"
dropFields = ["OBJECTID_1", "FEATURE__2", ]

# Execute DeleteField
arcpy.DeleteField_management(Av_Feature_Morphometry, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

# ----------------------------------------------------------------------------------------
## DELETE UNWANTED COLUMNS - TRANSECT_MORPHOMETRY ##

# Set local variables
dropFields = ["AV_HEIGHT", "AV_ASYMM", "AV_SLOPE", "AV_WIDTH", "TOTAL_VOL"]

# Execute DeleteField
arcpy.DeleteField_management(transect_data, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

# ----------------------------------------------------------------------------------------
## DELETE UNWANTED COLUMNS - FEATURE_OUTLINE ##

# Set local variables
dropFields = ["LENGTH", "SINUOSITY", "CL_Z", "AV_BASE_TE", "AV_HEIGHT", "AV_ASYMM", "AV_SLOPE", "AV_WIDTH", "TOTAL_VOL"]

# Execute DeleteField
arcpy.DeleteField_management(outline_data, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

# ----------------------------------------------------------------------------------------
## DELETE UNWANTED COLUMNS - FEATURE_MORPHOMETRY ##

# Set local variables
dropFields = ["OBJECTID", "FIELD1", "OBJECTID1", "OBJECTID_1", "CL_Z", "AVG_SLOPE", "FEATURE__1", "AV_BASE_TE", "HEIGHT", "CS_AR", "CS_VOL", "ASYMMETRY", "WIDTH", "SHAPE_LENG", "AV_HEIGHT", "AV_ASYMM", "AV_SLOPE", "AV_WIDTH", "TOTAL_VOL"]

# Execute DeleteField
arcpy.DeleteField_management(feature_data, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")
# ----------------------------------------------------------------------------------------
## DELETE UNWANTED FILES

arcpy.Delete_management(FGDB)
arcpy.Delete_management(output_folder+"\\"+"Feature_MorphometryCopy")
arcpy.Delete_management(output_folder+"\\"+"OutlineCopy")
arcpy.Delete_management(output_folder+"\\"+"Av_Feature_Morphometry1.shp")

os.remove(output_folder+"\\"+"Transect_Morphometry.xls")
os.remove(output_folder+"\\"+"Feature_Data.csv")

print("Unnecessary files removed")
arcpy.AddMessage("Unnecessary files removed")

# -----------------------------------------------------------------------------------------

print("Script completed with no errors")
arcpy.AddMessage("Script completed with no errors")

# -----------------------------------------------------------------------------------------


