### GIS_SCRIPT FOR TRANSECT ANALYSIS ###

## Import modules ##

import arcpy
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

print("Modules imported")
arcpy.AddMessage("Modules imported")

## Define variables
env_workspace = arcpy.GetParameterAsText(5)
outlines = arcpy.GetParameterAsText(0)
crestlines = arcpy.GetParameterAsText(1)
DEM = arcpy.GetParameterAsText(2)
Transect_Intervals = arcpy.GetParameterAsText(3)
Transect_Length = arcpy.GetParameterAsText(4)
output_folder = env_workspace

print("Variables defined")
arcpy.AddMessage("Variables defined")

## Create file geodatabase (fdgb) ##

FGDB = arcpy.management.CreateFileGDB(output_folder, "fgdb", "CURRENT")

print("File geodatabase create")
arcpy.AddMessage("File geodatabase created")

## Calculate Sinuosity ##

arcpy.AddField_management(crestlines, "SINUOSITY", "DOUBLE")

import math

sinuosity_expression = '''
import math
def getSinuosity(shape):
    length = shape.length
    d = math.sqrt((shape.firstPoint.X - shape.lastPoint.X) ** 2 +
                  (shape.firstPoint.Y - shape.lastPoint.Y) ** 2)
    return d / length
'''

arcpy.CalculateField_management(crestlines,
                                "SINUOSITY",
                                'getSinuosity(!shape!)',
                                'PYTHON_9.3',
                                 sinuosity_expression)


## Clip crestlines to extent of outlines ##

# Define variables
Clipped_CLs = output_folder+"\\"+"Clipped_CLs.shp"

# Run 'Clip' tool
arcpy.Clip_analysis(crestlines, outlines, Clipped_CLs)

print("Crestlines clipped")
arcpy.AddMessage("Crestlines clipped")

## Add ID & geometry fields ##

arcpy.AddField_management(outlines, "Feature_ID", "SHORT")
arcpy.AddField_management(outlines, "AREA", "LONG")
arcpy.AddField_management(Clipped_CLs, "CL_ID", "SHORT")
arcpy.AddField_management(Clipped_CLs, "LENGTH", "LONG")


print("New fields have been added")
arcpy.AddMessage("New fields have been added")

## Calculate ID & geometry fields ##

arcpy.management.CalculateField(outlines, "Feature_ID", '!FID!', "PYTHON_9.3")
arcpy.management.CalculateField(Clipped_CLs, "CL_ID", '!FID!', "PYTHON_9.3")
arcpy.management.CalculateField(outlines, "AREA", '!shape.area@meters!', "PYTHON_9.3")
arcpy.management.CalculateField(Clipped_CLs, "LENGTH", '!shape.length@meters!', "PYTHON_9.3")

print("ID fields have been calculated")
arcpy.AddMessage("ID fields have been calculated")
print("Geometry fields have been calculated")
arcpy.AddMessage("Geometry fields have been calculated")
print("Sinuosity has been calculated")
arcpy.AddMessage("Sinuosity has been calculated")

## Join outlines shapefile to crestlines, based on spatial location, closest to ##

# Define variables
CL_JOIN1 = output_folder+"\\"+"MERGED_CL_OL1.shp"

# Execute function
arcpy.SpatialJoin_analysis(Clipped_CLs, outlines, CL_JOIN1, "#", "#", "#", "CLOSEST", "#", "DISTANCE")

print("Join completed")
arcpy.AddMessage("Join completed")

## Delete unwanted columns ##

# Set local variables
inFeatures = CL_JOIN1
CL_JOIN = output_folder+"\\"+"Feature_Morphometry.shp"
dropFields = ["Join_Count", "DISTANCE", "TARGET_FID", "Id", "CL_ID", "Id_1"]

# Execute CopyFeatures to make a new copy of the feature class
#  Use CopyRows if you have a table
arcpy.CopyFeatures_management(inFeatures, CL_JOIN)

# Execute DeleteField
arcpy.DeleteField_management(CL_JOIN, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

## Run 'Transect' tool ##

print("Transect toolbox imported")
arcpy.AddMessage("Running transect tool")

############################### Import transect script #######################################

# Import system modules
import arcpy
from arcpy import env
import math

arcpy.env.overwriteOutput = True

# Set environments
arcpy.env.overwriteOutput = True
arcpy.env.XYResolution = "0.00001 Meters"
arcpy.env.XYTolerance = "0.0001 Meters"

# Set local variables
env.workspace = env_workspace
Lines = CL_JOIN
SplitType = "Split at approximate distance"
DistanceSplit = float(Transect_Intervals)
TransecLength = Transect_Length
TransecLength_Unit = "METERS"
OutputTransect = CL_JOIN = output_folder+"\\"+"Transects.shp"


# Def splitline module
###START SPLIT LINE CODE IN A SAME DISTANCE### Source: http://nodedangles.wordpress.com/2011/05/01/quick-dirty-arcpy-batch-splitting-polylines-to-a-specific-length/
def splitline(inFC, FCName, alongDist):
    OutDir = env.workspace
    outFCName = FCName
    outFC = OutDir + "/" + outFCName

    def distPoint(p1, p2):
        calc1 = p1.X - p2.X
        calc2 = p1.Y - p2.Y

        return math.sqrt((calc1 ** 2) + (calc2 ** 2))

    def midpoint(prevpoint, nextpoint, targetDist, totalDist):
        newX = prevpoint.X + ((nextpoint.X - prevpoint.X) * (targetDist / totalDist))
        newY = prevpoint.Y + ((nextpoint.Y - prevpoint.Y) * (targetDist / totalDist))
        return arcpy.Point(newX, newY)

    def splitShape(feat, splitDist):
        # Count the number of points in the current multipart feature
        #
        partcount = feat.partCount
        partnum = 0
        # Enter while loop for each part in the feature (if a singlepart feature
        # this will occur only once)
        #
        lineArray = arcpy.Array()

        while partnum < partcount:
            # Print the part number
            #
            # print "Part " + str(partnum) + ":"
            part = feat.getPart(partnum)
            # print part.count

            totalDist = 0

            pnt = part.next()
            pntcount = 0

            prevpoint = None
            shapelist = []

            # Enter while loop for each vertex
            #
            while pnt:

                if not (prevpoint is None):
                    thisDist = distPoint(prevpoint, pnt)
                    maxAdditionalDist = splitDist - totalDist

                    print("thisDist, totalDist, maxAdditionalDist")

                    if (totalDist + thisDist) > splitDist:
                        while (totalDist + thisDist) > splitDist:
                            maxAdditionalDist = splitDist - totalDist
                            # print thisDist, totalDist, maxAdditionalDist
                            newpoint = midpoint(prevpoint, pnt, maxAdditionalDist, thisDist)
                            lineArray.add(newpoint)
                            shapelist.append(lineArray)

                            lineArray = arcpy.Array()
                            lineArray.add(newpoint)
                            prevpoint = newpoint
                            thisDist = distPoint(prevpoint, pnt)
                            totalDist = 0

                        lineArray.add(pnt)
                        totalDist += thisDist
                    else:
                        totalDist += thisDist
                        lineArray.add(pnt)
                        # shapelist.append(lineArray)
                else:
                    lineArray.add(pnt)
                    totalDist = 0

                prevpoint = pnt
                pntcount += 1

                pnt = part.next()

                # If pnt is null, either the part is finished or there is an
                #   interior ring
                #
                if not pnt:
                    pnt = part.next()
                    if pnt:
                        print("Interior Ring:")
            partnum += 1

        if (lineArray.count > 1):
            shapelist.append(lineArray)

        return shapelist

    if arcpy.Exists(outFC):
        arcpy.Delete_management(outFC)

    arcpy.Copy_management(inFC, outFC)

    # origDesc = arcpy.Describe(inFC)
    # sR = origDesc.spatialReference

    # revDesc = arcpy.Describe(outFC)
    # revDesc.ShapeFieldName

    deleterows = arcpy.UpdateCursor(outFC)
    for iDRow in deleterows:
        deleterows.deleteRow(iDRow)

    try:
        del iDRow
        del deleterows
    except:
        pass

    inputRows = arcpy.SearchCursor(inFC)
    outputRows = arcpy.InsertCursor(outFC)
    fields = arcpy.ListFields(inFC)

    numRecords = int(arcpy.GetCount_management(inFC).getOutput(0))
    OnePercentThreshold = numRecords // 100

    # printit(numRecords)

    iCounter = 0
    iCounter2 = 0

    for iInRow in inputRows:
        inGeom = iInRow.shape
        iCounter += 1
        iCounter2 += 1
        if (iCounter2 > (OnePercentThreshold + 0)):
            # printit("Processing Record "+str(iCounter) + " of "+ str(numRecords))
            iCounter2 = 0

        if (inGeom.length > alongDist):
            shapeList = splitShape(iInRow.shape, alongDist)

            for itmp in shapeList:
                newRow = outputRows.newRow()
                for ifield in fields:
                    if (ifield.editable):
                        newRow.setValue(ifield.name, iInRow.getValue(ifield.name))
                newRow.shape = itmp
                outputRows.insertRow(newRow)
        else:
            outputRows.insertRow(iInRow)

    del inputRows
    del outputRows

    # printit("Done!")


###END SPLIT LINE CODE IN A SAME DISTANCE###

# Create "General" file geodatabase
WorkFolder = env.workspace
General_GDB = WorkFolder + "\General.gdb"
arcpy.CreateFileGDB_management(WorkFolder, "General", "CURRENT")
env.workspace = General_GDB

# Unsplit Line
LineDissolve = "LineDissolve"
arcpy.Dissolve_management(Lines, LineDissolve, "", "", "SINGLE_PART")
LineSplit = "LineSplit"

# Split Line
if SplitType == "Split at approximate distance":
    splitline(LineDissolve, LineSplit, DistanceSplit)
else:
    arcpy.SplitLine_management(LineDissolve, LineSplit)

# Add fields to LineSplit
FieldsNames = ["LineID", "Direction", "Azimuth", "X_mid", "Y_mid", "AziLine_1", "AziLine_2", "Distance"]
for fn in FieldsNames:
    arcpy.AddField_management(LineSplit, fn, "DOUBLE")

# Calculate Fields
CodeBlock_Direction = """def GetAzimuthPolyline(shape):
 radian = math.atan((shape.lastpoint.x - shape.firstpoint.x)/(shape.lastpoint.y - shape.firstpoint.y))
 degrees = radian * 180 / math.pi
 return degrees"""

CodeBlock_Azimuth = """def Azimuth(direction):
 if direction < 0:
  azimuth = direction + 360
  return azimuth
 else:
  return direction"""
CodeBlock_NULLS = """def findNulls(fieldValue):
    if fieldValue is None:
        return 0
    elif fieldValue is not None:
        return fieldValue"""
arcpy.CalculateField_management(LineSplit, "LineID", "!OBJECTID!", "PYTHON_9.3")
arcpy.CalculateField_management(LineSplit, "Direction", "GetAzimuthPolyline(!Shape!)", "PYTHON_9.3",
                                CodeBlock_Direction)
arcpy.CalculateField_management(LineSplit, "Direction", "findNulls(!Direction!)", "PYTHON_9.3", CodeBlock_NULLS)
arcpy.CalculateField_management(LineSplit, "Azimuth", "Azimuth(!Direction!)", "PYTHON_9.3", CodeBlock_Azimuth)
arcpy.CalculateField_management(LineSplit, "X_mid", "!Shape!.positionAlongLine(0.5,True).firstPoint.X", "PYTHON_9.3")
arcpy.CalculateField_management(LineSplit, "Y_mid", "!Shape!.positionAlongLine(0.5,True).firstPoint.Y", "PYTHON_9.3")
CodeBlock_AziLine1 = """def Azline1(azimuth):
 az1 = azimuth + 90
 if az1 > 360:
  az1-=360
  return az1
 else:
  return az1"""
CodeBlock_AziLine2 = """def Azline2(azimuth):
 az2 = azimuth - 90
 if az2 < 0:
  az2+=360
  return az2
 else:
  return az2"""
arcpy.CalculateField_management(LineSplit, "AziLine_1", "Azline1(!Azimuth!)", "PYTHON_9.3", CodeBlock_AziLine1)
arcpy.CalculateField_management(LineSplit, "AziLine_2", "Azline2(!Azimuth!)", "PYTHON_9.3", CodeBlock_AziLine2)
arcpy.CalculateField_management(LineSplit, "Distance", TransecLength, "PYTHON_9.3")

# Generate Azline1 and Azline2
spatial_reference = arcpy.Describe(Lines).spatialReference
Azline1 = "Azline1"
Azline2 = "Azline2"
arcpy.BearingDistanceToLine_management(LineSplit, Azline1, "X_mid", "Y_mid", "Distance", TransecLength_Unit,
                                       "AziLine_1", "DEGREES", "GEODESIC", "LineID", spatial_reference)
arcpy.BearingDistanceToLine_management(LineSplit, Azline2, "X_mid", "Y_mid", "Distance", TransecLength_Unit,
                                       "AziLine_2", "DEGREES", "GEODESIC", "LineID", spatial_reference)

# Create Azline and append Azline1 and Azline2
Azline = "Azline"
arcpy.CreateFeatureclass_management(env.workspace, "Azline", "POLYLINE", "", "", "", spatial_reference)
arcpy.AddField_management(Azline, "LineID", "DOUBLE")
arcpy.Append_management([Azline1, Azline2], Azline, "NO_TEST")

# Dissolve Azline
Azline_Dissolve = "Azline_Dissolve"
arcpy.Dissolve_management(Azline, Azline_Dissolve, "LineID", "", "SINGLE_PART")

# Add Fields to Azline_Dissolve
FieldsNames2 = ["x_start", "y_start", "x_end", "y_end"]
for fn2 in FieldsNames2:
    arcpy.AddField_management(Azline_Dissolve, fn2, "DOUBLE")

# Calculate Azline_Dissolve fields
arcpy.CalculateField_management(Azline_Dissolve, "x_start", "!Shape!.positionAlongLine(0,True).firstPoint.X",
                                "PYTHON_9.3")
arcpy.CalculateField_management(Azline_Dissolve, "y_start", "!Shape!.positionAlongLine(0,True).firstPoint.Y",
                                "PYTHON_9.3")
arcpy.CalculateField_management(Azline_Dissolve, "x_end", "!Shape!.positionAlongLine(1,True).firstPoint.X",
                                "PYTHON_9.3")
arcpy.CalculateField_management(Azline_Dissolve, "y_end", "!Shape!.positionAlongLine(1,True).firstPoint.Y",
                                "PYTHON_9.3")

# Generate output file
arcpy.XYToLine_management(Azline_Dissolve, OutputTransect, "x_start", "y_start", "x_end", "y_end", "", "",
                          spatial_reference)

# Delete General.gdb
arcpy.Delete_management(General_GDB)

############################################# END OF TRANSECT SCRIPT #############################

print("'Transect' tool finished")
arcpy.AddMessage("Transect tool finished")

## Run 'Create Points on Lines' tool ##

############################################ CREATE POINTS ON LINE SCRIPT #############################

arcpy.env.overwriteOutput = True

line = OutputTransect
create_from = "Beginning"
choice = "PERCENTAGE"
use_field = "NO"
field = ""
distance = float(0.5)
end_points = "NO"
CLPoutput = output_folder+"\\"+"CL_POINTS.shp"
CL_POINTS = CLPoutput

if "in_memory" in CLPoutput:
    mem_name = CLPoutput.split("\\")[-1]
else:
    mem_name = "mem_point"

mem_point = arcpy.CreateFeatureclass_management("in_memory", mem_name, "POINT", "", "DISABLED", "DISABLED", line)
arcpy.AddField_management(mem_point, "LineOID", "TEXT")
arcpy.AddField_management(mem_point, "Value", "TEXT")

result = arcpy.GetCount_management(line)
features = int(result.getOutput(0))

arcpy.SetProgressor("step", "Creating Points on Lines...", 0, features, 1)

fields = ["SHAPE@", "OID@"]

if use_field == "YES":
    fields.append(field)

reverse = False
if create_from == "END":
   reverse = True
   reversed_points = []

with arcpy.da.SearchCursor(line, (fields)) as search:
    with arcpy.da.InsertCursor(mem_point, ("SHAPE@", "LineOID", "Value")) as insert:
        for row in search:
            try:
                line_geom = row[0]
                length = float(line_geom.length)
                count = distance
                oid = str(row[1])
                start = arcpy.PointGeometry(line_geom.firstPoint)
                end = arcpy.PointGeometry(line_geom.lastPoint)

                if reverse == True:
                   for part in line_geom:
                       for p in part:
                           reversed_points.append(p)

                   reversed_points.reverse()
                   array = arcpy.Array([reversed_points])
                   line_geom = arcpy.Polyline(array)

                if use_field == "YES":
                    count = float(row[2])
                    distance = float(row[2])

                if choice == "DISTANCE":
                    point = line_geom.positionAlongLine(count, False)
                    insert.insertRow((point, oid, count))

                elif choice == "INTERVAL":
                    while count <= length:
                      point = line_geom.positionAlongLine(count, False)
                      insert.insertRow((point, oid, count))
                      count += distance

                elif choice == "PERCENTAGE":
                    point = line_geom.positionAlongLine(count, True)
                    insert.insertRow((point, oid, count))

                elif choice == "START/END POINTS":
                    insert.insertRow((start, oid, 0))
                    insert.insertRow((end, oid, str(length)))

                if end_points == "START":
                    insert.insertRow((start, oid, 0))

                elif end_points == "END":
                    insert.insertRow((end, oid, str(length)))

                elif end_points == "BOTH":
                    insert.insertRow((start, oid, 0))
                    insert.insertRow((end, oid, str(length)))

                arcpy.SetProgressorPosition()

            except Exception as e:
                arcpy.AddMessage(str(e.message))

if "in_memory" in CLPoutput:
    arcpy.SetParameter(8, mem_point)
else:
    arcpy.CopyFeatures_management(mem_point, CLPoutput)
    arcpy.Delete_management(mem_point)

arcpy.ResetProgressor()
arcpy.GetMessages()

############################################ END OF CPOL SCRIPT #######################################

print("'Create Points on Lines' tool finished")
arcpy.AddMessage("Create Points on Lines tool finished")

## Add elevation (z) information to CL_Points ##

# Run Add Surface Information 3D Tool
arcpy.CheckOutExtension("3D")
arcpy.AddSurfaceInformation_3d(CL_POINTS, DEM, "Z")

print("Elevation data added")
arcpy.AddMessage("Elevation data added")

## Convert shapefile to feature class

CL_POINTS_FC = arcpy.FeatureClassToFeatureClass_conversion(CL_POINTS, FGDB, "CL_POINTS_FC")

## Change 'Z' field name ##

# Define variables
fc = CL_POINTS_FC
field = "Z"
new_name = "CL_Z"
new_alias = "CL_Z"
new_type = "Double"
new_length = 5
new_is_nullable = "NULLABLE"
clear_alias = "FALSE"

# Run 'AlterField_management' tool
arcpy.AlterField_management (fc, field, new_name, new_alias, new_type, new_length, new_is_nullable, clear_alias)

print("Field name changed")
arcpy.AddMessage("Field name changed")

## Clip transects to extent of outlines ##

# Define variables
Clipped_Transects = output_folder+"\\"+"Clipped_Transects.shp"

# Run 'Clip' tool
arcpy.Clip_analysis(OutputTransect, outlines, Clipped_Transects)

print("Transects clipped")
arcpy.AddMessage("Transects clipped")

## Explode clipped transects ##

# Define variables
Exploded_Transects = output_folder+"\\"+"Exploded_Transects.shp"

# Run 'Multipart to Singlepart' tool
arcpy.MultipartToSinglepart_management (Clipped_Transects, Exploded_Transects)

print("'Clipped_Transects' exploded")
arcpy.AddMessage("Clipped transects exploded")

## Join the CL_Points shapefile to the Exploded_Transects shapefile ##

# Define variables
EXPLODED_JOIN = output_folder+"\\"+"EXPLODED_JOIN.shp"

# Execute function
arcpy.SpatialJoin_analysis(Exploded_Transects, CL_POINTS_FC, EXPLODED_JOIN, "#", "#", "#", "INTERSECT", 0, "#")

print("Join completed")
arcpy.AddMessage("Join completed")

## Delete unwanted transect segments based on intersection with CL_POINTS ##

with arcpy.da.UpdateCursor(EXPLODED_JOIN, "Join_Count") as cursor:
    for row in cursor:
        if row[0] < 1:
            cursor.deleteRow()

print("Unwanted transect segments removed")
arcpy.AddMessage("Unwanted transect segments removed")

## Add average slope information to EXPLODED_JOIN ##

# Run 'Add Surface Information 3D Tool'
arcpy.CheckOutExtension("3D")
arcpy.AddSurfaceInformation_3d(EXPLODED_JOIN, DEM, "AVG_SLOPE")

print("Slope data added")
arcpy.AddMessage("Slope data added")

## Run 'Create Points on Lines' tool ##

print("Toolbox imported")
arcpy.AddMessage("Create Points on Lines tool imported")

############################################ CREATE POINTS ON LINE SCRIPT #############################

arcpy.env.overwriteOutput = True

line = EXPLODED_JOIN
create_from = "#"
choice = "START/END POINTS"
use_field = "#"
field = "#"
distance = "#"
end_points = "#"
TPoutput = output_folder+"\\"+"TRANSECT_POINTS.shp"
TRANSECT_POINTS = TPoutput

if "in_memory" in TPoutput:
    mem_name = TPoutput.split("\\")[-1]
else:
    mem_name = "mem_point"

mem_point = arcpy.CreateFeatureclass_management("in_memory", mem_name, "POINT", "", "DISABLED", "DISABLED", line)
arcpy.AddField_management(mem_point, "LineOID", "TEXT")
arcpy.AddField_management(mem_point, "Value", "TEXT")

result = arcpy.GetCount_management(line)
features = int(result.getOutput(0))

arcpy.SetProgressor("step", "Creating Points on Lines...", 0, features, 1)

fields = ["SHAPE@", "OID@"]

if use_field == "YES":
    fields.append(field)

reverse = False
if create_from == "END":
   reverse = True
   reversed_points = []

with arcpy.da.SearchCursor(line, (fields)) as search:
    with arcpy.da.InsertCursor(mem_point, ("SHAPE@", "LineOID", "Value")) as insert:
        for row in search:
            try:
                line_geom = row[0]
                length = float(line_geom.length)
                count = distance
                oid = str(row[1])
                start = arcpy.PointGeometry(line_geom.firstPoint)
                end = arcpy.PointGeometry(line_geom.lastPoint)

                if reverse == True:
                   for part in line_geom:
                       for p in part:
                           reversed_points.append(p)

                   reversed_points.reverse()
                   array = arcpy.Array([reversed_points])
                   line_geom = arcpy.Polyline(array)

                if use_field == "YES":
                    count = float(row[2])
                    distance = float(row[2])

                if choice == "DISTANCE":
                    point = line_geom.positionAlongLine(count, False)
                    insert.insertRow((point, oid, count))

                elif choice == "INTERVAL":
                    while count <= length:
                      point = line_geom.positionAlongLine(count, False)
                      insert.insertRow((point, oid, count))
                      count += distance

                elif choice == "PERCENTAGE":
                    point = line_geom.positionAlongLine(count, True)
                    insert.insertRow((point, oid, count))

                elif choice == "START/END POINTS":
                    insert.insertRow((start, oid, 0))
                    insert.insertRow((end, oid, str(length)))

                if end_points == "START":
                    insert.insertRow((start, oid, 0))

                elif end_points == "END":
                    insert.insertRow((end, oid, str(length)))

                elif end_points == "BOTH":
                    insert.insertRow((start, oid, 0))
                    insert.insertRow((end, oid, str(length)))

                arcpy.SetProgressorPosition()

            except Exception as e:
                arcpy.AddMessage(str(e.message))

if "in_memory" in TPoutput:
    arcpy.SetParameter(8, mem_point)
else:
    arcpy.CopyFeatures_management(mem_point, TPoutput)
    arcpy.Delete_management(mem_point)

arcpy.ResetProgressor()
arcpy.GetMessages()

############################################ END OF CPOL SCRIPT #######################################

print("'Create Points on Lines' tool finished")
arcpy.AddMessage("Create Points on Lines tool finished")

## Add elevation (z) information to TRANSECT_POINTS ##

# Run Add Surface Information 3D Tool
arcpy.CheckOutExtension("3D")
arcpy.AddSurfaceInformation_3d(TRANSECT_POINTS, DEM, "Z")

print("Elevation data added")
arcpy.AddMessage("Elevation data added")

## Convert shapefile to feature class

TRANSECT_POINTS_FC = arcpy.FeatureClassToFeatureClass_conversion(TRANSECT_POINTS, FGDB, "TRANSECT_POINTS_FC")

## Change 'Z' field name ##

# Define local variables
fc = TRANSECT_POINTS_FC
field = "Z"
new_name = "Transect_Z"
new_alias = "Transect_Z"
new_type = "Double"
new_length = 5
new_is_nullable = "NULLABLE"
clear_alias = "FALSE"

# Run AlterField_management tool
arcpy.AlterField_management (fc, field, new_name, new_alias, new_type, new_length, new_is_nullable, clear_alias)

print("Field name changed")
arcpy.AddMessage("Field name changed")

## Join the 'EXPLODED_JOIN' to the 'TRANSECT_POINTS' shapefile ##

# Define variables
Transect_Points_Join = output_folder+"\\"+"Transect_Points_Join.shp"

# Execute function
arcpy.SpatialJoin_analysis(TRANSECT_POINTS_FC, EXPLODED_JOIN, Transect_Points_Join, "#", "#", "#", "CLOSEST", "#", "DISTANCE")

print("Join completed")
arcpy.AddMessage("Join completed")

## Join the original 'CL_JOIN' to the 'Transect_Points_Join' shapefile ##

# Define variables
TP_CL_Join = output_folder+"\\"+"TP_CL_Join.shp"

# Execute function
arcpy.SpatialJoin_analysis(Transect_Points_Join, CL_JOIN, TP_CL_Join, "#", "#", "#", "CLOSEST", "#", "DISTANCE")

print("Join completed")
arcpy.AddMessage("Join completed")

## Extract asymmetry data ##

# Run the 'Split Lines at Point' tool to split the clipped transects at the 'CL_POINTS'

# Define local variables
SPLIT_TRANSECTS = output_folder+"\\"+"SPLIT_TRANSECTS.shp"

# Execute tool
arcpy.SplitLineAtPoint_management(EXPLODED_JOIN, CL_POINTS_FC, SPLIT_TRANSECTS, '#')

print("Transect split completed")
arcpy.AddMessage("Transect split completed")

## Join the 'TP_CL_Join' shapefile to the 'SPLIT_TRANSECTS' shapefile ##

# Define variables
Final_Join = output_folder+"\\"+"Final_Join.shp"

# Execute function
arcpy.SpatialJoin_analysis(SPLIT_TRANSECTS, TP_CL_Join, Final_Join, "#", "#", "#", "CLOSEST", "#", "DISTANCE")

print("Join completed")
arcpy.AddMessage("Join completed")

# -----------------------------------------------------------------------------------------

## Delete unwanted columns ##

# Set local variables
inFeatures = Final_Join
Transect_Data1 = output_folder+"\\"+"Transect_Data1.shp"
dropFields = ["Join_Count", "DISTANCE", "TARGET_FID", "Join_Count", "TARGET_FID", "ORIG_FID", "LineOID", "Value", "Join_Count_1", "DISTANCE", "TARGET_FID_1", "Join_Count", "DISTANCE", "TARGET_FID", "LineOID_1", "Join_Count", "TARGET_FID", "x_start_1", "y_start_1", "x_end_1", "y_end_1", "ORIG_FID_1", "LineOID_12", "Value_12", "Avg_Slope_1", "Join_Count_12", "DISTANCE_1", "TARGET_FID_12", "X_start_12", "X_end_12", "Y_start_12", "Y_end_12", "Id"]

# Execute CopyFeatures to make a new copy of the feature class
#  Use CopyRows if you have a table
arcpy.CopyFeatures_management(inFeatures, Transect_Data1)

# Execute DeleteField
arcpy.DeleteField_management(Transect_Data1, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

# -------------- JOIN OUTLINES TO TRANSECT --------------------------------- #
# Define variables
Transect_Data2 = output_folder+"\\"+"Transect_Data2.shp"

# Execute function
arcpy.SpatialJoin_analysis(Transect_Data1, outlines, Transect_Data2, "JOIN_ONE_TO_MANY", "#", "#", "HAVE_THEIR_CENTER_IN", "#", "")

print("Join completed")
arcpy.AddMessage("Join completed")

## Delete unwanted columns ##

# Set local variables
inFeatures = Transect_Data2
Transect_Data3 = output_folder+"\\"+"Transect_Data3.shp"
dropFields = ["Join_Count", "TARGET_FID", "JOIN_FID", "Join_Cou_1", "TARGET_F_1", "Join_Cou_2", "Join_Cou_3", "DISTANCE_2", "TARGET_F_3", "Join_Cou_4", "TARGET_F_4", "TARGET_F_2", "CL_Z_1", "Avg_Slop_1", "Id", "AREA"]

# Execute CopyFeatures to make a new copy of the feature class
#  Use CopyRows if you have a table
arcpy.CopyFeatures_management(inFeatures, Transect_Data3)

# Execute DeleteField
arcpy.DeleteField_management(Transect_Data3, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")

## Convert final shapefile to feature class

Transect_Data_FC = arcpy.FeatureClassToFeatureClass_conversion(Transect_Data3, FGDB, "Transect_Data_FC")

#-----------------------------------------------------------------------------------------
## Re-merge split transects ##

# Define variables
Transect_Data = output_folder+"\\"+"Transect_Data.shp"

arcpy.management.Dissolve(Transect_Data3, Transect_Data, ['x_start', 'y_start', 'x_end', 'y_end', 'CL_Z', 'Avg_Slope', 'Feature_ID'], "", "SINGLE_PART", "DISSOLVE_LINES")


# -----------------------------------------------------------------------------------------

## Export Transect_Data attribute data to excel file ##

# Define local variables
in_table = Transect_Data_FC
out_xls = output_folder+"\\"+"Transect_Data.xls"

# Run conversion
arcpy.TableToExcel_conversion(in_table, out_xls)

print("Transect data exported to excel")
arcpy.AddMessage("Transect data exported to excel")

# -----------------------------------------------------------------------------------------
## Calculate morphometry

# Read and write excel file
read_file = pd.read_excel(out_xls)
read_file.to_csv(output_folder+"\\"+"Transect_Data_1.csv", index = None, header=True)
df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Data_1.csv"))

# Calculate average base terrain
df.insert(9, 'Av_Base_Te', '')
df['Av_Base_Te'] = df['Transect_Z'].rolling(window=2).mean()

print("Average base terrain calculated")
arcpy.AddMessage("Average base terrain calculated")

# Calculate height
df.insert(10, 'Height', '')
df['Height'] = df['CL_Z'] - df['Av_Base_Te']

print("Height calculated")
arcpy.AddMessage("Height calculated")

# Calculate asymmetry
df.insert(12, 'Width', '')
df['Width'] = df['Shape_Length'].rolling(window=2).sum()
df.insert(11, 'Asymmetry', '')
df['Asymmetry'] = df['Shape_Length'] / df['Width']

print("Asymmetry calculated")
arcpy.AddMessage("Asymmetry calculated")

# Remove duplicate rows (split transects)
df.drop(df.index[df['Value_1'] == 0], inplace=True)

print("Duplicate rows removed")
arcpy.AddMessage("Duplicate rows removed")

# Calculate cross-sectional area
df.insert(11, 'CS_Ar', '')
df['CS_Ar'] = 0.5 * df['Width'] * df['Height']

print("Cross-sectional area calculated")
arcpy.AddMessage("Cross-sectional area calculated")

# Calculate cross-sectional volume
df.insert(12, 'CS_Vol', '')
df['CS_Vol'] = df['CS_Ar'] * DistanceSplit

print("Cross-sectional volume calculated")
arcpy.AddMessage("Cross-sectional volume calculated")

# Rename ObjectID column to OID
df.rename(columns={'OBJECTID': 'FID'})

# Remove unwanted columns
df = df.drop(['Transect_Z', 'Value_1', 'Shape_Length'], axis=1)

print("Unwanted columns removed")
arcpy.AddMessage("Unwanted columns removed")

# Export to final csv file
df.to_csv(output_folder+"\\"+"Transect_Morphometry.csv")

print("Calculated transect data exported to csv")
arcpy.AddMessage("Calculated transect data exported to csv")

# ----------------------------------------------------------------------------------------
# JOIN CALCULATED MORPH DATA TO TRANSECT SHAPEFILE

# set local variables
inFeatures = Transect_Data
joinTable = output_folder+"\\"+"Transect_Morphometry.csv"
joinField = "x_start"

Transect_MorphometryCopy = arcpy.CopyRows_management(joinTable, output_folder+"\\"+"Transect_MorphometryCopy")

Transect_MorphometryFL = arcpy.MakeFeatureLayer_management(inFeatures, "Transect_MorphometryFL")
arcpy.JoinField_management(Transect_MorphometryFL, joinField, Transect_MorphometryCopy, joinField)
Transect_MorphometryJOIN = arcpy.CopyFeatures_management(Transect_MorphometryFL, output_folder+"\\"+"Transect_Morphometry")

print("Calculated Morphometrics Joined to Transect Shapefile")
arcpy.AddMessage("Calculated Morphometrics Joined to Transect Shapefile")

## Delete unwanted columns ##

# Set local variables
inFeatures = Transect_MorphometryJOIN
Transect_Morphometry = output_folder+"\\"+"Transect_Morphometry.shp"
dropFields = ["OBJECTID", "FIELD1", "X_START_1", "Y_START_1", "X_END_1", "Y_END_1", "OBJECTID_1", "CL_Z_1", "AVG_SLOP_1", "FEATURE__1"]

# Execute DeleteField
arcpy.DeleteField_management(Transect_Morphometry, dropFields)

print("Unwanted fields removed")
arcpy.AddMessage("Unwanted fields removed")
# ----------------------------------------------------------------------------------------
# Delete unwanted files

arcpy.Delete_management(CL_POINTS)
arcpy.Delete_management(Clipped_CLs)
arcpy.Delete_management(FGDB)
arcpy.Delete_management(Clipped_Transects)
arcpy.Delete_management(EXPLODED_JOIN)
arcpy.Delete_management(Exploded_Transects)
arcpy.Delete_management(Final_Join)
arcpy.Delete_management(SPLIT_TRANSECTS)
arcpy.Delete_management(TP_CL_Join)
arcpy.Delete_management(Transect_Data1)
arcpy.Delete_management(TRANSECT_POINTS)
arcpy.Delete_management(Transect_Points_Join)
arcpy.Delete_management(OutputTransect)
arcpy.Delete_management(Transect_Data2)
arcpy.Delete_management(Transect_Data3)
arcpy.Delete_management(CL_JOIN1)
arcpy.Delete_management(Transect_Data)
arcpy.Delete_management(output_folder+"\\"+"transect_morphometrycopy")

import os
os.remove(output_folder+"\\"+"Transect_Data.xls")
import os
os.remove(output_folder+"\\"+"Transect_Data_1.csv")

print("Unnecessary files removed")
arcpy.AddMessage("Unnecessary files removed")

# ----------------------------------------------------------------------------------------
df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='Height', color='red', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Height (m)")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Height.png")

df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='Width', color='turquoise', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Width (m)")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Width.png")

df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='Asymmetry', color='gold', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Asymmetry")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Asymmetry.png")

df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='Avg_Slope', color='blue', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Average Slope (degrees)")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Average_Slope.png")

df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='CS_Ar', color='green', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Cross-Sectional Area ($m^2$)")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Cross_Sectional_Area.png")

df = pd.DataFrame(pd.read_csv(output_folder+"\\"+"Transect_Morphometry.csv"))
fig = df.hist(column='CS_Vol', color='maroon', edgecolor='white')
ax = fig[0]
for x in ax:
    x.set_xlabel("Cross-Sectional Volume ($m^3$)")
    x.set_ylabel("Frequency")
import matplotlib.pyplot as plt
plt.savefig(output_folder+"\\"+"Cross_Sectional_Volume.png")
# -----------------------------------------------------------------------------------------
# Delete unwanted files

import os
os.remove(output_folder+"\\"+"Transect_Morphometry.csv")

print("Unnecessary files removed")
arcpy.AddMessage("Unnecessary files removed")

# -----------------------------------------------------------------------------------------
# Script successfully completed

print("Script completed with no errors")
arcpy.AddMessage("Script completed with no errors")

# -----------------------------------------------------------------------------------------







































