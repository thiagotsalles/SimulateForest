"""
This script builds an ArcGIS toolbox containing 5 tools. They aim to create a
group of points that represent trees of an even-aged forest, with attribute
tables containing DBH, height, and volume values.
"""

# Imports and config
import arcpy
import random
import numpy as np
from scipy.stats import norm, gamma, weibull_min
from arcpy import env
env.overwriteOutput = True  # Allows to overwrite files


class Toolbox(object): # Toolbox containing 5 tools
    def __init__(self):
        self.label = "Simulate forest"
        self.alias = "simulate_forest"
        self.tools = [CreateTreesandDiameters,
                      ProjectDiameters,
                      EstimateHeight,
                      EstimateVolume,
                      dbhClasses]


class CreateTreesandDiameters(object): # Tool number 1
    def __init__(self):
        self.label = "Create trees with DBH"
        self.description = (  # Tool description
            "Creates a mesh of points based on a surface and spacing values. "
            "Also generates DBH values for these points based on a probability distribution.")
        self.canRunInBackground = False

    # Definition of inputs and outputs
    def getParameterInfo(self):
        forest_surface = arcpy.Parameter(
            displayName="Forest surface",
            name="forest_surface",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")
        forest_surface.filter.list = ["Polygon"]

        spacing = arcpy.Parameter(
            displayName="Spacing",
            name="spacing",
            datatype="String",
            parameterType="Required",
            direction="Input")
        spacing.value = "3 x 3"  # Default value as example

        distribution = arcpy.Parameter(
            displayName="Distribution",
            name="distribution",
            datatype="String",
            parameterType="Required",
            direction="Input")
        distribution.filter.list = ["Weibull", "Normal", "Gamma"]  # Available distributions
        distribution.value = "Weibull"  # Default value as example

        age = arcpy.Parameter(
            displayName="Age",
            name="age",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        age.value = 2.0  # Default value as example

        location_prm = arcpy.Parameter(
            displayName="Location parameter",
            name="location_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        location_prm.value = 2.33630  # Default value as example

        scale_prm = arcpy.Parameter(
            displayName="Scale parameter",
            name="scale_prm",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        scale_prm.value = 9.53052  # Default value as example

        form_prm = arcpy.Parameter(
            displayName="Form parameter",
            name="form_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        form_prm.value = 3.50449  # Default value as example

        live_trees = arcpy.Parameter(
            displayName="% of live trees",
            name="live_trees",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        live_trees.filter.type = "Range"
        live_trees.filter.list = [1.0, 100.0]  # Range limit
        live_trees.value = 100.0  # Default value as example

        out_feature = arcpy.Parameter(  # Output
            displayName="Output feature",
            name="out_feature",
            datatype="Feature Class",
            parameterType="Required",
            direction="Output")

        points_ha = arcpy.Parameter(
            displayName=(
                "Number of random points per hectare created to guide DBH allocation"),
            name="points_ha",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            category="Spatial distribution parameters")
        points_ha.value = 130  # Default value

        min_dist = arcpy.Parameter(
            displayName=("Minimum distance between random points"),
            name="min_dist",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            category="Spatial distribution parameters")
        min_dist.value = 25.0  # Default value

        heterogeneity = arcpy.Parameter(
            displayName=("Level of heterogeneity applied to DBH allocation"),
            name="heterogeneity",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            category="Spatial distribution parameters")
        heterogeneity.filter.type = "Range"
        heterogeneity.filter.list = [0.0, 1.0]  # Range limit
        heterogeneity.value = 0.30  # Default value

        params = [forest_surface,
                  spacing,
                  distribution,
                  age,
                  location_prm,
                  scale_prm,
                  form_prm,
                  live_trees,
                  out_feature,
                  points_ha,
                  min_dist,
                  heterogeneity]
        return params

    # License conditions
    def isLicensed(self):
        return True

    # Field validation according to what is inserted as parameter
    def updateParameters(self, parameters):
        distribution = parameters[2]
        location_prm = parameters[4]
        scale_prm = parameters[5]
        form_prm = parameters[6]

        if distribution.value == "Normal":
            location_prm.enabled = True
            scale_prm.enabled = True
            form_prm.value = None
            form_prm.enabled = False
        elif distribution.value == "Weibull":
            location_prm.enabled = True
            scale_prm.enabled = True
            form_prm.enabled = True
        elif distribution.value == "Gamma":
            location_prm.value = None
            location_prm.enabled = False
            scale_prm.enabled = True
            form_prm.enabled = True
        else:
            location_prm.value = None
            location_prm.enabled = False
            scale_prm.value = None
            scale_prm.enabled = False
            form_prm.value = None
            form_prm.enabled = False
        return

    # Warning and error messagens
    def updateMessages(self, parameters):
        spacing = parameters[1]
        distribution = parameters[2]
        location_prm = parameters[4]
        form_prm = parameters[6]

        # Field validation for conditional required values
        if distribution.value == "Weibull":
            if not form_prm.value and (form_prm.value != 0):
                form_prm.setErrorMessage(
                    "Shape parameter is required for selected distribution.")
        elif distribution.value == "Normal":
            if not location_prm.value and (location_prm.value != 0):
                location_prm.setErrorMessage(
                    "Location parameter is required for selected distribution.")
        elif distribution.value == "Gamma":
            if not form_prm.value and (form_prm.value != 0):
                form_prm.setErrorMessage(
                    "Form parameter is required for selected distribution")

        # Checks if spacing is in the right format. It should be a string with
        # three items: a float, an 'x' and another float. Blank spaces make no
        # difference.
        def convert_to_float(x):
            try:
                return float(x)
            except ValueError:
                return x

        spacing_txt = spacing.valueAsText
        if spacing.value:
            if hasattr(spacing_txt, "split"):
                # If the value of spacing can be splitted, it creates a list
                # with each character of the given spacing. The character is
                # converted to float if it's possible.
                spacings = [convert_to_float(k) for k in spacing_txt.split()]
            else:
                spacing.setErrorMessage(
                    'Spacing must be in format: "number x number"')

            condition_1 = len(spacings) != 3
            condition_2 = str(spacings[1]).lower() != "x"
            condition_3 = type(spacings[0]) is not float
            condition_4 = type(spacings[2]) is not float
            if condition_1 or condition_2 or condition_3 or condition_4:
                spacing.setErrorMessage(
                    'Spacing must be in format: "number x number"')
        return

    # Tool execution
    def execute(self, parameters, messages):
        forest_surface = parameters[0].valueAsText
        spacing = parameters[1].valueAsText
        distribution = parameters[2].valueAsText
        age = parameters[3].value
        location_prm = parameters[4].value
        scale_prm = parameters[5].value
        form_prm = parameters[6].value
        live_trees = parameters[7].value
        out_feature = parameters[8].valueAsText
        points_ha = parameters[9].value
        min_dist = parameters[10].value
        heterogeneity = parameters[11].value

        # Values for extent, origin point and y axis.
        desc_surface = arcpy.Describe(forest_surface)
        xtent_list = str(desc_surface.extent).split()[:4]
        xtent = " ".join(k for k in xtent_list)
        origin_point = " ".join(xtent_list[:2])
        y_axis = xtent_list[0] + " " + str(float(xtent_list[1]) + 10)

        # Sets the coordinate system of the main layer as the one for the outputs.
        env.outputCoordinateSystem = desc_surface.SpatialReference

        # Creates a field called StandArea and calculates the area for each
        # feature of the forest_surface layer.
        arcpy.AddField_management(forest_surface, "StandArea", "DOUBLE")
        with arcpy.da.UpdateCursor(forest_surface, ["SHAPE@AREA", "StandArea"]
        ) as cursor:
            for row in cursor:
                row[1] = round(row[0] / 10000.0, 2)
                cursor.updateRow(row)

        # Creates and calculates the field for TotalArea of the stands.
        stands_area_list = arcpy.da.TableToNumPyArray(
            forest_surface, "StandArea", "", True).astype(np.float32).tolist()

        total_area = round(sum(stands_area_list), 2)
        arcpy.AddField_management(forest_surface, "TotalArea", "DOUBLE")
        with arcpy.da.UpdateCursor(forest_surface, ["TotalArea"]) as cursor:
            for row in cursor:
                row[0] = total_area
                cursor.updateRow(row)

        # Assigns the values of spacing to separate variables.
        spacings = [float(x) for x in spacing.split() if x.lower() != 'x']
        spacing_1 = spacings[0]
        spacing_2 = spacings[1]

        # Creates a fishnet with labels (points) for the given spacing and
        # extent. The points are the representation of the trees.
        arcpy.CreateFishnet_management(
            out_feature, origin_point, y_axis, spacing_1, spacing_2,
            "", "", "", "LABELS", xtent)

        # Erases fishnet, leaving only the feature containing the points
        # (labels), wich ends with _label.
        desc_out = arcpy.Describe(out_feature)

        # Creates an intersection feature based on forest_extension and the
        # points of the trees (_label).
        path = desc_out.Path
        catalog_path = desc_out.catalogPath
        basename = desc_out.baseName
        data_type = arcpy.Describe(catalog_path).dataType
        if data_type not in ("FeatureLayer", "FeatureClass"):
            labels = path + "\\" + basename + "_label.shp"
        else:
            labels = path + "\\" + basename + "_label"

        intersect = "in_memory" + "\\" + 'intersect'
        arcpy.Intersect_analysis(
            [forest_surface, labels], intersect, "ALL", "", "POINT")

        # Saves the result of the intersection, that was in the memory, as a
        # file. The file has its fields mapped so that it contains only the
        # area and ID fields.
        forest_surface_id = ('FID_' + desc_surface.baseName)
        fmps1 = arcpy.FieldMappings()
        fmps1.addTable(intersect)
        field_list = [str(forest_surface_id), 'StandArea', 'TotalArea']
        for field in fmps1.fields:
            if field.name not in field_list:
                fmps1.removeFieldMap(fmps1.findFieldMapIndex(field.name))

        out_name = arcpy.Describe(out_feature).file
        arcpy.FeatureClassToFeatureClass_conversion(
            intersect, path, out_name, field_mapping=fmps1)

        # Erases the _label feature and the intersect in memory, leaving only
        # the feature containing the intersect (path + out_name).
        arcpy.Delete_management(labels)
        arcpy.Delete_management(intersect)

        # Creates and calculates the Spacing field.
        arcpy.AddField_management(out_feature, "Spacing", "TEXT")
        with arcpy.da.UpdateCursor(out_feature, ["Spacing"]) as cursor:
            for row in cursor:
                row[0] = spacing
                cursor.updateRow(row)

        # Creates and calculates the field containing the number of trees per
        # hectare at time zero (creation of the trees, without any mortality).
        n_trees_t0 = int(arcpy.GetCount_management(out_feature).getOutput(0))
        N_ha = round(n_trees_t0 / total_area, 0)
        arcpy.AddField_management(out_feature, "N_t0", "FLOAT")
        with arcpy.da.UpdateCursor(out_feature, ["N_t0"]) as cursor:
            for row in cursor:
                row[0] = N_ha
                cursor.updateRow(row)

        # Adds and calculates a field for age.
        arcpy.AddField_management(out_feature, "Age_t1", "FLOAT")
        with arcpy.da.UpdateCursor(out_feature, ["Age_t1"]) as cursor:
            for row in cursor:
                row[0] = age
                cursor.updateRow(row)

        # Adds a field for dbh.
        arcpy.AddField_management(out_feature, "dbh_t1", "FLOAT")

        # Definition of distribution parameters.
        if location_prm == "#" or not location_prm:
            location_prm = 0.0

        location = location_prm
        scale = scale_prm
        form = form_prm

        # Definition of the models.
        def Weibull(location, scl, form, p):
            return weibull_min.ppf(p, form, loc=location, scale=scl)

        def Normal(mean, standard_deviation, p):
            return norm.ppf(p, loc=mean, scale=standard_deviation)

        def Gamma(scl, form, p):
            return gamma.ppf(p, a=form, scale=scl)

        # Calculation of dbh.
        oid_list = [k[0] for k in
                    sorted(arcpy.da.SearchCursor(out_feature, ['OID@']))]

        if distribution == 'Weibull':
            dbh_list = [Weibull(location, scale, form, random.random())
                       for k in oid_list]
        elif distribution == 'Normal':
            dbh_list = [Normal(location, scale, random.random())
                       for k in oid_list]
        elif distribution == 'Gamma':
            dbh_list = [Gamma(scale, form, random.random())
                       for k in oid_list]

        # Base folder adress.
        folder = arcpy.Describe(out_feature).path
        desc = arcpy.Describe(folder)
        while desc.dataType != 'Folder':
            folder = desc.path
            desc = arcpy.Describe(folder)

        # Creation of a base surface with spatial correlation for allocating all DBH.
        # Step 1: Create a point on each vertice of the forest surface extent.
        out_cell_size = (spacing_1 * spacing_2)**(0.5)
        vertices_names = ['1', '2', '3', '4']
        xtent_list_float = [float(k) for k in xtent_list]
        vertices_coords = [(xtent_list_float[0] - 3*out_cell_size,
                            xtent_list_float[1] - 3*out_cell_size),
                           (xtent_list_float[0] - 3*out_cell_size,
                            xtent_list_float[3] + 3*out_cell_size),
                           (xtent_list_float[2] + 3*out_cell_size,
                            xtent_list_float[1] - 3*out_cell_size),
                           (xtent_list_float[2] + 3*out_cell_size,
                            xtent_list_float[3] + 3*out_cell_size)]

        vertices_list = [(k, i[0], i[1]) for k, i in
                        zip(vertices_names, vertices_coords)]

        fields_and_datatype = {'names': ('point', 'coord_x', 'coord_y'),
                               'formats': ('a1', np.dtype(float), np.dtype(float))}

        vertices_array = np.rec.fromrecords(
            vertices_list, dtype=fields_and_datatype)

        out_table = "in_memory" + "\\" + 'out_table'
        arcpy.da.NumPyArrayToTable(vertices_array, out_table)
        arcpy.MakeXYEventLayer_management(
            out_table, 'coord_x', 'coord_y', 'xtent_vertices')

        # Step 2: Creates random points for each 1 ha of the forest surface,
        # maintaining a minimum distance between points.
        out_dissolve = folder + "\\" + 'out_dissolve.shp'
        arcpy.Dissolve_management(forest_surface, out_dissolve)
        n_rand_points = int(total_area * points_ha)
        arcpy.CreateRandomPoints_management(
            folder, 'rand_points', out_dissolve,
            number_of_points_or_field=n_rand_points,
            minimum_allowed_distance=min_dist)

        # Step 3: Merge the extent points and random points into a feature
        # called surface_points.
        surface_points = folder + "\\" + 'surface_points.shp'
        rand_points = folder + "\\" + 'rand_points.shp'
        fmps2 = arcpy.FieldMappings()
        fmps2.addTable('xtent_vertices')
        fmps2.addTable(rand_points)
        for field in fmps2.fields:
            if field.name not in ["dbh_t1"]:
                fmps2.removeFieldMap(fmps2.findFieldMapIndex(field.name))

        arcpy.Merge_management(
            ['xtent_vertices', rand_points], surface_points, fmps2)

        # Step 4: Adds a field for values in the surface_points and erases
        # unecessary features.
        arcpy.AddField_management(surface_points, "value", "DOUBLE")
        arcpy.Delete_management(rand_points)
        arcpy.Delete_management(out_dissolve)
        arcpy.Delete_management(out_table)

        # Step 5: Creates values for the points in surface_points, according
        # to the diameter distribution, and generate the raster surface using
        # Natural Neighbour interpolation. Erases the surface_points.
        if distribution == 'Weibull':
            srf_points_dbh = {
                k[0]: Weibull(location, scale, form, random.random())
                for k in arcpy.da.SearchCursor(surface_points, ['OID@'])}

        elif distribution == 'Normal':
            srf_points_dbh = {
                k[0]: Normal(location, scale, random.random())
                for k in arcpy.da.SearchCursor(surface_points, ['OID@'])}

        elif distribution == 'Gamma':
            srf_points_dbh = {
                k[0]: Gamma(scale, form, random.random())
                for k in arcpy.da.SearchCursor(surface_points, ['OID@'])}

        with arcpy.da.UpdateCursor(surface_points, ['OID@', 'value']
        ) as cursor:
            for row in cursor:
                row[1] = srf_points_dbh[row[0]]
                cursor.updateRow(row)

        nn_raster = folder + "\\" + 'NN_raster'
        outNN = arcpy.sa.NaturalNeighbor(
            surface_points, 'value', out_cell_size)

        outNN.save(nn_raster)
        arcpy.Delete_management(surface_points)

        # Extracts the values from the raster to a table, together with the
        # corresponding OIDs from the point feature of the trees. Erases the raster.
        raster_table = nn_raster + "_table"
        arcpy.ExtractValuesToTable_ga(out_feature, nn_raster, raster_table)
        arcpy.Delete_management(nn_raster)

        # Gives some heterogeneity to the values of the table from the raster.
        if heterogeneity > 0:
            with arcpy.da.UpdateCursor(raster_table, ['VALUE']) as cursor:
                for row in cursor:
                    value = abs(row[0])
                    row[0] = np.random.triangular(
                        value * (1 - heterogeneity), value, value * (1 + heterogeneity))
                    cursor.updateRow(row)

        # Makes a list of the values and OIDs from the raster_table. OIDs are
        # then sorted based on the values. This will serve as base for
        # allocating the dbh.
        table_list = [(k[1], k[0]) for k in
                      arcpy.da.SearchCursor(raster_table, ['SrcID_Feat', 'Value'])]

        table_list.sort()
        arcpy.Delete_management(raster_table)

        # dbh_list is also sorted.
        dbh_list.sort()

        # Now that the list from the raster (table_list) and the list of dbh
        # (dbh_list) are sorted, assigns OIDs from table_list to the values of
        # dbh. This assignment corresponds to the allocation of the dbh
        # according to the previous raster.
        final_dbh_oids = {k[1]: round(v, 2) for k, v in
                          zip(table_list, dbh_list)}

        # Fills the dbh field.
        with arcpy.da.UpdateCursor(out_feature, ['OID@', "dbh_t1"]) as cursor:
            for row in cursor:
                if row[1] in (None, 0):
                    row[1] = final_dbh_oids[row[0]]
                    cursor.updateRow(row)
                else:
                    pass
                    cursor.updateRow(row)

        # Dead trees at the creation of dbh.
        arcpy.AddMessage("Removing dead trees")

        # Takes a random sample of trees, based on the percentage of living trees.
        desc = arcpy.Describe(out_feature)
        data_type = arcpy.Describe(desc.catalogPath).dataType
        n_trees_t1 = int(arcpy.GetCount_management(out_feature).getOutput(0))
        ndead_trees = int(round((n_trees_t1 * (100 - live_trees)) / 100.0, 0))
        random_sample = random.sample(oid_list[1:], ndead_trees)
        random_sample_str = ", ".join(str(k) for k in random_sample)

        # Selects rows in the tree layer, based on the random sample that was
        # previously defined, than sets the dbh of the selected rows to Null or
        # zero depending on the data type.
        if live_trees < 100:
            selection_expression = (
                "" + desc.OIDFieldName + ' IN (' + random_sample_str + ')' + "")

            arcpy.MakeFeatureLayer_management(out_feature, "temp_layer")
            arcpy.SelectLayerByAttribute_management(
                "temp_layer", "NEW_SELECTION", selection_expression)

            if data_type not in ("FeatureLayer", "FeatureClass"):
                with arcpy.da.UpdateCursor("temp_layer", ["dbh_t1"]) as cursor:
                    for row in cursor:
                        row[0] = 0.00
                        cursor.updateRow(row)
            else:
                with arcpy.da.UpdateCursor("temp_layer", ["dbh_t1"]) as cursor:
                    for row in cursor:
                        row[0] = None
                        cursor.updateRow(row)

        # Calculates the number of living trees / Calculates the total area /
        # Calculates the number of trees per hectare / Adds the field for N/ha
        # and calculates it.
        N = len([k[0] for k in
                arcpy.da.SearchCursor(out_feature, ["dbh_t1"])
                if k[0] not in (None, 0)])

        N_ha = round(N / total_area, 0)
        arcpy.AddField_management(out_feature, "N_" + "dbh_t1", "FLOAT")
        with arcpy.da.UpdateCursor(out_feature, ["N_" + "dbh_t1"]) as cursor:
            for row in cursor:
                row[0] = N_ha
                cursor.updateRow(row)
        return


class ProjectDiameters(object):  # Tool number 2
    def __init__(self):
        self.label = "Project DBH and survival"
        self.description = (  # Tool description
            "Projects DBH and tree survival for a future age.")
        self.canRunInBackground = False

    # Definition of inputs and outputs
    def getParameterInfo(self):
        trees = arcpy.Parameter(
            displayName="Layer containing trees",
            name="trees",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")
        trees.filter.list = ["Point"]

        eq_type = arcpy.Parameter(
            displayName="Equation type",
            name="eq_type",
            datatype="String",
            parameterType="Required",
            direction="Input")
        eq_type.filter.type = "ValueList"
        eq_type.filter.list = ["Lundqvist_Korf", "Schumacher", "Richards"]  #  Available equations
        eq_type.value = "Richards"  # Default value as example

        b1_dbh = arcpy.Parameter(
            displayName=u"\u03b2\u2081 parameter",
            name="b1_dbh",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        b1_dbh.value = 0.63327  # Default value as example

        b2_dbh = arcpy.Parameter(
            displayName=u"\u03b2\u2082 parameter",
            name="b2_dbh",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b2_dbh.value = 0.17115  # Default value as example

        b3_dbh = arcpy.Parameter(
            displayName=u"\u03b2\u2083 parameter",
            name="b3_dbh",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b3_dbh.value = 0.08032  # Default value as example

        initial_age_field = arcpy.Parameter(
            displayName="Initial age",
            name="initial_age_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        initial_age_field.parameterDependencies = [trees.name]

        dbh_init_field = arcpy.Parameter(
            displayName="DBH at initial age",
            name="dbh_init_field",
            datatype="Field",
            parameterType="Required",
            direction="Input",)
        dbh_init_field.parameterDependencies = [trees.name]

        years = arcpy.Parameter(
            displayName="Years until final age",
            name="years",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        interval = arcpy.Parameter(
            displayName="Interval between projections",
            name="interval",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        interval.value = 1  # Default value

        time_first_proj = arcpy.Parameter(
            displayName="Time at first projection",
            name="time_first_proj",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        time_first_proj.filter.type = "Range"
        time_first_proj.filter.list = [2.00, float('inf')]

        ntrees_init_field = arcpy.Parameter(
            displayName="Trees per hectare at initial age",
            name="ntrees_init_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        ntrees_init_field.parameterDependencies = [trees.name]

        surv_type = arcpy.Parameter(
            displayName="Survival type",
            name="surv_type",
            datatype="String",
            parameterType="Required",
            direction="Input",
            category="Survival")
        surv_type.type = "ValueList"
        surv_type.filter.list = ["Percentage", "Pienaar and Shiver model"]
        surv_type.value = "Pienaar and Shiver model"  # Default value

        surv_percent = arcpy.Parameter(
            displayName="Survival percent at final age",
            name="surv_percent",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            category="Survival")
        surv_percent.value = 100  # Default value
        surv_percent.filter.type = "Range"
        surv_percent.filter.list = [1.00, 100.00]

        total_area = arcpy.Parameter(
            displayName="Total area",
            name="total_area",
            datatype=["Field", "Double"],
            parameterType="Optional",
            direction="Input",
            category="Survival")
        total_area.parameterDependencies = [trees.name]

        b1_n = arcpy.Parameter(
            displayName=u"\u03b2\u2081 parameter",
            name="b1_n",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            category="Survival")
        b1_n.value = 0.00382  # Default value as example

        b2_n = arcpy.Parameter(
            displayName=u"\u03b2\u2082 parameter",
            name="b2_n",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            category="Survival")
        b2_n.value = 1.70542  # Default value as example

        het_index = arcpy.Parameter(
            displayName="Heterogeneity index",
            name="het_index",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            category="Heterogeneity")
        het_index.value = 0.5  # Default value
        het_index.filter.type = "Range"
        het_index.filter.list = [0.00, 1.00]

        params = [trees,
                  eq_type,
                  b1_dbh,
                  b2_dbh,
                  b3_dbh,
                  initial_age_field,
                  dbh_init_field,
                  years,
                  interval,
                  time_first_proj,
                  ntrees_init_field,
                  surv_type,
                  surv_percent,
                  total_area,
                  b1_n,
                  b2_n,
                  het_index]
        return params

    # License conditions
    def isLicensed(self):
        return True

    # Field validation according to what is inserted as parameter
    def updateParameters(self, parameters):
        eq_type = parameters[1]
        b1_dbh = parameters[2]
        b2_dbh = parameters[3]
        b3_dbh = parameters[4]
        dbh_init_field = parameters[6]
        time_first_proj = parameters[9]
        surv_type = parameters[11]
        surv_percent = parameters[12]
        total_area = parameters[13]
        b1_n = parameters[14]
        b2_n = parameters[15]

        # Converts time_first_proj to int.
        if time_first_proj.value:
            time_first_proj.value = int(time_first_proj.value)

        # Validation according to what is inserted as eq_type.
        if eq_type.value == "Schumacher":
            b1_dbh.enabled = True
            b2_dbh.value = None
            b2_dbh.enabled = False
            b3_dbh.value = None
            b3_dbh.enabled = False
        elif eq_type.value == "Lundqvist_Korf":
            b1_dbh.enabled = True
            b2_dbh.enabled = True
            b3_dbh.value = None
            b3_dbh.enabled = False

        elif eq_type.value == "Richards":
            b1_dbh.enabled = True
            b2_dbh.enabled = True
            b3_dbh.enabled = True
        else:
            b1_dbh.value = None
            b1_dbh.enabled = False
            b2_dbh.value = None
            b2_dbh.enabled = False
            b3_dbh.value = None
            b3_dbh.enabled = False

        # Checks if the two last characters of dbh_init_field name are int,
        # sums it + 1 and sets them as default value for time_first_proj.
        my_list = []
        if dbh_init_field.value:
            for x in list(dbh_init_field.valueAsText[-2:]):
                try:
                    my_list.append(int(x))
                except:
                    pass

        cond_1 = dbh_init_field.value
        cond_2 = not time_first_proj.value
        cond_3 = len(my_list) > 0
        if cond_3:
            cond_4 = int("".join(str(x) for x in my_list))
        cond_5 = time_first_proj.value
        if cond_1 and cond_2 and cond_3:
            time_first_proj.value = cond_4 + 1

        # If time_first_proj <= time at dbh_init_field, changes time_first_proj
        # to a greater value than time at dbh_init_field.
        if cond_1 and cond_3 and cond_5:
            if time_first_proj.value <= cond_4:
                time_first_proj.value = cond_4 + 1

        # Field validation for survival
        if surv_type.value == "Percentage":
            surv_percent.enabled = True
            b1_n.value = None
            b1_n.enabled = False
            b2_n.value = None
            b2_n.enabled = False
            if surv_percent.value == 100:
                total_area.enabled = False
                total_area.value = None
            else:
                total_area.enabled = True
        elif surv_type.value == "Pienaar and Shiver model":
            surv_percent.value = None
            surv_percent.enabled = False
            b1_n.enabled = True
            b2_n.enabled = True
        else:
            surv_percent.value = None
            surv_percent.enabled = False
            b1_n.value = None
            b1_n.enabled = False
            b2_n.value = None
            b2_n.enabled = False
        return

    # Warning and error messages
    def updateMessages(self, parameters):
        eq_type = parameters[1]
        b2_dbh = parameters[3]
        b3_dbh = parameters[4]
        surv_type = parameters[11]
        surv_percent = parameters[12]
        total_area = parameters[13]
        b1_n = parameters[14]
        b2_n = parameters[15]

        # Error messages for parameters from DBH equations.
        if eq_type.value:
            if eq_type.value in ("Lundqvist_Korf", "Richards"):
                if not b2_dbh.value and (b2_dbh.value != 0):
                    b2_dbh.setErrorMessage(
                        u"\u03b2\u2082 parameter is required for selected model.")
            if eq_type.value == "Richards":
                if not b3_dbh.value and (b3_dbh.value != 0):
                    b3_dbh.setErrorMessage(
                        u"\u03b2\u2083 parameter is required for selected model.")

        # Error messages for survival calculation inputs.
        if surv_type.value != "Percentage":
            if not total_area.value:
                total_area.setErrorMessage(
                    "Total area is required to calculate trees per hectare.")
            if not b1_n.value and (b1_n.value != 0):
                b1_n.setErrorMessage(
                    u"\u03b2\u2081 parameter is required for selected model.")
            if not b2_n.value and (b2_n.value != 0):
                b2_n.setErrorMessage(
                    u"\u03b2\u2083 parameter is required for selected model.")
        else:
            if not surv_percent.value:
                surv_percent.setErrorMessage(
                    "Survival percent is required to calculate trees per hectare.")
            if surv_percent.value < 100:
                if not total_area.value:
                    total_area.setErrorMessage(
                        "Total area is required to calculate trees per hectare.")
        return

    # Tool execution
    def execute(self, parameters, messages):
        trees = parameters[0].valueAsText
        eq_type = parameters[1].valueAsText
        b1_dbh = parameters[2].value
        b2_dbh = parameters[3].value
        b3_dbh = parameters[4].value
        initial_age_field = parameters[5].valueAsText
        dbh_init_field = parameters[6].valueAsText
        years = parameters[7].value
        interval = parameters[8].value
        time_first_proj = parameters[9].value
        ntrees_init_field = parameters[10].valueAsText
        surv_type = parameters[11].valueAsText
        surv_percent = parameters[12].value
        total_area = parameters[13].valueAsText
        b1_n = parameters[14].value
        b2_n = parameters[15].value
        het_index = parameters[16].value

        arcpy.AddMessage("Preparing data")

        # Definition of the dbh models.
        def Lundqvist_Korf(dbh1, a1, a2, b1, b2):
            return dbh1 * np.exp(b1 / a1 ** b2 - b1 / a2 ** b2)

        def Schumacher(dbh1, a1, a2, b1):
            return dbh1 * np.exp(b1 / a1 - b1 / a2)

        def Richards(dbh1, a1, a2, b1, b2, b3):
            return ((dbh1 * (1 - b1 * np.exp(-b2 * a2)) ** (1 / (1 - b3))) /
                    (1 - b1 * np.exp(-b2 * a1)) ** (1 / (1 - b3)))

        # Definition of survival model
        def Pienaar_Shiver(n1, a1, a2, b1, b2):
            return n1 * np.exp(-b1 * (a2 ** b2 - a1 ** b2))

        # At this point, some lists ans dictionaries have to be created to
        # guide the calculations in the 'for' loop that comes in sequence.
        # The loop is used because the projections may be performed more than once.

        # The type of the data.
        desc = arcpy.Describe(trees)
        data_type = arcpy.Describe(desc.catalogPath).dataType

        # The dictionary of the OID's and ages at time zero of the projections.
        init_ages = {k[0]: k[1] for k in
                     arcpy.da.SearchCursor(trees, ["OID@", initial_age_field])}

        # The dictionary containing the age for each OID at each time in the
        # projections (age_dic).
        increments_list = [0]
        increment = interval
        while increment <= years:
            increments_list.append(increment)
            increment += interval
        increments_array = np.array(increments_list)
        age_dic = {k: increments_array + init_ages[k] for k in init_ages}

        # The list with each time of the projections.
        time_list = [k + int(time_first_proj) for k in
                     range(len(increments_list) - 1)]

        # The list with the names of the dbh fields for each time in the projections.
        dbh_list = ["dbh_t" + str(k) for k in time_list]
        dbh_list.insert(0, dbh_init_field)

        # The list with the names of the 'trees per hectare' fields for each
        # time in the projections.
        n_list = ["N_" + k for k in dbh_list[1:]]
        n_list.insert(0, ntrees_init_field)

        # Gets the value for the number os trees at time zero (n0).
        with arcpy.da.SearchCursor(trees, [ntrees_init_field]) as cursor:
            for row in cursor:
                n0 = row[0]
                break

        # Gets the value for the age time zero (a0).
        a0_list = [k[0] for k in
                   arcpy.da.SearchCursor(trees, [initial_age_field])
                   if k[0] not in (None, 0)]

        a0 = sum(a0_list) / len(a0_list)

        # Gets the value for the age at final time (a_final).
        a_final_list = [k[0] + increments_list[-1] for k in
                        arcpy.da.SearchCursor(trees, [initial_age_field])
                        if k[0] not in (None, 0)]

        a_final = sum(a_final_list) / len(a_final_list)

        # Fits a logarithmic curve based on the number of trees at initial
        # age an the survival rate at final age.
        cond1 = surv_type == 'Percentage'
        cond2 = surv_percent < 100
        if cond1 and cond2:
            n_final = surv_percent / 100 * n0
            xy_list = [(np.log(a0), n0), (np.log(a_final), n_final)]

            def Beta1(xy_list):
                sum_xy = sum([k[0] * k[1] for k in xy_list])
                sum_x = sum([k[0] for k in xy_list])
                sum_y = sum([k[1] for k in xy_list])
                sum_x2 = sum([k[0] ** 2 for k in xy_list])
                n = len(xy_list)
                b1 = (sum_xy - sum_x * sum_y / n) / (sum_x2 - sum_x ** 2 / n)
                return b1

            def Beta0(xy_list, b1):
                average_x = sum([k[0] for k in xy_list]) / len(xy_list)
                average_y = sum([k[1] for k in xy_list]) / len(xy_list)
                b0 = average_y - b1 * average_x
                return b0

            b1_n = Beta1(xy_list)
            b0_n = Beta0(xy_list, b1_n)

        # Gets the value of the total area.
        if total_area:
            try:
                total_area = float(total_area)
            except ValueError:
                with arcpy.da.SearchCursor(trees, [total_area]) as cursor:
                    for row in cursor:
                        total_area = row[0]
                        break
        else:
            total_area = 0.0

        # The 'for' loop that executes the projections. It runs once for each
        # time given in the time_list.
        for i in range(len(time_list)):
            # Initial dbh field
            i_dbh_f = dbh_list[i]
            # Final dbh field
            dbh_field = dbh_list[i + 1]
            time = time_list[i]
            age_field = "Age_t" + str(time)
            # Trees per hectare field
            n1_field = n_list[i]
            n2_field = n_list[i + 1]

            # Calculations of projected dbh. It involves the creation of a
            # dictionary with values of OID, initial age, final age and initial
            # dbh and uses it combined with the models to estimate the final
            # (projected) dbh.
            arcpy.AddMessage("Calculating DBH for " + dbh_field)

            # Dictionary containing the values of OID, initial age, final age
            # and initial dbh.
            age_dbh = {k[0]: (age_dic[k[0]][i], age_dic[k[0]][i + 1], k[1]) for
                       k in arcpy.da.SearchCursor(trees, ["OID@", i_dbh_f])}

            # Calculation of the growth rate for each tree (dbh2 / dbh1).
            if eq_type == "Lundqvist_Korf":
                growth = {k: (age_dbh[k][2],
                              Lundqvist_Korf(age_dbh[k][2], age_dbh[k][0],
                                             age_dbh[k][1], b1_dbh, b2_dbh) / age_dbh[k][2])
                          if age_dbh[k][2] not in (None, 0)
                          else (age_dbh[k][2], age_dbh[k][2])
                          for k in age_dbh}
            elif eq_type == "Schumacher":
                growth = {k: (age_dbh[k][2],
                              Schumacher(age_dbh[k][2], age_dbh[k][0],
                                         age_dbh[k][1], b1_dbh) / age_dbh[k][2])
                          if age_dbh[k][2] not in (None, 0)
                          else (age_dbh[k][2], age_dbh[k][2])
                          for k in age_dbh}
            elif eq_type == "Richards":
                growth = {k: (age_dbh[k][2],
                              Richards(age_dbh[k][2], age_dbh[k][0],
                                       age_dbh[k][1], b1_dbh, b2_dbh, b3_dbh) / age_dbh[k][2])
                          if age_dbh[k][2] not in (None, 0)
                          else (age_dbh[k][2], age_dbh[k][2])
                          for k in age_dbh}

            # Definition of a formula to apply of a disturbance in the growth
            # rate of each tree.
            def Disturbance(growth_rate, het_index):
                # Recalculates the growth rate. The new growth rate will be
                # a random value from a normal distribution with the original
                # rate as the mean. The standard deviation is defined so that
                # 99.7% of the possible values of the distribution are between
                # 1.0 and 2 * growt_rate - 1. The values of the randomized
                # probability inside the normal distribution are limited so
                # that 100% of the results are within the condition above.
                deviation = (growth_rate - 1) / 3.0
                p_min = norm.cdf(1.0, loc=growth_rate, scale=deviation)
                p_min = (0.5 - p_min) * (1 - het_index) + p_min
                p_max = 1 - p_min
                new_rate = norm.ppf(random.uniform(p_min, p_max),
                                    loc=growth_rate, scale=deviation)
                return new_rate

            # Calculation of the projected dbh (prj_dbh). The 'if' statement is
            # there to maintain the values of dbh of the dead trees as zero or None.
            prj_dbh = {k: round(Disturbance(growth[k][1], het_index) * growth[k][0], 2)
                       if growth[k][0] not in (None, 0)
                       else growth[k][0] for k in growth}

            # Adds and calculates a field for age.
            arcpy.AddField_management(trees, age_field, "FLOAT")

            with arcpy.da.UpdateCursor(trees, ["OID@", age_field]) as cursor:
                for row in cursor:
                    row[1] = age_dic[row[0]][i + 1]
                    cursor.updateRow(row)

            # Adds and fills the field for dbh.
            arcpy.AddField_management(trees, dbh_field, "FLOAT")
            with arcpy.da.UpdateCursor(trees, ["OID@", dbh_field]) as cursor:
                for row in cursor:
                    row[1] = prj_dbh[row[0]]
                    cursor.updateRow(row)

            # Projection survival. When using the Pienaar and Shiver model, it
            # takes the average age of the trees (a1 and a2) as input.
            arcpy.AddMessage("Projecting survival")

            # Values for average initial age (a1) and average final age (a2).
            initial_age_list = [age_dbh[k][0] for k in age_dbh]
            final_age_list = [age_dbh[k][1] for k in age_dbh]
            a1 = sum(initial_age_list) / len(initial_age_list)
            a2 = sum(final_age_list) / len(final_age_list)

            # Calculation of the number of trees per hectare at initial age
            # (n1), trees per hectare at final age (n2) and the total number of
            # dead trees at final age (ndead_trees).
            if surv_type == "Percentage":
                with arcpy.da.SearchCursor(trees, [n1_field]) as cursor:
                    for row in cursor:
                        n1 = row[0]
                        break

                if surv_percent < 100:
                    n2 = b0_n + b1_n * np.log(a2)
                else:
                    n2 = n1
            elif surv_type == "Pienaar and Shiver model":
                n1 = Pienaar_Shiver(n0, a0, a1, b1_n, b2_n)
                n2 = Pienaar_Shiver(n1, a1, a2, b1_n, b2_n)

            ndead_trees = int(round((n1 - n2) * total_area, 2))

            # Takes a random sample of OIDs from living trees, based on the
            # number of dead trees (ndead_trees).
            living_trees = [k[0] for k in sorted(
                            arcpy.da.SearchCursor(trees, ['OID@', dbh_field]))
                            if k[1] not in (None, 0)]
            random_sample = random.sample(living_trees[1:], ndead_trees)
            random_sample_str = ", ".join(str(k) for k in random_sample)

            # Selects rows in the tree layer, based on the random sample that
            # was previously defined, than sets the dbh of the selected rows to
            # Null orzero depending on the data type.
            if random_sample_str is not "":
                selection_expression = (
                    "" + desc.OIDFieldName + ' IN (' + random_sample_str + ')' + "")

                arcpy.MakeFeatureLayer_management(trees, "temp_layer")
                arcpy.SelectLayerByAttribute_management(
                    "temp_layer", "NEW_SELECTION", selection_expression)

                if data_type not in ("FeatureLayer", "FeatureClass"):
                    with arcpy.da.UpdateCursor("temp_layer", [dbh_field]) as cursor:
                        for row in cursor:
                            row[0] = 0.00
                            cursor.updateRow(row)
                else:
                    with arcpy.da.UpdateCursor("temp_layer", [dbh_field]) as cursor:
                        for row in cursor:
                            row[0] = None
                            cursor.updateRow(row)

            arcpy.AddField_management(trees, n2_field, "FLOAT")
            with arcpy.da.UpdateCursor(trees, [n2_field]) as cursor:
                for row in cursor:
                    row[0] = round(n2, 0)
                    cursor.updateRow(row)
        return


class EstimateHeight(object):
    def __init__(self):
        self.label = "Estimar altura"
        self.description = (u"Estima a altura das árvores baseada no dap e "
                            u"idade.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        trees = arcpy.Parameter(
            displayName=u"Layer contendo as árvores",
            name="trees",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")
        trees.filter.list = ["Point"]

        eq_type = arcpy.Parameter(
            displayName="Tipo de modelo",
            name="eq_type",
            datatype="String",
            parameterType="Required",
            direction="Input")
        eq_type.filter.type = "ValueList"
        eq_type.filter.list = ["Lundqvist_Korf", "Schumacher",
                               "Gompertz", "Weibull"]
        eq_type.value = "Weibull"

        dependent_variable = arcpy.Parameter(
            displayName=u"Tipo de variável dependente",
            name="dependent_variable",
            datatype="String",
            parameterType="Required",
            direction="Input")
        dependent_variable.filter.type = "ValueList"
        dependent_variable.filter.list = ["dap", "dap * Idade"]
        dependent_variable.value = "dap * Idade"

        b0_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2080",
            name="b0_prm",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        b0_prm.value = 27.79485

        b1_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2081",
            name="b1_prm",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        b1_prm.value = 20.78744

        b2_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2082",
            name="b2_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b2_prm.value = 0.02690

        b3_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2083",
            name="b3_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b3_prm.value = 0.90071

        dbh_fields = arcpy.Parameter(
            displayName="dap",
            name="dbh_fields",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        dbh_fields.parameterDependencies = [trees.name]

        age_fields = arcpy.Parameter(
            displayName="Idade",
            name="age_fields",
            datatype="Field",
            parameterType="Optional",
            direction="Input",
            multiValue=True)
        age_fields.parameterDependencies = [trees.name]

        het_index = arcpy.Parameter(
            displayName=u"\xcdndice de heterogeneidade",
            name="var_index",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            category=u"Heterogeneidade")
        het_index.value = 0.20
        het_index.filter.type = "Range"
        het_index.filter.list = [0.00, 0.99]

        params = [trees,
                  eq_type,
                  dependent_variable,
                  b0_prm,
                  b1_prm,
                  b2_prm,
                  b3_prm,
                  dbh_fields,
                  age_fields,
                  het_index]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        eq_type = parameters[1]
        dependent_variable = parameters[2]
        b2_prm = parameters[5]
        b3_prm = parameters[6]
        age_fields = parameters[8]

        # Field validation
        if dependent_variable.value == "dap * Idade":
            age_fields.enabled = True
        elif dependent_variable.value == "dap":
            age_fields.enabled = False

        if eq_type.value == "Weibull":
            b2_prm.enabled = True
            b3_prm.enabled = True
        elif eq_type.value in ("Lundqvist_Korf", "Gompertz"):
            b2_prm.enabled = True
            b3_prm.value = None
            b3_prm.enabled = False
        elif eq_type.value == "Schumacher":
            b2_prm.value = None
            b2_prm.enabled = False
            b3_prm.value = None
            b3_prm.enabled = False

    def updateMessages(self, parameters):
        eq_type = parameters[1]
        dependent_variable = parameters[2]
        b2_prm = parameters[5]
        b3_prm = parameters[6]
        dbh_fields = parameters[7]
        age_fields = parameters[8]

        if dbh_fields.value:
            dbh_fields_list = dbh_fields.valueAsText.split(";")
        if age_fields.value:
            age_fields_list = age_fields.valueAsText.split(";")

        if dependent_variable.value == "dap * Idade":
            if not age_fields.value:
                age_fields.setErrorMessage(
                    u"Campo(s) de idade são necessários quando a "
                    u"variável dependente \xea 'dap * Idade'")
            if dbh_fields.value and age_fields.value:
                if len(dbh_fields_list) != len(age_fields_list):
                    dbh_fields.setErrorMessage(
                        u"Núamero de campos de dap e Idade deve ser o mesmo")
                    age_fields.setErrorMessage(
                        u"Núamero de campos de dap e Idade deve ser o mesmo")

        if eq_type.value:
            if eq_type.value is "Weibull" or "Lundqvist_Korf" or "Gompertz":
                if not b2_prm.value and (b2_prm.value != 0):
                    b2_prm.setErrorMessage(
                        u"Parâmetro \xe9 necessário para uso do modelo "
                        "selecionado.")
            if eq_type.value == "Weibull":
                if not b3_prm.value and (b3_prm.value != 0):
                    b3_prm.setErrorMessage(
                        u"Parâmetro \xe9 necessário para uso do modelo "
                        "selecionado.")
        return

    def execute(self, parameters, messages):
        trees = parameters[0].valueAsText
        eq_type = parameters[1].valueAsText
        dependent_variable = parameters[2].valueAsText
        b0_prm = parameters[3].value
        b1_prm = parameters[4].value
        b2_prm = parameters[5].value
        b3_prm = parameters[6].value
        dbh_fields = parameters[7].valueAsText
        age_fields = parameters[8].valueAsText
        het_index = parameters[9].value

        dbh_fields_list = dbh_fields.split(";")

        # Creates a dummy value for age_fields, so the calculations of height
        # can work when dependent variable is only dbh.
        if age_fields == "#" or not age_fields:
            age_fields = len(dbh_fields_list) * "OID@;"
            age_fields = age_fields[:-1]
        age_fields_list = age_fields.split(";")

        # Calculation of a 'seed' for each tree. The seed is an input value to
        # the calculation of pseudo-random numbers, used when establishing a
        # disturbance on the estimates of height. It could be any number, as
        # long as it is a unique value for each tree.
        seeds = {k[0]: k[1][0] + k[1][1] for k in
                 arcpy.da.SearchCursor(trees, ['OID@', 'SHAPE@XY'])}

        # Definition of the models. The models are conditioned by the dependent
        # variable being dbh or dbh*Age.
        def Lundqvist_Korf(b0, b1, b2, dbh, age):
            if dependent_variable == "dap * Idade":
                return b0*np.exp(-b1/(dbh*age)**b2)
            else:
                return b0*np.exp(-b1/dbh**b2)

        def Schumacher(b0, b1, dbh, age):
            if dependent_variable == "dap * Idade":
                return b0*np.exp(-b1/(dbh*age))
            else:
                return b0*np.exp(-b1/dbh)

        def Gompertz(b0, b1, b2, dbh, age):
            if dependent_variable == "dap * Idade":
                return b0*np.exp(-b1*np.exp(-b2*(dbh*age)))
            else:
                return b0*np.exp(-b1*np.exp(-b2*dbh))

        def Weibull(b0, b1, b2, b3, dbh, age):
            if dependent_variable == "dap * Idade":
                return b0-b1*np.exp(-b2*(dbh*age)**b3)
            else:
                return b0-b1*np.exp(-b2*dbh**b3)

        # Definition of a formula to apply of a disturbance on the height of
        # each tree.
        def Disturbance(seed, height, het_index):
            # Recalculates the height of the tree (new_height). The new_height
            # will be a value from a normal distribution with the original
            # height as the mean. The value of std is defined so that 99.7% of
            # the possible values of the distribution are between the minimum
            # and maximum desired new_height. The seed can be any number and it
            # is used as an input to generate a pseudo-random number 'p'. The
            # pseudo-random value generated for a particular seed is always the
            # same. According to the 'seeds' calculation performed previously,
            # each value of seed used here is unique for each tree. Therefore,
            # the pseudo-random number generated for each seed is always the
            # same. 'p' is the probability used to get the new value of height
            # (new_height) in the normal distribution and is limited so the
            # minimum possible new_height is height * (1 - het_index).
            ht_min = height * (1 - het_index)
            std = (height - ht_min) / 3.0
            p_min = norm.cdf(ht_min, loc=height, scale=std)
            p_max = 1 - p_min
            random.seed(seed)
            p = random.uniform(p_min, p_max)
            new_height = norm.ppf(p, loc=height, scale=std)
            return new_height

        # The 'for' loop that estimates the heights. It runs once for each
        # item in dbh_fields_list.
        for dbhfield, agefield in zip(dbh_fields_list, age_fields_list):
            # Dictionary of OIDs, age and dbh.
            dbh_dic = {k[0]: (k[1], k[2]) for k in
                       arcpy.da.SearchCursor(trees,
                                             ["OID@", agefield, dbhfield])}

            arcpy.AddMessage("Estimando alturas de " + dbhfield)

            # Dictionary of OIDs and estimated heights.
            if eq_type == "Lundqvist_Korf":
                heights = {k: Lundqvist_Korf(b0_prm, b1_prm, b2_prm,
                                             dbh_dic[k][1], dbh_dic[k][0])
                           if dbh_dic[k][1] not in (None, 0)
                           else dbh_dic[k][1] for k in dbh_dic}
            elif eq_type == "Schumacher":
                heights = {k: Schumacher(b0_prm, b1_prm,
                                         dbh_dic[k][1], dbh_dic[k][0])
                           if dbh_dic[k][1] not in (None, 0)
                           else dbh_dic[k][1] for k in dbh_dic}
            elif eq_type == "Gompertz":
                heights = {k: Gompertz(b0_prm, b1_prm, b2_prm,
                                       dbh_dic[k][1], dbh_dic[k][0])
                           if dbh_dic[k][1] not in (None, 0)
                           else dbh_dic[k][1] for k in dbh_dic}
            elif eq_type == "Weibull":
                heights = {k: Weibull(b0_prm, b1_prm, b2_prm, b3_prm,
                                      dbh_dic[k][1], dbh_dic[k][0])
                           if dbh_dic[k][1] not in (None, 0)
                           else dbh_dic[k][1] for k in dbh_dic}

            # Estimation of heights. The 'if' statement is there to maintain
            # the values of height of the deadtrees as zero or None.
            new_heights = {k: round(Disturbance(seeds[k], heights[k],
                                                het_index), 2)
                           if heights[k] not in (None, 0)
                           else heights[k] for k in heights}

            # Adds field for height
            height_field = "Ht_" + dbhfield
            arcpy.AddField_management(trees, height_field, "FLOAT")

            # Fills field of height
            with arcpy.da.UpdateCursor(trees, ['OID@',
                                               height_field]) as cursor:
                for row in cursor:
                    row[1] = new_heights[row[0]]
                    cursor.updateRow(row)
        return


class EstimateVolume(object):
    def __init__(self):
        self.label = "Estimar volume"
        self.description = (u"Estima o volume das árvores baseado no dap e "
                            "altura.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        trees = arcpy.Parameter(
            displayName=u"Layer contendo as árvores",
            name="trees",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")
        trees.filter.list = ["Point"]

        eq_type = arcpy.Parameter(
            displayName="Tipo de estimativa",
            name="eq_type",
            datatype="String",
            parameterType="Required",
            direction="Input")
        eq_type.filter.type = "ValueList"
        eq_type.filter.list = ["Modelo Schumacher & Hall", "Fator de forma"]
        eq_type.value = "Modelo Schumacher & Hall"

        b0_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2080",
            name="b0_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b0_prm.value = -10.44771

        b1_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2081",
            name="b1_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b1_prm.value = 1.86852

        b2_prm = arcpy.Parameter(
            displayName=u"Parâmetro \u03b2\u2082",
            name="b2_prm",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        b2_prm.value = 1.17904

        ff = arcpy.Parameter(
            displayName="Fator de forma",
            name="ff",
            datatype="Double",
            parameterType="Optional",
            direction="Input")

        dbh_fields = arcpy.Parameter(
            displayName="dap",
            name="dbh_fields",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        dbh_fields.parameterDependencies = [trees.name]

        h_fields = arcpy.Parameter(
            displayName="Altura",
            name="h_fields",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        h_fields.parameterDependencies = [trees.name]

        params = [trees,
                  eq_type,
                  b0_prm,
                  b1_prm,
                  b2_prm,
                  ff,
                  dbh_fields,
                  h_fields]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        eq_type = parameters[1]
        b0_prm = parameters[2]
        b1_prm = parameters[3]
        b2_prm = parameters[4]
        ff = parameters[5]

        # Form validation
        if eq_type.value == "Modelo Schumacher & Hall":
            ff.value = None
            ff.enabled = False
            b0_prm.enabled = True
            b1_prm.enabled = True
            b2_prm.enabled = True
        else:
            ff.enabled = True
            b0_prm.value = None
            b0_prm.enabled = False
            b1_prm.value = None
            b1_prm.enabled = False
            b2_prm.value = None
            b2_prm.enabled = False
        return

    def updateMessages(self, parameters):
        eq_type = parameters[1]
        b0_prm = parameters[2]
        b1_prm = parameters[3]
        b2_prm = parameters[4]
        ff = parameters[5]
        dbh_fields = parameters[6]
        h_fields = parameters[7]

        if dbh_fields.value:
            dbh_fields_list = dbh_fields.valueAsText.split(";")
        if h_fields.value:
            h_fields_list = h_fields.valueAsText.split(";")

        if dbh_fields.value and h_fields.value:
            if len(dbh_fields_list) != len(h_fields_list):
                dbh_fields.setErrorMessage(
                    u"Núamero de campos de dap e altura deve ser o mesmo")
                h_fields.setErrorMessage(
                    u"Núamero de campos de dap e altura deve ser o mesmo")

        if eq_type.value == "Modelo Schumacher & Hall":
            if not b0_prm.value and (b0_prm.value != 0):
                b0_prm.setErrorMessage(
                    u"Parâmetro \xe9 necessário para uso do modelo "
                    "selecionado.")
            if not b1_prm.value and (b1_prm.value != 0):
                b1_prm.setErrorMessage(
                    u"Parâmetro \xe9 necessário para uso do modelo "
                    "selecionado.")
            if not b2_prm.value and (b2_prm.value != 0):
                b2_prm.setErrorMessage(
                    u"Parâmetro \xe9 necessário para uso do modelo "
                    "selecionado.")
        else:
            if not ff.value:
                ff.setErrorMessage(
                    u"Fator de forma \xe9 necessário para cálculo do "
                    "volume.")
        return

    def execute(self, parameters, messages):
        trees = parameters[0].valueAsText
        eq_type = parameters[1].valueAsText
        b0_prm = parameters[2].value
        b1_prm = parameters[3].value
        b2_prm = parameters[4].value
        ff = parameters[5].value
        dbh_fields = parameters[6].valueAsText
        h_fields = parameters[7].valueAsText

        dbh_fields_list = dbh_fields.split(";")
        h_fields_list = h_fields.split(";")

        # Definition of the models.
        def Schumacher_Hall(b0, b1, b2, dbh, height):
            return np.exp(b0 + b1 * np.log(dbh) + b2 * np.log(height))

        def Form_factor(ff, dbh, height):
            return (np.pi * dbh ** 2 / 40000) * height * ff

        # Adds fields for volume and calculates it.
        for dbhfield, hfield in zip(dbh_fields_list, h_fields_list):
            volume_field = "V_" + dbhfield
            arcpy.AddMessage(u"Calculando volumes para " + volume_field)
            arcpy.AddField_management(trees, volume_field, "FLOAT")
            if eq_type == "Modelo Schumacher & Hall":
                volumes = {k[0]: round(Schumacher_Hall(b0_prm, b1_prm, b2_prm,
                                                       k[1], k[2]), 4)
                           if k[1] not in (None, 0) else k[1] for k in
                           arcpy.da.SearchCursor(trees,
                                                 ['OID@', dbhfield, hfield])}
            elif eq_type == "Fator de forma":
                volumes = {k[0]: round(Form_factor(ff, k[1], k[2]), 4)
                           if k[1] not in (None, 0) else k[1] for k in
                           arcpy.da.SearchCursor(trees,
                                                 ['OID@', dbhfield, hfield])}

            with arcpy.da.UpdateCursor(trees,
                                       ['OID@', volume_field]) as cursor:
                for row in cursor:
                    row[1] = volumes[row[0]]
                    cursor.updateRow(row)
        return


class dbhClasses(object):
    def __init__(self):
        self.label = "Classificar dap"
        self.description = (u"Gera o centro de classe de diâmetro para "
                            "cada dap.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        trees = arcpy.Parameter(
            displayName=u"Layer contendo as árvores",
            name="trees",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")
        trees.filter.list = ["Point"]

        class_range = arcpy.Parameter(
            displayName="Intervalo de classe",
            name="class_range",
            datatype="Double",
            parameterType="Required",
            direction="Input")
        class_range.value = 2

        minimum_dbh = arcpy.Parameter(
            displayName=u"dap mínimo",
            name="minimum_dbh",
            datatype="Double",
            parameterType="Optional",
            direction="Input")
        minimum_dbh.value = 4

        dbh_fields = arcpy.Parameter(
            displayName="dap",
            name="dbh_fields",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        dbh_fields.parameterDependencies = [trees.name]

        params = [trees,
                  class_range,
                  minimum_dbh,
                  dbh_fields]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        trees = parameters[0].valueAsText
        class_range = parameters[1].value
        minimum_dbh = parameters[2].value
        dbh_fields = parameters[3].valueAsText

        dbh_fields_list = dbh_fields.split(";")

        # List containing all smallest values of dbh.
        small_dbh_list = []
        big_dbh_list = []
        for field in dbh_fields_list:
            dbh_list = [k[0] for k in arcpy.da.SearchCursor(trees, [field])
                        if k[0] not in (None, 0)]

            small_dbh_list.append(min(dbh_list))
            big_dbh_list.append(max(dbh_list))

        # Value for minimum_dbh if it is left empty and value for smallest_dbh.
        if minimum_dbh == "#" or not minimum_dbh:
            minimum_dbh = 1
            smallest_dbh = min(small_dbh_list)
        else:
            smallest_dbh = minimum_dbh

        # Creates a list with the upper limits of the diameter classes.
        first_class = int(smallest_dbh + class_range)
        last_class = int(max(big_dbh_list) + class_range + 1)
        c_range = int(class_range)
        range_of_limits = range(first_class, last_class, c_range)
        class_limits = [k for k in range_of_limits]

        # Definition of the function that classifies dbh.
        def dbh_class(dbh):
            if dbh < minimum_dbh:
                return minimum_dbh + class_range / 2.0
            else:
                for k in class_limits:
                    if dbh < k and (dbh >= k - class_range):
                        return k - class_range / 2.0

        for field in dbh_fields_list:
            arcpy.AddMessage(u"Classificando os dap para " + field + "_C")
            arcpy.AddField_management(trees, field + "_C", "FLOAT")
            with arcpy.da.UpdateCursor(trees, [field, field + "_C"]) as cursor:
                for row in cursor:
                    if row[0] in (None, 0):
                        row[1] = row[0]
                        cursor.updateRow(row)
                    else:
                        row[1] = dbh_class(row[0])
                        cursor.updateRow(row)
        return
