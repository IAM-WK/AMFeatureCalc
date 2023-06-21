#@ OpService ops
from ij import IJ, WindowManager
from ij.plugin import Duplicator
from ij.measure import ResultsTable
from net.imagej.ops.image.cooccurrenceMatrix import MatrixOrientation2D
import csv, copy


# Parameters for evaluation ------------------------
# The values of the parameters are set here. Only the first entry in the list is used for calculation.
# The other entries are intended for a parameter variation, so that different settings can be investigated.
ParameterScript = {"ThresholdRadius": [50],  # which pixels are included for the local threshold (given in pixels)
				   "HaralickDistance": [1, 3, 6, 9, 12],  # pixels with the distance X are set in a neighborhood relationship (given in pixels)
				   "HaralickDirection": [MatrixOrientation2D.HORIZONTAL, MatrixOrientation2D.VERTICAL, MatrixOrientation2D.DIAGONAL, MatrixOrientation2D.ANTIDIAGONAL],
				   "MorphoLibPoreShape": ["all", "round", "deformed"],  # should the round pores, the deformed pores or all pores be kept?
				   "MorphoLibJMinPoreSize": [5],  # Pores smaller than this value are filtered because they are at the resolution limit (given in pixels)
				   "MorphoLibJDilateErode": [True, False],  # should erode and dilate be used to determine the pore labels?
				   "SlidingImageHeight": [12, 3, 21],  # this is the height of the image section - only relevant for quasi-continuous curve (given in layers)
				   "SlidingDistance": [1]  # with this distance value, the image section is rasterized over the entire image - only relevant for quasi-continuous curve (given in pixels)
				   }
Thresholdmethods = ["Mean"]  # other methods are also available and can be added to the list. See: https://imagej.net/plugins/auto-local-threshold
MorphoLibJMethodAverage = True  # it is possible to choose between "median" or "average". If false, then the median is used.
ImageResolution = 130  # resolution of the image (given in pixels per mm (pixel/mm))
LayerHeight = 0.26  # only relevant for quasi-continuous curve (given in mm)
# Here you can define which parameter is to be varied. For this parameter, iteration is performed over the list 'ParameterScript'
# If no parameter variation is to be performed, then set ChangeParameter to None --> ChangeParameter = None
ChangeParameter = "SlidingImageHeight"
# ChangeParameter = None
# If only one value for a feature is to be calculated for the entire image, then set EntireImage to True --> EntireImage = True
# Otherwise a quasi continuous curve is calculated for the Haralick Entropy and porosity
EntireImage = False
# Path and names -----------------------------------
FilePath = "D:\Users\XXX"  # path to the image to be evaluated
FileNames = ["imgXXX"]  # list of file names for the images to be evaluated
FileFormat = ".tif"  # image file extension
SavePath = "D:\Users\XXX"  # path to the location where the evaluation is to be saved
SaveNameCSV = "Porosity_HaralickEntropy_MorphoLibJ"  # filename for the evaluation. Furthermore, the image name and maybe the parameter name to be varied are appended to the filename
# --------------------------------------------------

# Image algorithms methods -------------------------
def CropImage(Image, PosHeight, Height):
	width = Image.width
	Image.setRoi(0, PosHeight, width, Height)
	return Image.crop()

def Binarise(Image, method, radius):
	TempImage = Duplicator().run(Image)
	IJ.run(TempImage, "Auto Local Threshold", "method={} radius={} parameter_1=0 parameter_2=0 white".format(method, radius))
	return TempImage

def BinariseInvert(Image, method, radius):
	TempImage = Duplicator().run(Image)
	IJ.run(TempImage, "Auto Local Threshold", "method={} radius={} parameter_1=0 parameter_2=0".format(method, radius))
	return TempImage

def PoresList(Image, DilateErode):
	global ImageResolution
	if DilateErode:
		IJ.run("Options...", "iterations=1 count=1")  # Set default background for the Dilate operation
		IJ.run(Image, "Dilate", "")  # Reduce pore area to interrupt connections
	IJ.run(Image, "Set Scale...", "distance={} known=1 unit=mm".format(ImageResolution))
	IJ.run(Image, "Connected Components Labeling", "connectivity=8 type=[16 bits]")
	if DilateErode:
		IJ.run(IJ.getImage(), "Dilate Labels", "radius=1")
	IJ.run(IJ.getImage(), "Analyze Regions", "pixel_count area perimeter circularity centroid equivalent_ellipse ellipse_elong.")
	MorphoTabResult = ResultsTable.getResultsTable(IJ.getImage().getTitle() + "-Morphometry")
	# Save data into vector
	Data = {}
	for column in range(MorphoTabResult.getLastColumn() + 1):
		Data[MorphoTabResult.getColumnHeading(column)] = MorphoTabResult.getColumn(column)
	# Close Windows
	WindowManager.getWindow(IJ.getImage().getTitle() + "-Morphometry").close()
	IJ.getImage().close()
	if DilateErode:
		IJ.getImage().close()
	return Data

# Data processing methods -------------------------
def FilterPoresList(Data, Filterfeature, FilterValueMin, FilterValueMax):
	pointer = 0
	while True:
		try:
			if Data[Filterfeature][pointer] < float(FilterValueMin):
				del Data[Filterfeature][pointer]
				for keyname in Data:
					if keyname != Filterfeature:
						del Data[keyname][pointer]
			elif Data[Filterfeature][pointer] > float(FilterValueMax):
				del Data[Filterfeature][pointer]
				for keyname in Data:
					if keyname != Filterfeature:
						del Data[keyname][pointer]
			else:
				pointer += 1
		except IndexError:
			break
	return Data

def ChangePoresValues(Data, Changefeature, FilterValueMax):
	pointer = 0
	while True:
		try:
			if Data[Changefeature][pointer] > FilterValueMax:
				Data[Changefeature][pointer] = FilterValueMax
			pointer += 1
		except IndexError:
			break
	return Data

def DetermineFeatureAverage(Data, ROIMin, ROIMax, ROIArea):
	ListFeatureAverage = {}
	ExcludeFeature = ["PixelCount", "Centroid.X", "Centroid.Y", "Ellipse.Center.X", "Ellipse.Center.Y"]
	for Feature in Data:
		if Feature not in ExcludeFeature:
			ListFeatureAverage[Feature] = 0.0
	pointer = 0
	PoreAmount = 0.0
	for value in Data["Centroid.Y"]:
		if value <= ROIMax and value >= ROIMin:
			for Feature in ListFeatureAverage:
				ListFeatureAverage[Feature] += Data[Feature][pointer]
			PoreAmount += 1
		pointer += 1
	if PoreAmount != 0:
		for Feature in ListFeatureAverage:
			ListFeatureAverage[Feature] = ListFeatureAverage[Feature] / PoreAmount  # mean value
	ListFeatureAverage["PoreAmount"] = PoreAmount / ROIArea
	return ListFeatureAverage

def DetermineFeatureMedian(Data, ROIMin, ROIMax, ROIArea):
	ListFeatureMedian = {}
	ExcludeFeature = ["PixelCount", "Centroid.X", "Centroid.Y", "Ellipse.Center.X", "Ellipse.Center.Y"]
	for Feature in Data:
		if Feature not in ExcludeFeature:
			ListFeatureMedian[Feature] = []
	pointer = 0
	PoreAmount = 0.0
	for value in Data["Centroid.Y"]:
		if value <= ROIMax and value >= ROIMin:
			for Feature in ListFeatureMedian:
				ListFeatureMedian[Feature].append(Data[Feature][pointer])
			PoreAmount += 1
		pointer += 1
	for Feature in ListFeatureMedian:
		ListFeatureMedian[Feature].sort()
		if len(ListFeatureMedian[Feature]) % 2 != 0:
			m = int((len(ListFeatureMedian[Feature]) + 1) / 2 - 1)
			ListFeatureMedian[Feature] = ListFeatureMedian[Feature][m]
		else:
			m1 = int(len(ListFeatureMedian[Feature]) / 2 - 1)
			m2 = int(len(ListFeatureMedian[Feature]) / 2)
			ListFeatureMedian[Feature] = (ListFeatureMedian[Feature][m1] + ListFeatureMedian[Feature][m2]) / 2
	ListFeatureMedian["PoreAmount"] = PoreAmount / ROIArea
	return ListFeatureMedian

def DeterminePorosity(Image):
	ListBin = Image.getProcessor().getHistogram()
	white = ListBin[255]
	black = ListBin[0]
	porosity = round(black/(float(white + black)), 3)
	return porosity

def HaralickFeatureASM(Image, GreyValues, Distance, Direction):
	ASM = ops.run("haralick.asm", Image, GreyValues, Distance, Direction).get()
	return ASM

def HaralickFeatureDifferenceEntropy(Image, GreyValues, Distance, Direction):
	DifferenceEntropy = ops.run("haralick.differenceEntropy", Image, GreyValues, Distance, Direction).get()
	return DifferenceEntropy

def HaralickFeatureEntropy(Image, GreyValues, Distance, Direction):
	Entropy = ops.run("haralick.entropy", Image, GreyValues, Distance, Direction).get()
	return Entropy

def HaralickFeatureTextureHomogeneity(Image, GreyValues, Distance, Direction):
	TextureHomogeneity = ops.run("haralick.textureHomogeneity", Image, GreyValues, Distance, Direction).get()
	return TextureHomogeneity

# Write data in a readable format ---------------------
def WriteCSV(Tab, FileName, CSVNameTag):
	# Write values in CSV-file
	global SavePath, SaveNameCSV
	PathAndFilename = SavePath + "\\" + SaveNameCSV + "_" + FileName + "_" + CSVNameTag
	with open(PathAndFilename + ".csv", mode='w') as csv_file:
		csv_writer = csv.writer(csv_file, delimiter=';', lineterminator='\n')
		for row in Tab:
			csv_writer.writerow(row)

def WriteCSVTab(DataVector):
	Tab = [[]]
	CreateRow = True
	for keyname in DataVector["Column_number"]:
		# Header
		Tab[0].append(keyname)
		if CreateRow:  # Create row if not present
			for DataPoint in DataVector[keyname]:
				Tab.append([DataPoint])
		else:
			Pos = 1
			for DataPoint in DataVector[keyname]:
				Tab[Pos].append(DataPoint)
				Pos += 1
		CreateRow = False
	return Tab

def WriteCSVTabPara(Tab):
	global ParameterScript
	Tab[0].append("")
	Tab[1].append("")
	for ParameterName in ParameterScript:
		Tab[0].append(ParameterName)
		Tab[1].append(ParameterScript[ParameterName][0])
	return Tab

# Image sliding and determination of all features ------------
def DetermineCropPosition(HaralickData, StartPos, Height):
	global ParameterScript
	Pos = StartPos
	while Pos > 0:
		HaralickData["ROIPosition"].append(Pos)
		HaralickData["Position"].append(int(Pos + Height / 2))  # This is in the centre of the analysed image section
		Pos = Pos - ParameterScript["SlidingDistance"][0]
	return HaralickData

def ImageSliding(Img):
	# Slide over the entire image and generate data
	global ParameterScript, Thresholdmethods, MorphoLibJMethodAverage, ImageResolution, LayerHeight
	SlidingImageHeightPixel = int(LayerHeight * ImageResolution * ParameterScript["SlidingImageHeight"][0])
	SlidingImageArea = (Img.width * SlidingImageHeightPixel)/(ImageResolution * ImageResolution)
	ROIPos = Img.height - SlidingImageHeightPixel
	DataVector = DetermineCropPosition({"Column_number": ["ROIPosition", "Position"], "ROIPosition": [], "Position": []}, ROIPos, SlidingImageHeightPixel)

	# Convert image
	IJ.run(Img, "8-bit", "")

	# For different thresholdmethods
	for method in Thresholdmethods:
		print("Start with method {} ...".format(method))

		# Binarise image
		TempImg = Binarise(Img, method, ParameterScript["ThresholdRadius"][0])

		# MorphoLibJ - determine pore labels and filter pores
		TempImgInvert = BinariseInvert(Img, method, ParameterScript["ThresholdRadius"][0])
		MorphoLibJData = PoresList(TempImgInvert, ParameterScript["MorphoLibJDilateErode"][0])
		MorphoLibJData = FilterPoresList(MorphoLibJData, "PixelCount", ParameterScript["MorphoLibJMinPoreSize"][0], (Img.width / 4) * (Img.height / 4))  # filter small pores
		MorphoLibJData = ChangePoresValues(MorphoLibJData, "Circularity", 1.0)  # Change value of circularity
		if ParameterScript["MorphoLibPoreShape"][0] == "all":
			pass
		elif ParameterScript["MorphoLibPoreShape"][0] == "round":
			MorphoLibJData = FilterPoresList(MorphoLibJData, "Circularity", 0.5, 1.0)  # filter deformed pores
		elif ParameterScript["MorphoLibPoreShape"][0] == "deformed":
			MorphoLibJData = FilterPoresList(MorphoLibJData, "Circularity", 0.0, 0.5)  # filter round pores

		# Create new data columns - name of the data columns
		DataVector["Column_number"].append("Porosity(-) " + method)
		DataVector["Column_number"].append("Entropy(-) " + method)
		DataVector["Column_number"].append("Area(mm^2) " + method)
		DataVector["Column_number"].append("Perimeter " + method)
		DataVector["Column_number"].append("Circularity(-) " + method)
		DataVector["Column_number"].append("Ellipse.Radius1 " + method)
		DataVector["Column_number"].append("Ellipse.Radius2 " + method)
		DataVector["Column_number"].append("Ellipse.Orientation " + method)
		DataVector["Column_number"].append("Ellipse.Elong(-) " + method)
		DataVector["Column_number"].append("PoreAmount(1/mm^2) " + method)

		# Create new data columns in DataVector
		for column in DataVector["Column_number"]:
			if column not in DataVector:
				DataVector[column] = []

		# Determine porosity, entropy and features from MorphoLibJData for actual "ROI-position"
		for ROIPos in DataVector["ROIPosition"]:
			print(ROIPos)
			CropImg = CropImage(TempImg, ROIPos, SlidingImageHeightPixel)  # Crop image

			# Show image with ROI and show position of the cropped image
			# TempImg.show()

			# Determine porosity
			DataVector["Porosity(-) " + method].append(str(DeterminePorosity(CropImg)).replace(".", ","))

			# Determine entropy
			DataVector["Entropy(-) " + method].append(str(HaralickFeatureEntropy(CropImg, 2, ParameterScript["HaralickDistance"][0], ParameterScript["HaralickDirection"][0])).replace(".", ","))

			# Determine features from MorphoLibJData
			SlidingImageMin = float(ROIPos) / ImageResolution
			SlidingImageMax = float(ROIPos + SlidingImageHeightPixel) / ImageResolution
			if MorphoLibJMethodAverage:
				AllFeature = DetermineFeatureAverage(MorphoLibJData, SlidingImageMin, SlidingImageMax, SlidingImageArea)
			else:
				AllFeature = DetermineFeatureMedian(MorphoLibJData, SlidingImageMin, SlidingImageMax, SlidingImageArea)
			DataVector["Area(mm^2) " + method].append(str(AllFeature["Area"]).replace(".", ","))
			DataVector["Perimeter " + method].append(str(AllFeature["Perimeter"]).replace(".", ","))
			DataVector["Circularity(-) " + method].append(str(AllFeature["Circularity"]).replace(".", ","))
			DataVector["Ellipse.Radius1 " + method].append(str(AllFeature["Ellipse.Radius1"]).replace(".", ","))
			DataVector["Ellipse.Radius2 " + method].append(str(AllFeature["Ellipse.Radius2"]).replace(".", ","))
			DataVector["Ellipse.Orientation " + method].append(str(AllFeature["Ellipse.Orientation"]).replace(".", ","))
			DataVector["Ellipse.Elong(-) " + method].append(str(AllFeature["Ellipse.Elong"]).replace(".", ","))
			DataVector["PoreAmount(1/mm^2) " + method].append(str(AllFeature["PoreAmount"]).replace(".", ","))

		# Close image
		TempImg.close()
	# Output of the data
	return DataVector

# Only entire image and determination of all features --------
def ImageEntire(Img):
	# Use the entire Image and generate data
	global ParameterScript, Thresholdmethods, MorphoLibJMethodAverage, ImageResolution
	ImageArea = (Img.width * Img.height)/(ImageResolution * ImageResolution)
	DataVector = {"Column_number": ["ROIPosition", "Position"], "ROIPosition": [0], "Position": [int(Img.height / 2)]}

	# Convert image
	IJ.run(Img, "8-bit", "")

	# For different thresholdmethods
	for method in Thresholdmethods:
		print("Start with method {} ...".format(method))

		# Binarise image
		TempImg = Binarise(Img, method, ParameterScript["ThresholdRadius"][0])

		# MorphoLibJ - determine pore labels and filter pores
		TempImgInvert = BinariseInvert(Img, method, ParameterScript["ThresholdRadius"][0])
		MorphoLibJData = PoresList(TempImgInvert, ParameterScript["MorphoLibJDilateErode"][0])
		MorphoLibJData = FilterPoresList(MorphoLibJData, "PixelCount", ParameterScript["MorphoLibJMinPoreSize"][0], ((Img.width * Img.height)/1))  # filter small pores
		MorphoLibJData = ChangePoresValues(MorphoLibJData, "Circularity", 1.0)  # Change value of circularity
		if ParameterScript["MorphoLibPoreShape"][0] == "all":
			pass
		elif ParameterScript["MorphoLibPoreShape"][0] == "round":
			MorphoLibJData = FilterPoresList(MorphoLibJData, "Circularity", 0.5, 1.0)  # filter deformed pores
		elif ParameterScript["MorphoLibPoreShape"][0] == "deformed":
			MorphoLibJData = FilterPoresList(MorphoLibJData, "Circularity", 0.0, 0.5)  # filter round pores

		# Create new data columns - name of the data columns
		DataVector["Column_number"].append("Porosity(-) " + method)
		DataVector["Column_number"].append("Entropy(-) " + method)
		DataVector["Column_number"].append("Area(mm^2) " + method)
		DataVector["Column_number"].append("Perimeter " + method)
		DataVector["Column_number"].append("Circularity(-) " + method)
		DataVector["Column_number"].append("Ellipse.Radius1 " + method)
		DataVector["Column_number"].append("Ellipse.Radius2 " + method)
		DataVector["Column_number"].append("Ellipse.Orientation " + method)
		DataVector["Column_number"].append("Ellipse.Elong(-) " + method)
		DataVector["Column_number"].append("PoreAmount(1/mm^2) " + method)

		# Create new data columns in DataVector
		for column in DataVector["Column_number"]:
			if column not in DataVector:
				DataVector[column] = []

		# Determine porosity
		DataVector["Porosity(-) " + method].append(str(DeterminePorosity(TempImg)).replace(".", ","))

		# Determine entropy
		DataVector["Entropy(-) " + method].append(str(HaralickFeatureEntropy(TempImg, 2, ParameterScript["HaralickDistance"][0], ParameterScript["HaralickDirection"][0])).replace(".", ","))

		# Determine features from MorphoLibJData
		if MorphoLibJMethodAverage:
			AllFeature = DetermineFeatureAverage(MorphoLibJData, 0.0, float(Img.height), ImageArea)
		else:
			AllFeature = DetermineFeatureMedian(MorphoLibJData, 0.0, float(Img.height), ImageArea)
		DataVector["Area(mm^2) " + method].append(str(AllFeature["Area"]).replace(".", ","))
		DataVector["Perimeter " + method].append(str(AllFeature["Perimeter"]).replace(".", ","))
		DataVector["Circularity(-) " + method].append(str(AllFeature["Circularity"]).replace(".", ","))
		DataVector["Ellipse.Radius1 " + method].append(str(AllFeature["Ellipse.Radius1"]).replace(".", ","))
		DataVector["Ellipse.Radius2 " + method].append(str(AllFeature["Ellipse.Radius2"]).replace(".", ","))
		DataVector["Ellipse.Orientation " + method].append(str(AllFeature["Ellipse.Orientation"]).replace(".", ","))
		DataVector["Ellipse.Elong(-) " + method].append(str(AllFeature["Ellipse.Elong"]).replace(".", ","))
		DataVector["PoreAmount(1/mm^2) " + method].append(str(AllFeature["PoreAmount"]).replace(".", ","))

		# Close image
		TempImg.close()
	# Output of the data
	return DataVector

# Execute script -----------------------------------
ParameterScriptBackup = copy.deepcopy(ParameterScript)
for FileName in FileNames:
	Img = IJ.openImage(FilePath + "\\" + FileName + FileFormat)
	print(FilePath + "\\" + FileName + FileFormat)
	while True:
		if EntireImage:
			Data = ImageEntire(Img)
		else:
			Data = ImageSliding(Img)
		# Transfer data into table form
		CSVTab = WriteCSVTab(Data)
		# Write the parameters in the table
		CSVTab = WriteCSVTabPara(CSVTab)
		if ChangeParameter is not None:
			CSVNameTag = ChangeParameter + str(ParameterScript[ChangeParameter][0])
			WriteCSV(CSVTab, FileName, CSVNameTag)  # Write file
			print("Calculation with parameter set finished")
			# Change Parameter
			if len(ParameterScript[ChangeParameter]) != 1:
				del ParameterScript[ChangeParameter][0]
			else:
				break
		else:
			CSVNameTag = ""
			WriteCSV(CSVTab, FileName, CSVNameTag)  # Write file
			break
	print("Calculation with image finished")
	ParameterScript = copy.deepcopy(ParameterScriptBackup)
print("Calculation completely finished")
