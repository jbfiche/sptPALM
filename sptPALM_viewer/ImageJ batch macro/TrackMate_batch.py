import sys
import os
from ij import IJ
from ij import WindowManager
from glob import glob
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO, DisplaySettings
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.visualization.table import TrackTableView
from fiji.plugin.trackmate.visualization.table import AllSpotsTableView
from fiji.plugin.trackmate.action import ExportTracksToXML
from java.io import File

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')

# ----------------------------------------------
# Indicate the parameters for the spot detection
# ----------------------------------------------
data_folder = "C:\Users\sCMOS-1\Desktop\Test_data"
object_diameter = 5 # RADIUS of the searched objects in pixels
quality_threshold = 2.5
use_median_filetering = True
sub_pixel_localization = True
perform_filtering = False
filtering_threshold = None

# ----------------------------------------------------
# Indicate the parameters for the track reconstruction
# ----------------------------------------------------
linking_max_distance = 2 # Distance in pixels
gap_closing_distance = 5 # Distance in pixels
max_frame_gap = 3

# --------------------------------------------------------
# Look for all the tif files within the selected directory
# --------------------------------------------------------
data_files = glob(os.path.join(data_folder,"**", "*.tif"))

for data_file in data_files:

	imp = IJ.openImage(data_file)
	imp.show()
	print('Analyzing : ' + data_file)
	
	# Define the image properties
	# ---------------------------
	NFrames = imp.getNFrames()
	NSlices = imp.getNSlices()	
	print('The image is composed of ' + str(NFrames) + ' frames and '+ str(NSlices) + ' slices')

	if NFrames == 1 and NSlices == 1:
		print('ERROR : the data is composed of a single image - analysis is aborted')
		continue
	elif NFrames == 1 and NSlices > 1:
		updateProperties = "".join(["channels=1 slices=1 frames=",str(NSlices)," unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1"])
		IJ.run("Properties...", updateProperties);
	elif NFrames > 1 and NSlices == 1:
		updateProperties = "".join(["channels=1 slices=1 frames=",str(NFrames)," unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1"])
		IJ.run("Properties...", updateProperties);
	else:
		continue
	
	#----------------------------
	# Create the model object now
	#----------------------------
	
	# Some of the parameters we configure below need to have
	# a reference to the model at creation. So we create an
	# empty model now.
	model = Model()
	
	# Send all messages to ImageJ log window.
	model.setLogger(Logger.IJ_LOGGER)
	
	#------------------------
	# Prepare settings object
	#------------------------
	
	settings = Settings(imp)
	
	# Configure detector - We use the Strings for the keys
	settings.detectorFactory = LogDetectorFactory()
	settings.detectorSettings = {
	    'DO_SUBPIXEL_LOCALIZATION' : sub_pixel_localization,
	    'RADIUS' : float(object_diameter),
	    'TARGET_CHANNEL' : 1,
	    'THRESHOLD' : float(quality_threshold),
	    'DO_MEDIAN_FILTERING' : use_median_filetering,
	}  
	
	# Configure spot filters - Classical filter on quality
	if perform_filtering and (filtering_threshold is not None):
		filter1 = FeatureFilter('QUALITY', float(filtering_threshold), True)
		settings.addSpotFilter(filter1)
	
	# Configure tracker
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
	settings.trackerFactory.getDefaultSettings() # to use as a template
	settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
	settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
	settings.trackerSettings['LINKING_MAX_DISTANCE'] = float(linking_max_distance)
	settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = float(gap_closing_distance)
	settings.trackerSettings['MAX_FRAME_GAP'] = int(max_frame_gap)
	
	# Add ALL the feature analyzers known to TrackMate. They will 
	# yield numerical features for the results, such as speed, mean intensity etc.
	settings.addAllAnalyzers()
	
	#-------------------
	# Instantiate plugin
	#-------------------

	trackmate = TrackMate(model, settings)
	
	#--------
	# Process
	#--------
	ok = trackmate.checkInput()
	if not ok:
	    sys.exit(str(trackmate.getErrorMessage()))
	
	ok = trackmate.process()
	if not ok:
	    sys.exit(str(trackmate.getErrorMessage()))
	
	#------------------------------
	# Save the results as csv or xml
	#-------------------------------
	filename = os.path.splitext(os.path.basename(data_file))[0]
	filepath = os.path.dirname(data_file)
	output_filename_track = File(filepath, filename + '_track.csv')
	output_filename_track_xml = File(filepath, filename + '_track.xml')
	output_filename_spots = File(filepath, filename + '_all_spots.csv')
	
	sm = SelectionModel(trackmate.getModel())
	ds = DisplaySettings()
	trackTableView = TrackTableView(trackmate.getModel(), sm, ds)
	trackTableView.getSpotTable().exportToCsv(output_filename_track)
	ExportTracksToXML.export(model, settings, output_filename_track_xml)
	
	spotsTableView = AllSpotsTableView(trackmate.getModel(), sm, ds)
	spotsTableView.exportToCsv(output_filename_spots.getAbsolutePath())
	
	# ---------------
	# Close the image
	# ---------------
	imp.close()

# Echo results with the logger we set at start:
model.getLogger().log( str( model ) )
