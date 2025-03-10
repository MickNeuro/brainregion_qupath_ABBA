/**
 * Script: Performs cell detection on the root node annotation for both AF647 and DAPI channels,
 * calculates cell counts, fractions, area of AF647 cells, total region area, and mean intensity
 * separately for left and right brain regions, and writes results to a single CSV file.
 * Original Authors: Mick de Koning (cell detection), Modified and Combined: March 10, 2025
 * 
 * This version works with ABBA-aligned brain regions and splits measurements
 * within each region into left and right hemispheres based on the midline position.
 */

// Required imports
import qupath.lib.objects.PathObjectTools
import qupath.lib.regions.RegionRequest
import qupath.imagej.tools.IJTools
import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.roi.ROIs
import qupath.lib.roi.GeometryTools
import qupath.lib.objects.PathObjects
import qupath.lib.roi.interfaces.ROI
import qupath.lib.geom.Point2

// Ensure the image type is set to fluorescence 
setImageType('FLUORESCENCE')

// Clear previous detections
clearDetections()

// Retrieve the root annotation
def rootAnnotationName = 'root'  // Annotation name is lowercase 'root'
def root = getAnnotationObjects().find { it.getName().equalsIgnoreCase(rootAnnotationName) }

// Check if root annotation exists
if (root == null) {
    print("Root annotation named '${rootAnnotationName}' not found. Please check the available annotations.")
    return
}

// Get current image data and server
def imageData = getCurrentImageData()
def server = imageData.getServer()

// Find the bounding box of all annotations to determine the aligned brain region
def allAnnotations = getAnnotationObjects().findAll { it != root }
double minX = Double.MAX_VALUE
double maxX = Double.MIN_VALUE

allAnnotations.each { annotation ->
    def roi = annotation.getROI()
    def bounds = roi.getBoundsX()
    if (bounds < minX) minX = bounds
    double rightEdge = bounds + roi.getBoundsWidth()
    if (rightEdge > maxX) maxX = rightEdge
}

// Calculate the midline of the aligned brain region
double brainMidlineX = (minX + maxX) / 2
print("Calculated brain midline at x-coordinate: " + brainMidlineX)

// Function to determine if a cell is in the left or right hemisphere
String determineHemisphere(double x, double midlineX) {
    return x < midlineX ? "Left" : "Right"
}

// Helper function to measure mean intensity for a given ROI on a specified channel
def measureIntensityForROI = { roi, channelName, downsample = 1.0 ->
    try {
        def channels = server.getMetadata().getChannels()
        def channelNames = channels.collect { it.getName() }
        def channelIndex = channelNames.indexOf(channelName) + 1
        
        if (channelIndex == 0) {
            print "Channel '${channelName}' not found."
            return 0
        }
        
        def request = RegionRequest.createInstance(server.getPath(), downsample, roi)
        def pathImage = IJTools.convertToImagePlus(server, request)
        def imp = pathImage.getImage()
        
        // Set the desired channel (ImageJ channels are 1-based)
        imp.setC(channelIndex)
        def stats = imp.getProcessor().getStatistics()
        return stats.mean
    } catch (Exception e) {
        print("Error measuring intensity: " + e.getMessage())
        return 0
    }
}

// Duplicate root annotation for detection
def regionMeasured = PathObjectTools.transformObject(root, null, true)
regionMeasured.setName("Measured Region")
addObject(regionMeasured)
setSelectedObject(regionMeasured)

// Define cell detection parameters for AF647 and DAPI channels
def detectionParamsAF647 = '{"detectionImage": "AF647", "requestedPixelSizeMicrons": 0.5, "backgroundRadiusMicrons": 8.0, "medianRadiusMicrons": 0.0, "sigmaMicrons": 1.5, "minAreaMicrons": 10.0, "maxAreaMicrons": 400.0, "threshold": 150.0, "watershedPostProcess": true, "cellExpansionMicrons": 0.0, "includeNuclei": true, "smoothBoundaries": true, "makeMeasurements": true}'
def detectionParamsDAPI = '{"detectionImage": "DAPI", "requestedPixelSizeMicrons": 2.0, "backgroundRadiusMicrons": 8.0, "medianRadiusMicrons": 0.0, "sigmaMicrons": 1.5, "minAreaMicrons": 10.0, "maxAreaMicrons": 400.0, "threshold": 200.0, "watershedPostProcess": true, "cellExpansionMicrons": 0.0, "includeNuclei": true, "smoothBoundaries": true, "makeMeasurements": true}'

// Run cell detection for AF647 channel
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', detectionParamsAF647)
def af647Cells = getDetectionObjects()
clearDetections() // Clear detections before running the next channel

// Run cell detection for DAPI channel
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', detectionParamsDAPI)
def dapiCells = getDetectionObjects()

// Remove the temporary annotation used for detection but keep the detected cells
removeObject(regionMeasured, true)

// Add detected cells to the original root annotation
root.addChildObjects(af647Cells + dapiCells) // Using addChildObjects instead of deprecated addPathObjects

// Update the view
fireHierarchyUpdate()

// Data structures for storing results
def regionData = [:]

// Process each annotation region
def annotations = getAnnotationObjects()

// Create a vertical line at the midline for ROI splitting
def imageHeight = server.getHeight()
def midlinePoints = [new Point2(brainMidlineX, 0), new Point2(brainMidlineX, imageHeight)]

annotations.each { region ->
    if (region != root) {
        def regionROI = region.getROI()
        def regionName = region.getName()
        
        // Initialize data structure for this region if not exists
        if (!regionData.containsKey(regionName)) {
            regionData[regionName] = [
                'Left': [
                    'countAF647': 0,
                    'countDAPI': 0,
                    'areaAF647': 0.0,
                    'areaTotal': 0.0,
                    'intensityAF647': 0.0
                ],
                'Right': [
                    'countAF647': 0,
                    'countDAPI': 0, 
                    'areaAF647': 0.0,
                    'areaTotal': 0.0,
                    'intensityAF647': 0.0
                ]
            ]
        }
        
        // Split the region ROI into left and right parts using the midline
        try {
            // Get region geometry
            def geometry = regionROI.getGeometry()
            
            // Create bounding boxes for left and right sides of the brain
            def imageBounds = [0, 0, server.getWidth(), server.getHeight()]
            def leftBox = GeometryTools.createRectangle(0, 0, brainMidlineX, imageHeight)
            def rightBox = GeometryTools.createRectangle(brainMidlineX, 0, server.getWidth() - brainMidlineX, imageHeight)
            
            // Intersect region with left and right boxes
            def leftGeom = geometry.intersection(leftBox)
            def rightGeom = geometry.intersection(rightBox)
            
            // Create ROIs for the left and right parts
            def leftROI = null
            def rightROI = null
            
            if (!leftGeom.isEmpty()) {
                leftROI = GeometryTools.geometryToROI(leftGeom, regionROI.getImagePlane())
                regionData[regionName]['Left']['areaTotal'] = leftGeom.getArea()
                
                // Measure intensity for left part
                if (leftROI != null) {
                    regionData[regionName]['Left']['intensityAF647'] = measureIntensityForROI(leftROI, "AF647", 1.0)
                }
            }
            
            if (!rightGeom.isEmpty()) {
                rightROI = GeometryTools.geometryToROI(rightGeom, regionROI.getImagePlane())
                regionData[regionName]['Right']['areaTotal'] = rightGeom.getArea()
                
                // Measure intensity for right part
                if (rightROI != null) {
                    regionData[regionName]['Right']['intensityAF647'] = measureIntensityForROI(rightROI, "AF647", 1.0)
                }
            }
            
            // Process AF647 cells
            af647Cells.each { cell ->
                def cellROI = cell.getROI()
                double cellX = cellROI.getCentroidX()
                double cellY = cellROI.getCentroidY()
                
                if (regionROI.contains(cellX, cellY)) {
                    String hemisphere = determineHemisphere(cellX, brainMidlineX)
                    regionData[regionName][hemisphere]['countAF647'] += 1
                    regionData[regionName][hemisphere]['areaAF647'] += cellROI.getArea()
                }
            }
            
            // Process DAPI cells
            dapiCells.each { cell ->
                def cellROI = cell.getROI()
                double cellX = cellROI.getCentroidX() 
                double cellY = cellROI.getCentroidY()
                
                if (regionROI.contains(cellX, cellY)) {
                    String hemisphere = determineHemisphere(cellX, brainMidlineX)
                    regionData[regionName][hemisphere]['countDAPI'] += 1
                }
            }
            
        } catch (Exception e) {
            print("Error processing region ${regionName}: " + e.getMessage())
        }
    }
}

// Retrieve the image name
def imageName = server.getMetadata().getName().replaceAll('\\.\\w+$', '')

// Determine the desktop path dynamically
def userHome = System.getProperty("user.home")
def desktopPath = new File(userHome, "Desktop")
def outputFile = new File(desktopPath, imageName + '_cellcount_intesity.csv')

// Helper function to properly format CSV fields (escaping commas in region names)
def formatCSVField = { field ->
    // Check if the field contains commas, quotes, or newlines
    if (field ==~ /.*[,"\n].*/) {
        // Escape any quotes by doubling them and enclose in quotes
        return '"' + field.replaceAll('"', '""') + '"'
    }
    return field
}

// Write all measurements to a single CSV file
outputFile.withWriter { writer ->
    writer.println('Region,Side,Cell Count AF647,Cell Count DAPI,Fraction AF647,Fraction DAPI,Area AF647,Total Region Area,Mean Intensity AF647')
    
    regionData.each { regionName, hemisphereData ->
        hemisphereData.each { hemisphere, data ->
            def totalCount = data['countAF647'] + data['countDAPI']
            def fractionAF647 = totalCount > 0 ? data['countAF647'] / totalCount : 0
            def fractionDAPI = totalCount > 0 ? data['countDAPI'] / totalCount : 0
            
            // Format the region name to handle commas properly
            def escapedRegionName = formatCSVField(regionName)
            
            writer.println("${escapedRegionName},${hemisphere},${data['countAF647']},${data['countDAPI']},${fractionAF647},${fractionDAPI},${data['areaAF647']},${data['areaTotal']},${data['intensityAF647']}")
        }
    }
}

print('Cell counts, fractions, areas, and intensity measurements have been written to ' + imageName + '_combined_measurements.csv on the Desktop.')