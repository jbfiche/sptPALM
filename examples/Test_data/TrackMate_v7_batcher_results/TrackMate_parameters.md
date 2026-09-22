For the batcher, you need to provide an .xml files from a previous TrackMate analysis. The following parameters were used fro TrackMate.

- DoG detector
- Estimated blob diameter = 5 pixels
- Threshold = 3 (but will probably be higher on 'real' data) with the median option
- Simple LAP tracker
- Linking max distance = 2 pixels
- Gap-closing max distance = 5 pixels
- Gap-closing max frame gap = 3

Keep in mind that these parameters were used for stack of images with no-predefined pixel size (in FIJI / Image / Properties, the X/Y dimensions are kept in pixels). If the pixel size is defined (nm, µm), do not forget to adjust the sptPALM_viewer parameters accordingly. 

For this example, each trackMate file corresponds to the analysis of a single movie of 1.000 frames. For sptPALM_viewer, the default parameters were used.