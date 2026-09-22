When loading the movies in FIJI, you can concatenate them using **stacks/tools/concatenate**. After checking the movies have been selected in the right order, **untick the 4D option** and launch the procedure. However, note that for the Vizualization tool, it is recommended that the TrackMate file correspond to the selected movie. Else, the frame number might not match and it will take longer to load the data.

The last step is to define the properties of the movie (either simply keep 1 pixel for X & Y or define the actual pixel size but then you will have to adjust the parameters used for the size of the objects.) 

- DoG detector
- Estimated blob diameter = 5 pixels
- Threshold = 3 (but will probably be higher on 'real' data) with the median option
- Simple LAP tracker
- Linking max distance = 2 pixels
- Gap-closing max distance = 5 pixels
- Gap-closing max frame gap = 3

For this example, the trackMate file corresponds to the analysis of 10 movies of 1.000 frames concatenated into a single stack of images. For sptPALM_viewer, the default parameters were used.