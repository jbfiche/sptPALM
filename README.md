sptPALM_viewer is a Matlab software used to analyze and visualize single particle tracking data. The software does not work on the raw images but on the tracks coordinates reconstructed by either **MTT** or **TrackMate** software. 

sptPALM_viewer was written in Matlab2019a and was tested on either Windows 10 or Linux (Ubuntu 18.04.3) running computers. Compiled version for windows is available on demand. The compiled version will necessitate to download compiler for the Mathorks website (installation procedure is attached to the compiled software). This step should require <20 minutes to download and install. 

To launch the GUI, either clone the GitHub repository and execute the **sptPALM_viewer.m** script or launch the compiled version. On normal computer, launching the compiled software will take <1minute. 

Two windows should appear :

- sptPALM_Control_Panel - the GUI
- sptPALM_Display_Panel - the window where all the plots/images will be displayed

The control panel is divided into six sections :

1. **Loading data** where the user is indicating the acquisition parameters as well as the type of data to analyze (MTT or TrackMate) 
2. **Analysis parameters** is where the parameters for the calculation of the MSD/Dinst are selected
3. **Diffusion analysis** is where an analysis can be launched and previous results plotted
4. The **Statistics** section is indicating general information regarding the number of files loaded, the number of tracks analyzed etc.
5. The **Visualization tool** can be used to overlay the acquired images with the reconstructed tracks
6. The **Simulation of sptPALM experiments** section can be used to simulate simple Brownian motion experiments with either one or two populations with specific diffusion properties. 

# Launch a first analyzis

1. To start an analysis, indicate first the **acquisition time in ms** and the **pixel size in µm** for the experiment. Indicate the name of the **Results file** where all the results will be saved (the default name is *MTT_sptPALM_analysis.mat*). 

   

2. Select the format of the track files by selecting either MTT or TrackMate. Then click on **Load MTT/TrackMate files** and indicate the folder where the data are saved. Doing so:
   - all the .mat (for MTT) or .xml (for TrackMate) in the indicated folder will be loaded. 

   - The path to the folder will be indicated at the top of the sptPALM_Control_Panel window

   - The number of files and the number of tracks loaded will be indicated in the Statistics section

   - The cumulative distribution of the *tracks length* will be plotted and saved (**Cumulative_Distribution_LengthStep.png**) in the same folder as well as the *distribution of all trajectories* duration (**Trajectories_duration.png**)

     

3. Set the analysis parameters. By default (see **analysis parameters** section) :
   - The *minimum track length is 7 frames*

   - A *maximum of 3 consecutive frames* for the blinking

   - The calculation method for the instantaneous Diffusion is set to *average MSD* by default - this method is faster and allows for a quick survey of the data. However, the *fit method* is more accurate, particularly for the low values of D. For the fit, by default, the *first 4 points of the MSD* are used.

     

4. Check that the parameter **Results file name** is correct. If you want the program to plot the trajectories at the end of the analysis, check the **Plot trajectories** option (by default, this option is unchecked since this step can be time/memory consuming). Then launch the analysis by pushing the **Analyze trajectories** button. The program will ask whether we want to define a ROI. If yes :

   - Select a movie or a single image associated to the data. If a movie is selected, the program will calculate an average image.

   - The average positions of the tracks are overlaid to the selected image. A ROI can be manually drawn by clicking on the image in order to define the vertices of a polygon. To close the polygon, click-left on the last position and the ROI will be defined automatically.

   - When needed, several ROIs can be defined. 

     

   Note that when a ROI is defined, the track density will be computed within the ROI and displayed in the **Statistics** section. 

   

5. The distribution of instantaneous diffusion coefficient is plotted and the program asks whether we want to fit this distribution with a **single Gaussian model** or a **two Gaussian model**. In the latter case, the user will click where the separation between the two populations is expected. 

   

6. At the end of the analysis, the following documents can be found in the folder :

   - A plot of the diffusion coefficient distribution with the Gaussian fit (Diffusion_distribution_XXX_Method.png)
   - A plot of the MSD (MSD_Curves.png)
   - A file containing all the parameters used for the analysis (Parameters_analysis.txt)
   - A file with all the diffusion coefficients (Saved_Diffusion_Coeff.txt)
   - A .mat file with all the results of the analysis that could be reloaded later in the program for further analysis.



NB: for the demo, loading the files should take less than 1 min for MTT files and a few minutes for TrackMate files. For the analysis, the *Average method* should take < 1min wihout defining any ROI and plotting the trajectories. The *fit* method is however slower and should take a few minutes (<4min). 

# Load previous analysis

To load a previous analysis :

1. Press the **Load previous analysis** button and select the file you want to load. If the selected file has the expected format, all the previous data will loaded and displayed. 
2. A previous analysis can also be loaded by pressing the **Load MTT/TrackMate files** button. In that case, if in the selected folder there is a result file with the same name that **Results file name**, the program will ask whether we want to load the previous analysis or simply start again.
3. The results of the previous analysis can be selected and plotted using the **Plot previous analysis** button. Simply use the menu to select the type of data:
   - **D distribution** plot the distribution of log10(Dinst)
   - **MSD curve** the MSD curve
   - **Trajectories plot (all)** plot the trajectories without specific color coding. If ROIs was defined, only the trajectories within the ROIs are plotted 
   - **Trajectories plot (population)** will plot the trajectories of the two populations obtained fitting the Dinst distribution. The slow population will be plotted in blue, the fast in orange. 

Note that you can launch a new analysis but **ALL** the previous results will be lost if the name of the **Results file name** is not updated. 

# Visualization tool

This tool is working as follows :

1. Select the movie you want to display
2. Select the associated MTT/TrackMate data (should be the same number of frames)
3. Indicate the **Minimum track length** you want to display (1 means that all data will be displayed. 7 will mean that only the tracks that are at least 7 frames long will be displayed)
4. You can adjust the contrast using the **Upper and lower limit** slide-bars
5. You can select the frame by either using the **Select frame** slide bar or directly entering the frame number
6. You can create a rectangular ROI by clicking on the **Create ROI** button.

For each image, the tracks are overlaid with a color code indicating whether the track is starting (blue) or finishing (red). It is also possible to save the data as a .avi movie. In that case, indicates the **first and last** frame numbers as well as the title of the output movie. 

# Simulation of sptPALM experiments

A simulation of a sptPALM experiment can be launched as follows :

1. Indicate the average diffusion coefficient (in µm^2/s)

2. If you are working with a mixture of two populations, indicate the fraction of protein #1 as well (in %)

3. Indicate the number of frames you want to saved (by default 1000 - the program will create several movies, each with a maximum of 1000 frames)

4. Emissions parameters can be modified by pressing the **Change emission parameters** button. Those parameters are used to define the emission properties of the fluorescent proteins as well as the parameters used to analyze the simulated data. 

5. The **acquisition parameters** can also be tuned in order to match the properties of the detector used for the acquisition. The size of the simulated images is also defined there.

6. When **Launch simulation** is pressed, the program is :

   - Calculating all the events and trajectories for each frame using the parameters indicated in **emission parameters**
   - Running MSD and instantaneous diffusion calculation with the parameters indicated in **emission parameters**
   - Finally, it computes all the images using the parameters indicated in **acquisition parameters**.

   

For more information regarding the simulations, check the document sptPALM_simulation.pdf in the folder Doc_simulation.
