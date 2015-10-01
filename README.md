ML_AudMotion
============
Paige Lab Movement Analysis programs
Designed to open analyse and export auditory motion experiment files
One program for free-field stimuli a second for headphones

---
*1. Program name is 'SLab_Movement'*

*2. Load Data file *

2.1 File dialog box prompts for you to identify the file for analysis. Default directory is GSM on share drive.

2.2 Auto-detect data file type and Cycles or ramp (default)

2.3 EOG and Elmar (automatic by datastructure)

2.4 Parse XML to Matlab Structure file - this uses function of Falkena, Wanner, Smirnov & Mo


*3. Display Data - Four or three plots*

3.1 Raw Position- Arm(Magenta), Eye(Blue), Gaze(Green) Head Yaw (Yellow)

3.2 Velocity
* Time = 500ms is the start of the arm motion.
* Arm motion start time is determined from velocity representation, first time at which the abs value of arm velocity is 20% of the max abs arm velocity.
* This plot is not shown for Cycles data

3.3 Candidate saccades - Plot shows in white the fast component of velocity and there is a magenta horizontal line representing the threshold. In red are saccaded above the threshold

3.4 Processed - desaccaded etc
* very good correspondance between start of movement from velocity and the by-eye start of arm movement.
* Yellow dashed line is the start of head motion, also as 10% of max head motion.
* blue dahed is start of gaze motion, set at gaze velocity great than 20 deg/sec

3.5 Note that the eye signal is taken from the 'good' eye. For EOG there is only one eye. For Elmar, the Good Eye is the eye signal from the hemifield that the speaker is in. This has the most impact on calibration, when the bad eye can saturate and throw off the scaling

3.6 In cycles mode the third plot is ommitted and the fourth plot is twice as tall as for ramps

*4. Desaccading*

4.1 Skirt - use 20ms 'skirt' around and saccade points 

4.2 Threshold is dynamic - difference between fast changes in velocity and gradual.
Operationally, the velocity is calculated twice, once with a window of 20 points, then again with 500 points. The absolute difference plot shows well-defined peaks when saccades occur. The default threshold is set at 30% of the maximum sacade velocity, and this is adjustable on a per trial basis. This is shown graphically in the  

*5. Export*

5.1 Peak pursuit, fastest segment velocity

5.2 Time to peak is the middle of the fastest segment

5.3 New metric: Peak and time to peak in 50, look at the middle 50% of the movement and check out the fastest segment there

5.5 Export for cycles saves the datastructure that has columns of Speaker, eyes, desaccaded eyes with NaN for saccades, head movement
