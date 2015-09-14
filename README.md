ML_AudMotion
============
Project plan to integrate Movement Analysis into a few programs

Specs:
- Can read data formatted in 'both' formats
- Initial criteria are to reproduce Basic and Head-free ramp analysis
- Seperate programs for Ramps/Cycles and Headphone Ramps

---
*1. Program name is 'SLab_Movement'*

%2. Auto-detect data file type and Cycles or ramp (default)
%2.1 Head-Free and Head-Fixed ramps and cycles
%2.2 EOG and Elmar (automatic)

3. File dialog box prompts for you to identify the file for analysis
Default directory is GSM on share drive.

*4. Load file*
4.1 Parse XML to Matlab Structure file - this uses function of Falkena, Wanner, Smirnov & Mo

*5. Display Data - Four plots*
5.1 Raw Position- Arm(Magenta), Eye(Blue), Gaze(Green) Head Yaw (Yellow)
5.2 Velocity
* Time = 500ms is the start of the arm motion.
Arm motion start time is determined from velocity representation, first time at which the abs value of arm velocity is 20% of the max abs arm velocity.
5.3 Candidate saccades - Plot shows in white the fast component of velocity and there is a magenta horizontal line representing the threshold. In red are saccaded above the threshold

5.4 Processed - desaccaded etc
* very good correspondance between start of movement from velocity and the by-eye start of arm movement.
* Yellow dashed line is the start of head motion, also as 10% of max head motion.
* blue dahed is start of gaze motion, set at gaze velocity great than 20 deg/sec
5.5 Note that the eye signal is taken from the 'good' eye. For EOG there is only one eye. For Elmar, the Good Eye is the eye signal from the hemifield that the speaker is in. This has the most impact on calibration, when the bad eye can saturate and throw off the scaling
5.6 In cycles mode the third plot is ommitted and the fourth plot is twice as tall as for ramps

*6. Desaccading*
% 6.1 Skirt - use 20ms 'skirt' around and saccade points 
6.2 Threshold is dynamic - difference between fast changes in velocity and gradual.
Operationally, the velocity is calculated twice, once with a window of 20 points, then again with 500 points. The absolute difference plot shows well-defined peaks when saccades occur. The default threshold is set at 30% of the maximum sacade velocity, and this is adjustable on a per trial basis. This is shown graphically in the  

7. Export
7.1 Peak pursuit, fastest segment velocity
7.2 Time to peak is the middle of the fastest segment
7.3 New metric: Peak and time to peak in 50, look at the middle 50% of the movement and check out the fastest segment there
7.4 Make sure that deleted trials works properly here
7.5 Export for cycles saves the datastructure that has columns of Speaker, eyes, desaccaded eyes with zeros for saccades, head movement
