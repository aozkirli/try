Dataset for Visual serial dependence in an audiovisual stimulus

----------------------
Filename   convention
----------------------
XX-#-multi-#-DD_M_YYYY_-ZZZZZZ.mat

XX - Sbj initials
# - File number (ordered)
multi - Experiment name
# - Response type: 1-Visual-only Responses, 2-Auditory-only Responses, 3-Visual/Audio Responses, 9-Training (Not analyzed).
ZZZZZ- IGNORE / not important.

-----------------------
MATLAB DATA Explanation
-----------------------
The matrix of each file, called "Values", is 104-by-8.
Each row represent a trial.

The columns are as follows:
Stimulus Value -- Stimulus Onset Timestamp -- Response Start Time -- Response End Time -- Random Response Start Position -- Sbj's Response Position -- Response Modality -- Noise Sample Start Position

Stimulus Value refer to which stimulus is used on the current trial (both Visual and Auditory are equally mapped in a 180deg space)

Stimulus Onset Timestamp is the time the stimulus is presented

Response Start Time is the start of the response phase
Response End Time is the moment Sbj made a response by clicking the mouse

Random Response Start Position is the randomized initial start position during the response phase

Sbj's Response Position is the value selected by the sbj after clicking the mouse

Response Modality is either 1-Visual or 2-Audio
Note: In Exp 1 and 2, the values 0 because Sbj make Visual-only, or Auditory-only responses.
In Exp 3, the values are 1/2 to correspond to the response type.

Noise Sample Start Position is the randomized start position where the noise mask is played during the experiment.
The noise sample is derived from "Sounds-178sounds&Noise.mat"

