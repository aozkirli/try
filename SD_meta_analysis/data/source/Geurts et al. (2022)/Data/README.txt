### List of variables in data files ###
# Note: each file contains data for a single participant #
# See associated publication for a detailed description of the (pre)processing procedures #

# Variables with one value per trial #
StimOri = True orientation of the presented stimulus (degrees)
RespOri = Orientation reported by the observer (degrees)
RespErrDebias = Debiased orientation response error (degrees; see Methods for debiasing procedure)
Conf = Confidence values reported by the observer, z-scored per session
InclTrials = Logical variable indicating for each trial whether or not it was included in further analyses (1 = included, 0 = excluded)
RunN = Run number for each trial (RunN < 100: session 1; RunN > 100: session 2)
DecOri_1000vox = Decoded stimulus orientation estimate: mean of the decoded distribution (ROI: V1-V3, max. 1000 voxels, see Methods)
DecUnc_1000vox = Decoded sensory uncertainty: variance of the decoded distribution (ROI: V1-V3, max. 1000 voxels, see Methods)
DecOri_1500vox = Decoded stimulus orientation estimate: mean of the decoded distribution (ROI: V1-V3, max. 1500 voxels, see Methods)
DecUnc_1500vox = Decoded sensory uncertainty: variance of the decoded distribution (ROI: V1-V3, max. 1500 voxels, see Methods)
DecOri_2000vox = Decoded stimulus orientation estimate: mean of the decoded distribution (ROI: V1-V3, max. 2000 voxels, see Methods)
DecUnc_2000vox = Decoded sensory uncertainty: variance of the decoded distribution (ROI: V1-V3, max. 2000 voxels, see Methods)

# Variables with one value per fMRI datapoint #
TrialN = Trial number (to link above-listed variables to fMRI timeseries)
TSStim = Logical variable indicating when stimuli were presented (1 = stimulus presentation, 0 = other)
TSResp = Logical variable indicating timing of response windows (1 = response window, 0 = other)
TSBold_dACC = Average BOLD signal in dorsal anterior cingulate cortex (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_daIns = Average BOLD signal in dorsal anterior insula (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_dPCC = Average BOLD signal in dorsal posterior cingulate cortex (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_dPGACC = Average BOLD signal in dorsal perigenual anterior cingulate cortex (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_Precuneus = Average BOLD signal in precuneus (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_rlPFC = Average BOLD signal in left rostrolatera prefrontal cortex (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_SMA = Average BOLD signal in supplementary motor area (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_V123 = Average BOLD signal in stimulus-driven voxels (max. 2000, see Methods) in V1-V3 (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)
TSBold_vPCC = Average BOLD signal in ventral posterior cingulate cortex (controlled for head motion artefacts and CSF/WM signal; see Methods for voxel selection procedure)



