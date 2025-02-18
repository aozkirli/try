Within the subject folders, the data of the main task are stored in files with the name data_S${SubjectID}_Task_${TaskID}_Block_${BlockID}, where ${SubjectID} indicates the subject number, ${TaskID} indicates the 2IFC task (1 = orientation discrimination; 2 = size discrimination), and ${BlockID} indicates the block number (1 to 8).
Within each of these datafiles of the main task, trial data is stored in a cell structure. The columns of each cell represent the following data:
 
1: Subject number
2: Block number
3: IFC task (1 = orientation discrimination; 2 = size discrimination)
4: Trial number
5: Orientation of inducer (ifc) grating
6: Orientation of reference (ifc) stimulus
7: Size of reference (ifc) stimulus
8: Size of inducer (ifc) grating
9: Relative orientation of inducer (ifc) grating wrt test (adjustment) grating
10: Orientation of test grating
11: Adjustment response orientation
12: Adjustment response time
13: Correctness of 2IFC response (0 = incorrect; 1 = correct)
14: 2IFC response time