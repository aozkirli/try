================================================================================================================
% Datasets for two experiments on behavioral serial dependence, using orientation adjustment tasks.
% The datasets are used and described in:
% Title:  "Serial dependence does not originate from low-level visual processing"
% Authors: Gizay Ceylan, Michael H. Herzog & David Pascucci
% Year:    2021
% Journal: Cognition
The link and DOI to the related paper will be available soon.
================================================================================================================
Header description for the tabledata.mat file of Experiment 1 and 2
Each tabledata contains all the trials of all participants, concatenated in a long format (each row is one trial).
Please refer to the original publication for details on the method and meaning of the main variables.
HEADER:
'obs'        :   ID number of participants
'gender'     :   participants' gender
'age'        :   participants' age
'delta'      :   difference (previous minus present orientation)
'resp'       :   reported orientation (degrees, 0-180, 0=vertical,90=horizontal) 
'error'      :   response error (acute angle resp minus theta)
'errorc'     :   error 'c' corrected after data cleaning
'rt'         :   adjustment time
'rtc'        :   adjustment time after data cleaning
'theta'      :   orientation presented in each trial (degrees, 0-180, 0=vertical,90=horizontal)
'block'      :   block number
'condition'  :   condition variable (Experiment 1:'lowsf','highsf','mixed'; Experiment 2: 'gab','symd','mixed')
'stimtype'   :   Experiment 1:(1=Low SF; 2=High SF); Experiment 2:(1=Gabor; 2=Symmetric dot pattern)
'outliers'   :   outlier trials from data cleaning
'stimtype_bk':   stimulus type in the preceding trial
'case'       :   combinations in 'mixed' blocks (e.g., Exp 1: 1=Low->Low SF; 2=Low->High SF; 3=High->Low SF;4=High->High SF) 
'ntrials'    :   number of 'good' trials per subject, after data cleaning