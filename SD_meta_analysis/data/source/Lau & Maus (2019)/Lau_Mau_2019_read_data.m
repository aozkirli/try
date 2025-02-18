clc
% full path of the current script
scriptPath          = mfilename('fullpath');

% directory path
[currentDir, ~, ~]  = fileparts(scriptPath);

cd(currentDir)
fileID              = fopen('source.txt', 'r');

source              = fread(fileID, '*char')';
fclose(fileID);

fprintf(source)


% add inhouse toolbox
% this path must be replaced with the shared folder containing the
% collection of functions required
addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))

cd('Data\matfiles')
matfiles            = indir;
% The columns are as follows:
% Stimulus Value -- Stimulus Onset Timestamp -- Response Start Time -- Response End Time -- Random Response Start Position -- Sbj's Response Position -- Response Modality -- Noise Sample Start Position

clear id sizeValues
count       = 0;
tmp         = table;
for i = 1:numel(matfiles)
    filename    = matfiles(i).name;
    load(filename)
    tokens      = regexp(filename, '^([a-zA-Z]+)-(\d+)-multi-(\d+)', 'tokens');
    numid       = regexp(filename, '-(\d+)-', 'tokens');
    if size(Values,1)>100
        count               = count+1;
        id{count}           = tokens{1}{1};   

        sizeValues(count,:) = size(Values);
        theta               = Values(:,1);
        delta               = sdp_acute_ang(nbk(theta)-theta);
        resp                = Values(:,6);
        currresp            = Values(:,7);
        % Response type (cond): 1-Visual-only Responses, 2-Auditory-only Responses, 3-Visual/Audio Responses, 9-Training (Not analyzed).
        if str2double(tokens{1}{3})<=3
            cond                = str2double(tokens{1}{3}).*ones(size(resp));
        else
            cond                = nan(size(resp));
        end
        rt                  = Values(:,4)-Values(:,3);
        obsid               = repmat({id{count}},numel(cond),1);
        block               = str2double(numid{1}{1}).*ones(size(cond));
        % check visual response type
        if cond(1)==1 | currresp(1)==1
            visresp         = ones(size(resp));
        else
            visresp         = zeros(size(resp));
        end
        if size(Values,1)==104
            if cond(1)      <3
            experiment      = 1;
            else
            experiment      = 2;
            end
        else
            experiment      = 3;
        end
        experiment          = experiment.*ones(size(cond));
        resp(currresp==3)   = nan; %noresp trials, all 0
        tmp                 = vertcat(tmp,table(obsid,theta,resp,block,rt,cond,experiment,currresp,visresp,delta));
    end
end


% currresp = 1 should select only visual responses
tmp                         = tbl_subset(tmp,'visresp',1);

% Create a mapping from subject IDs to numbers
[uniqueIDs, ~, idx]         = unique(tmp.obsid);
tmp.obs                     = idx;

% check
% The idx variable now contains the numerical IDs
numIDs                      = idx;

% Display the mapping from subject IDs to numbers
mapping                     = table(uniqueIDs, (1:length(uniqueIDs))');
mapping.Properties.VariableNames = {'SubjectID', 'NumericID'};
disp(mapping);

tbl             = tmp;

fprintf('Converting to standard table format...\n\n');

tmp             = table;
tmp.obs         = tbl.obs;
tmp.theta       = tbl.theta;
tmp.resp        = tbl.resp;
tmp.error       = sdp_acute_ang(tbl.resp-tbl.theta);
tmp.rt          = tbl.rt;
tmp.block       = tbl.block;
tmp.delta       = tbl.delta; % relevant delta (previous target inducing attractive SD)
tbl             = tmp;


tbl_name      = 'Lau & Maus (2019)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\..\tables\single_experiments')
save Lau_Mau_2019.mat tbl tbl_name
writetable(tbl,'Lau_Mau_2019.csv','delimiter',';')