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

%% reading data using the code from osf

%% extract data
%this script needs to be run with the current directory equal to the directory containing this script

%relative paths from where this file is. MUST unzip the folders in "dependencies" first
% addpath(genpath('dependencies/')) 

%define subjects and file path
subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'}; 
fpath = [pwd '/raw_data/']; %set where the data are
clc
%create empty table
d = table();

%plot staircase for each sub?
plotStair = 0;

%starting loop over subjects
for s = 1:length(subs)
   
    % go to subject directory
    cd (strcat(fpath, subs{s}))
    
    % get file names
    files = dir('*stim*');
    
    %start block loop
    for block = 1:length(files)
        
        %load files
        load (files(block,1).name);  
%         e = getTaskParameters(myscreen,task);
        e.randVars  = task{1}.randVars;
        %NOTE: missing function in the repository (getTaskParameters), no
        %way to get RT, add fake ones (/all trials included)
        e.reactionTime = 9*ones(size(e.randVars.orientation));
%         e.reactionTime = task{1}.reactionTime;
        
%         %% staircase
%         if strcmp(task{1}.taskFilename, 'wmConEst2_thresh.m')
%             group_thresh(s,:) = [s stimulus.finalthresh];
%             
%             if plotStair == 1
%                 figure; plot(stimulus.stair.strength,'o-'); title(['subject' subs(s) 'threshold' stimulus.finalthresh])
%             end
%                   
%         end
        
        %% task blocks
        if strcmp(task{1}.taskFilename, 'wmConEst2.m')
            
            td = table();
            td.sub = cellstr(repmat(subs{s},length(e.randVars.orientation),1));
            td.ori = e.randVars.orientation';
            td.resp = e.randVars.resp';
            td.dist = e.randVars.dist';
            td.conf = e.randVars.conf';
            td.rt = e.reactionTime';
            td.posev = e.randVars.posev';
            td.delay = e.randVars.delay';
            td.block = (block-1) * ones(length(e.randVars.orientation),1);
                       
            d = [d;td]; 
        end
        
    
   
    end

end

%go back to main dir
cd ../..

% add a trial number column
d.trialnum = nan(height(d),1);
for s = 1:length(subs)
    d.trialnum(strcmp(subs{s},d.sub)) = 1:sum(strcmp(subs{s},d.sub));
    sum(strcmp(subs{s},d.sub));
end

%remove rows with no response
d(isnan(d.resp),:) = [];

%correct error in response coding
resp_fix = wrapTo360(d.resp.*2);
for i = 1:height(d)
    dist = rad2deg(circ_dist2(deg2rad(d.ori(i)),[deg2rad(resp_fix(i)),deg2rad(resp_fix(i)+180)]));
    [v,p] = min(abs(dist));
    dist_fix(i) = dist(p);
end

% insert correct values into dataframe
d.dist  = dist_fix';% NOTE THAT NEGATIVE DISTANCE VALUES INDICATES CCW ERROR
d.resp  = resp_fix;

d.ori   = mod(d.ori,180);
d.resp  = mod(d.resp,180);

% d(abs(d.dist) > 25,:) = [];

fprintf('Converting to standard table format...\n\n');
tmp             = table;

tmp.obs         = cellfun(@str2double, d.sub);
tmp.theta       = d.ori;
tmp.resp        = d.resp;
tmp.error       = sdp_acute_ang(d.resp-d.ori);
tmp.rt          = d.rt;
tmp.block       = d.block;
tmp.delta       = sdp_acute_ang(nbk(d.ori)-d.ori);
tmp.trial       = d.trialnum;
tmp.cond        = double(d.posev==1);

tbl             = tmp;

tbl_name            = 'Samaha et al. (2019)';
tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment      = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
fname         = ['Samaha_et_al_2019'];
save(fname,'tbl','tbl_name')
writetable(tbl,[fname '.csv'],'delimiter',';')

%% the following is the preprocessing script from osf, as a reference
%remove error trials
% d(abs(d.dist) > 25,:) = [];

%remove mean signed error for each sub
% fun = @(x) x-mean(x);
% tmp = varfun(fun,d,'groupingvariable','sub','inputvariable','dist');
% d.dist = tmp.Fun_dist;


% %starting loop over subjects to compute diff between current and prev ori
% prevdist= []; 
% for s = 1:length(subs)
%     td=d(strcmp(subs{s},d.sub),:);
%     
%     for i = 2:height(td)
%         dist = rad2deg(circ_dist2(deg2rad(td.ori(i)),[deg2rad(td.ori(i-1)),deg2rad(td.ori(i-1)+180)]));
%         [v,p] = min(abs(dist));
%         prevdisttemp(i) = dist(p);%starting loop over subjects
%     end
%     prevdisttemp(1)=nan;
%     prevdist=[prevdist, prevdisttemp];
%     prevdisttemp= [];
%     
% end
% d.prevdist= prevdist';

% %remove first trial of each block since no previous trial
% d.prevdist([1; diff(d.block)]~=0) = nan;
