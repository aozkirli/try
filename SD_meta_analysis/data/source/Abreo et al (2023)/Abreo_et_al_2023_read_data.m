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


% subject list (from their osf)
subs                = {'01','02', '03', '04','05', '06','07', '08', '09', '10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36'}; % Possible subjects to exclude: 12, 20, 23, (24 - boarderline)

% note: add rotation direction variable for future subjects

%create an empty table for my data to go in
data                = table();

for s = 1:length(subs)
    
    %change directory to where indivudal block files are
    cd( [currentDir '\subjects\' subs{s}] );

    %concatonate blocks
    lf = dir('*EXPO*');
    
    %create empty matrix for all orientations to be stored
    allOri = [];
    
     %loop over blocks
    for i = 1:length(lf)

        %load each block file
        load( lf(i).name );
        
        %figure out penultimate ori
        for t = 1:dat.ntrials
            psi = 7-sum(isnan(dat.allOri(:,t)));
            penStim(t) = dat.allOri(psi,t);
        end
        
        %if adding the 3 EEG subjects, as was done in the paper
        if any(strcmp(subs{s}, {'34','35','36'}))
            dat.subject = str2num(subs{s});
        end
        
        if ~isnumeric(dat.block)
            dat.block = str2double(dat.block);
        end
        
        %convert data into a table format - NOTE: NEGATIVE ERR AND ROT DIR MEANS CCW
        temptab = table(repmat(dat.subject, dat.ntrials,1),repmat(dat.block, dat.ntrials,1), dat.randSeq', dat.numStim', dat.targOri', dat.resp', dat.conf', dat.err', dat.RTresp', dat.RTconf', dat.itiSec', dat.respOriStart', dat.wedgeStartSize', dat.reward', dat.rotDir', penStim', 'variablenames', {'snum','block','randSeq', 'numStim', 'targOri', 'resp','conf','err','RTresp','RTconf','ITIsec','respOriStart', 'wedgeStartSize', 'reward','rotDir','penStim'});
        
        %save matrix of presented orientation
        allOri = [allOri, dat.allOri];
        
%         unique(temptab.block)
        
        %add temp table to data table 
        data = [data; temptab];
        
    end
    
    %save all ori for each subject
    allOriSub{s} = allOri;   %concatonate blocks
    
end

%compute relative difference between target and penultimate stim
fin = data.targOri;
pen = data.penStim;
diffs = [rad2deg(circ_dist(deg2rad(pen), deg2rad(fin))), rad2deg(circ_dist(deg2rad(pen), deg2rad(fin-180)))];
[v,i] = min(abs(diffs),[],2);
for ii = 1:length(diffs)
    reldiff(ii) = diffs(ii,i(ii));
end
% 
%add relative difference to table
data.relDiff = reldiff';

fprintf('Converting to standard table format...\n\n');

% 
%list of subs
subs          = unique(data.snum);
% 
% store in temporary table
tmp           = table; 
tmp.delta     = round(data.relDiff); % here I should recompute and use target 
tmp.theta     = data.targOri;
tmp.resp      = data.resp;
tmp.error     = data.err;
tmp.rt        = data.RTresp;
tmp.obs       = data.snum;
tmp.block     = data.block;
tmp.cond      = data.randSeq;
% 
% select the relevant condition
tbl           = tbl_subset(tmp,'cond',1);
tbl_name      = 'Abreo et al. (2023)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\..\tables\single_experiments')
save Abreo_et_al_2023.mat tbl tbl_name
writetable(tbl,'Abreo_et_al_2023.csv','delimiter',';')