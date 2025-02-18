  %%%%%% Analysis Script for Multi %%%%%
%REQUIRES CIRCSTAT TOOLBOX FOR ANALYSIS!%
%Please note that values are all in radians. Only labels are shown in
%degrees. If you want degrees, please use circ_rad2ang( ).
%cw is jy lr ml py wk 1:8
try
    %Find data folder in current directory. If you start at the data folder, you'll get an error.
    load('ReviewBeta');
    pathname = cd; DataFolder = fullfile(pathname,'Data');
    %Load the empirical betas from the unimodal analysis
    %load('2Oct18-rawSub-fittedCurveRaw.mat'); %Check if rawSub or meanSub
    %empVbeta = fittedCurve.Paras.N1.V.Mat(:,2); %In the older script (fitted to raw Data, there is not .Mat, just .V or .A)
    %empAbeta = fittedCurve.Paras.N1.A.Mat(:,2);%empVbeta;
    empAbeta = ReviewBeta.raw;
    clear fittedCurve
    %Data folder exist in this directory
    if exist(DataFolder,'dir')
        %Let's enter that directory
        cd(DataFolder); %Remember to change meanSub or rawSub to load betas!!!
        
%         %Please note that these  are cell values (everything is strings). So use {} curly fries!
%         UserPrompt = {'Enter subject initials. Use space for separate names.', 'Please enter session numbers. Use semi-colon : for range and comma , for individual values: e.g. 1:6,8,12,31:35. Use a space for files from different subjects.', 'Enter bin degrees.', ...
%             'What is the max N-trials?', 'Do you want a plot based on all data? Y/N. Note: These are individual figures', ...
%             'Do you want to run simulated results for each subjects? Y/N'};
%         Prompt_title = 'Analysis Multi Serial Dependence'; numlines = 1; %Number of input lines, 1 is more than enough.
%         defaultValue = {'', '', '20', '2', 'N', 'N'};
%         UserInput = inputdlg(UserPrompt, Prompt_title, numlines, defaultValue); %This is the cell output
        UserInput = {'cw is jy lr ml py wk'; '1:8 1:8 1:8 1:8 1:8 1:8 1:8'; '15'; '1'}; %Original
        %UserInput = {'ak ik ms qy si yu aw ji vs js'; '1:10 1:10 1:10 1:10 1:10 1:10 1:10 1:10 1:10 1:10'; '15'; '1'}; %No response condition
        %If at least one value is left blank, we'll print this error message.
    else
        error('There is no Data folder. Analysis aborted.');
    end
    
    %% Variables
    Nbacktrials = 1;%str2num(UserInput{4}); 
    AVmatrix = []; VisualMaxResponseTime = 15; AudioMaxResponseTime = 20;
    LB = -0.5*pi; UB = 0.5*pi; Degrees = 2*pi; median = 0; bindeg = circ_ang2rad(str2double(UserInput{3})); 
    OutlierRemove = 1; wrapDeg = 0.5*pi; maxwrapDeg = 2*wrapDeg; plotTicks = circ_ang2rad(0.5); %This determines size of increment in radians
    SimsPerTrial = 2000; permnum = 5000;
    %Check if we're analyzing Data2018 or Data folder:
    whichDir = dir; 
    if size({whichDir.name},2) - 2 == 100 %This is Data2018 folder because it has 100 items
        NumOfTrials = (1:102)'; ToneSizeFilenamePosition = 8; keyword = 'multi'; %'multi-3' for original Data folder & 'multi' for Data2018 folder
        SubsetNames = {'AA'}; %{'VV','AA','AXV','VXA','NV','NA'}; %SubsetNamesExtra = {'VNV','ANA','ANV','VNA','NV','NA'};
    else
        NumOfTrials = (1:104)'; ToneSizeFilenamePosition = 9; keyword = 'multi-3';
        SubsetNames = {'AA'}; %{'VV','AA','AXV','VXA'};%,'VXV','AXA'}; %SubsetNamesExtra = {'VAV','VVA','VAA','VAAV','AVA','AVV','AVVA'};
    end
    var1 = 1; var2 = 6; var3 = 1; var4 = 1; %columns
    SbjStart = 1; minCI = .025; maxCI = .975; 
    %Ticks and Tick labels
    XTICK = [0 pi/4 pi/2]; Xlim = [0 UB]; XTICKLABEL = {'0°','45°','90°'}; %{'-90°','-45°','0°','45°','90°'};
    YlimLB = -10/180 * pi; YlimUB = 10/180 * pi;  YTICK = [YlimLB -5/180*pi 0 5/180*pi YlimUB]; YTICKLABEL = {'-10°','-5°','0°','5°','10°'};
      
    %For graph titles
    if strcmp(keyword,'multi-1') == 1; CondType = 'Visual';
    elseif strcmp(keyword,'multi-2') == 1; CondType = 'Audio';
    elseif strcmp(keyword,'multi-3') == 1; CondType = 'A/V';
    else CondType = []; 
    end
    
    %Create a variable with Permnum size Nback names and Permnum names
    for ii = 1:permnum
        Nback{ii,1} = sprintf('N%d',ii);
        permname{ii,1} = sprintf('a%d',ii);
    end
    
    %% All analysis Lazy mode.
    if strcmpi(UserInput{1},'ALL') == 1
        allnames = dir; x = 1; %counter
        for i = 1:size(allnames,1)
            %If the string matches the filename keyword, save the name in variable
            keywordcheck = strfind(allnames(i).name,keyword);
            if keywordcheck > 0 %Found something that matches keyword
                namelist{x,1} = allnames(i).name; %save it
                x = x+1;
            end
        end
        
        %Converting the saved filenames into sbj names and sessions
        x = 1; %counter
        for i = 1:size(namelist,1)
            currstring = strsplit(strrep(namelist{i},'-',' '));
            if x == 1 %Start building names
                currname{x} = currstring{1,1}; x = x+1;
            else
                %Compare names
                if strcmp(currname{x-1},currstring{1,1}) == 0 %If this is not the same name
                    currname{x} = currstring{1,1}; %save new name
                    x = x+1;
                end
            end
        end
        %Find out the number of sessions, assuming all sbjs completed same number of sessions
        sessions = size(namelist,1)/size(currname,2);
        
        %Add in the session names
        sessionsize = sprintf('1:%d',sessions);
        for i = 1:size(currname,2); currsession{i} = sessionsize; end;
        
        %Join all the names and sessions together
        UserInput{1} = strjoin(currname); UserInput{2} = strjoin(currsession);
        clear x allnames namelist currname currstring currsession sessions
    end
    
    %Recap on on columns of Raw_Matrix: 1st col stim freq, 2nd col stim onset 3rd col response start, 4th col is resp end
    %5th col start freq of resp, 6th col selected freq
    [Raw_Matrix, SbjNames, SessionInput, ToneMatrixSize] = Loadfiles(UserInput{1},UserInput{2},keyword,ToneSizeFilenamePosition);
    
    %% Remove long response time trials
    ResponseT = [Raw_Matrix(:,4) - Raw_Matrix(:,2), Raw_Matrix(:,7)]; %Show RespT and modality
    ResponseT = [(ResponseT(:,1) > VisualMaxResponseTime + ResponseT(:,2) == 1),(ResponseT(:,2) > AudioMaxResponseTime + ResponseT(:,2) == 2)];
    Raw_Matrix(ResponseT(:,1),6) = NaN; Raw_Matrix(ResponseT(:,2),6) = NaN; %Convert the response for that trial to NaN
    if size(NumOfTrials,1) == 102 %Data2018
       Raw_Matrix(Raw_Matrix(:,7) == 3,6) = NaN; %We cannot analyze no response trials 
    end
    
    %% Mean Subtract in RAW Response
    temp_Raw = Raw_Matrix(:,[1 6 7]); %Stim, Resp, Modality currTrial    
    for currStim = min(temp_Raw(:,1)):max(temp_Raw(:,1)) %all 178 stimulus
        for currMod = 1:2 %for each modality
            %Get the matrix for the current stimulus number
            temp_Mat = temp_Raw(((temp_Raw(:,1) == currStim) + (temp_Raw(:,3) == currMod) == 2),[1 2]);
            %Mean values
            temp_Mat(:,2) = circ_ang2rad(temp_Mat(:,2).*2); %Change to rad
            temp_Mean = circ_mean(temp_Mat(:,2)); %Mean error in response
            Unity_to_Mean_Distance =  temp_Mean - circ_ang2rad(currStim*2); %Error from Unity 21
            temp_Mat(:,2) = temp_Mat(:,2) - Unity_to_Mean_Distance; %Mean subtract using distance from unity
            temp_Mat(:,2) = round(circ_rad2ang(temp_Mat(:,2))./2); %Convert to deg,rounded off
            %Wrapping
            temp_Mat(temp_Mat(:,2) < 0,2) = temp_Mat(temp_Mat(:,2) < 0, 2) + 180;
            temp_Mat(temp_Mat(:,2) > 179,2) = temp_Mat(temp_Mat(:,2) > 179, 2 ) - 180;
            temp_Raw(((temp_Raw(:,1) == currStim) + (temp_Raw(:,3) == currMod) == 2),[1 2]) = temp_Mat(:,[1 2]);
        end
    end
    
    %Replace the raw response values after cleanup
    Raw_Matrix(:,6) = temp_Raw(:,2);
    
    %% Form response errors and such
    GroupAnalysisMatrix = Raw_Matrix(:,var2) - Raw_Matrix(:,var1);
        
%     %% For N+1 analysis
%     GroupNOne = GroupAnalysisMatrix;
%     GroupNOne(:,2) = circshift(Raw_Matrix(:,var3),-1) - Raw_Matrix(:,var4); %Only N+1
%     GroupNOne(:,[3 4]) = [Raw_Matrix(:,7), circshift(Raw_Matrix(:,7),-1)]; %Column 3 = mod current trial, Col 4 = mod next trial
%     GroupNOne(:,1:2) = circ_ang2rad(GroupNOne(:,1:2).*(360/ToneMatrixSize(1))); %Convert to Deg
%     for ii = 1:5 %Wrap response error and N+1
%         temp_mat1 = GroupNOne(:,1); temp_mat2 = GroupNOne(:,2);
%         temp_mat1(temp_mat1 > wrapDeg) = temp_mat1(temp_mat1 > wrapDeg) - maxwrapDeg;
%         temp_mat2(temp_mat2 > wrapDeg) = temp_mat2(temp_mat2 > wrapDeg) - maxwrapDeg;
%         temp_mat1(temp_mat1 < -wrapDeg) = temp_mat1(temp_mat1 < -wrapDeg) + maxwrapDeg;
%         temp_mat2(temp_mat2 < -wrapDeg) = temp_mat2(temp_mat2 < -wrapDeg) + maxwrapDeg;
%         GroupNOne(:,[1 2]) = [temp_mat1, temp_mat2]; clear temp_mat1 temp_mat2;
%     end
    %End of N+1 analysis
    
    %Group Nback trials in columns
    for aa = 1:Nbacktrials
        GroupAnalysisMatrix(:,aa+1) = circshift(Raw_Matrix(:,var3),aa) - Raw_Matrix(:,var4);
    end
    %Add modality information
    for aa = 1:Nbacktrials+1 %1 more column because first modality column is for current trial
        GroupAnalysisMatrix(:,Nbacktrials+1+aa) = circshift(Raw_Matrix(:,7),aa-1);  %Since a/v analyzes only n-1, 3rd column of matrix is modality info
    end
    %Convert to degrees
    GroupAnalysisMatrix(:,1:Nbacktrials+1) = circ_ang2rad(GroupAnalysisMatrix(:,1:Nbacktrials+1).*(360/ToneMatrixSize(1)));
    
    %% Reviewer suggested analysis (1-back only)
    temp = GroupAnalysisMatrix(:,[1 2]);
    %temp(temp(:,1) <= 0 & temp(:,2) <= 0,:) = temp(temp(:,1) <= 0 & temp(:,2) <= 0,:) * -1;
    for t = 1:size(temp,1)
        if temp(t,1) < 0 && temp(t,2) < 0
            temp(t,:) = temp(t,:) * -1;
        elseif temp(t,1) > 0 && temp(t,2) < 0 %Convert this to positive values
            temp(t,:) = temp(t,:) * -1;
        end
    end
   
    %% Wrap data at wrapDeg
    for ff = 1:Nbacktrials+1
        Temp_currData = GroupAnalysisMatrix(:,ff);
        for loopTheWrap = 1:5 %Loop this a few times, more thorough
            Temp_currData(Temp_currData > wrapDeg) = Temp_currData(Temp_currData > wrapDeg) - maxwrapDeg;
            Temp_currData(Temp_currData < -wrapDeg) = Temp_currData(Temp_currData < -wrapDeg) + maxwrapDeg;
            %Review comments -- gotta wrap these too
            temp(temp > wrapDeg) = temp(temp > wrapDeg) - maxwrapDeg;
            temp(temp < -wrapDeg) = temp(temp < -wrapDeg) + maxwrapDeg;
        end
        GroupAnalysisMatrix(:,ff) = Temp_currData;
        clear Temp_currData
    end
    
    %Subtract the mean of the new calculations
    temp(:,1) = temp(:,1) - circ_nanmean(temp(:,1));
    GroupAnalysisMatrix(:,[1 2]) = temp;
    
%     %Below section commented out if removing mean errors from raw data
%     %Group error and SD
%     GroupMatrix.VMeanError = circ_nanmean(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 1,1));
%     GroupMatrix.AMeanError = circ_nanmean(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 2,1));
%     GroupMatrix.VSD = circ_nanstd(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 1,1));
%     GroupMatrix.ASD = circ_nanstd(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 2,1));
%     
%     %Remove other outliers and Mean error subtraction
%     for eachtrial = 1:size(GroupAnalysisMatrix,1) %For each iteration 
%         if GroupAnalysisMatrix(eachtrial,Nbacktrials+2) == 1 %For visual trials
%             %Mean error subtraction
%             GroupAnalysisMatrix(eachtrial,1) = GroupAnalysisMatrix(eachtrial,1) - GroupMatrix.VMeanError;
% %             if GroupAnalysisMatrix(eachtrial,1) >= circ_ang2rad(45) || GroupAnalysisMatrix(eachtrial,1) <= circ_ang2rad(-45)
% %             %if GroupAnalysisMatrix(eachtrial,1) >= GroupMatrix.VSD*3 || GroupAnalysisMatrix(eachtrial,1) <= GroupMatrix.VSD*-3
% %                  GroupAnalysisMatrix(eachtrial,1) = NaN; %Set that value as NaN
% %             end
%         else %For auditory trials
%             %Mean error subtraction
%             GroupAnalysisMatrix(eachtrial,1) = GroupAnalysisMatrix(eachtrial,1) - GroupMatrix.AMeanError;
% %             if GroupAnalysisMatrix(eachtrial,1) >= circ_ang2rad(45) || GroupAnalysisMatrix(eachtrial,1) <= circ_ang2rad(-45)
% %             %if GroupAnalysisMatrix(eachtrial,1) >= GroupMatrix.ASD*3 || GroupAnalysisMatrix(eachtrial,1) <= GroupMatrix.ASD*-3
% %                 GroupAnalysisMatrix(eachtrial,1) = NaN; %Set that value as NaN
% %             end
%         end
%     end
%     
%     %Group error and SD post subtraction
%     GroupMatrix.VPostMeanError = circ_nanmean(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 1,1));
%     GroupMatrix.APostMeanError = circ_nanmean(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 2,1));
% %     GroupMatrix.VPostSD = circ_nanstd(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 1,1));
% %     GroupMatrix.APostSD = circ_nanstd(GroupAnalysisMatrix(GroupAnalysisMatrix(:,Nbacktrials+2) == 2,1));
%     
%     %Save the Group analysis matrix
%     GroupMatrix.Analysis_Matrix = GroupAnalysisMatrix;
    
    %Sort data into smaller
    for currNback = 1:Nbacktrials
        clear currData
        %List of trials to skip
        SkippedTrials(currNback,:) = currNback:size(NumOfTrials,1):size(GroupAnalysisMatrix,1);
        TrialsToNan = SkippedTrials(1:currNback,:); TrialsToNan = TrialsToNan(:);
        currData = GroupAnalysisMatrix;
        %Append trials to skip (Nbacks), 4th column 1 = skip trial
        currData(TrialsToNan,size(currData,2)+currNback) = 1;
        
        %Sorting modalities, always comparing the first modality column (current trial)
        %with modalities in the next column (prev trials)
        %VV = currData(sum(currData(:,(Nbacktrials+2):(Nbacktrials+2+currNback)) == 1,2) == (currNback+1),[1 currNback+1 Nbacktrials+2 end]);
        AA = currData(sum(currData(:,(Nbacktrials+2):(Nbacktrials+2+currNback)) == 2,2) == (currNback+1),[1 currNback+1 Nbacktrials+2 end]);
        %sum of current modality and N-back modality = 2
        %VXA = currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 1) == 2),[1 currNback+1 (Nbacktrials+2):end]);
        %AXV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 2) == 2),[1 currNback+1 (Nbacktrials+2):end]);
        %NV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 3) == 2),[1 currNback+1 (Nbacktrials+2):end]);
        %NA =  currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 3) == 2),[1 currNback+1 (Nbacktrials+2):end]);
        %VXV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 1) == 2),[1 currNback+1 Nbacktrials+2 end]);
        %AXA = currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 2) == 2),[1 currNback+1 Nbacktrials+2 end]);
%         if size(NumOfTrials,1) == 102 %This is Data2018
%             if currNback == 2
%                 VNV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 ANA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 VNA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 ANV = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 NA = [ANA;VNA];
%                 NV = [VNV;ANV];
%             end
%         end
%         %Diffmodal calculation
%         if currNback == 1
%             diffModal = currData(currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+2+currNback)),[1 currNback+1 end]);
%         elseif currNback == 2
%             diffModal = currData((currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+1+currNback))) + (currData(:,(Nbacktrials+2+currNback)) == currData(:,(Nbacktrials+2))) == 2,[1 currNback+1 Nbacktrials+2 end]);
%         elseif currNback == 3
%             diffModal = currData((currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+currNback))) ...
%                 + (currData(:,(Nbacktrials+2)) == currData(:,(Nbacktrials+1+currNback))) ...
%                 + (currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+2+currNback))) == 3, [1 currNback+1 Nbacktrials+2 end]);
%         end
%         %Ignoreall modality
%         Ignoreall = currData(:,[1 currNback+1 Nbacktrials+2 end]);
        
%         %Plot the following and store them at N1, but only when Nback is up to N3
%         %Since this is not a conventional Nback analysis, don't need to remove trials
%         if currNback == 2 && (Nbacktrials >= 3)
%             VAV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%             VVA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2:end]);
%             VAA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2:end]);
%             AVA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%             AVV = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%             %Remove skipped trials
%             VAV(VAV(:,end) == 1,:) = [];
%             VVA(VVA(:,end) == 1,:) = [];
%             VAA(VAA(:,end) == 1,:) = []; 
%             AVA(AVA(:,end) == 1,:) = [];
%             AVV(AVV(:,end) == 1,:) = [];  
%             %Save data
%             GroupMatrix.N1.VAV = VAV(:,1:3);
%             GroupMatrix.N1.VVA = VVA(:,1:3);
%             GroupMatrix.N1.VAA = VAA(:,1:3);        
%             GroupMatrix.N1.AVA = AVA(:,1:3);
%             GroupMatrix.N1.AVV = AVV(:,1:3);
%         elseif currNback == 3 && (Nbacktrials >= 3)
%             VAAV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 1) == 4),[1 currNback+1 Nbacktrials+2 end]);
%             AVVA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                 (currData(:,(Nbacktrials+currNback)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+1)) == 1) + ...
%                 (currData(:,(Nbacktrials+currNback+2)) == 2) == 4),[1 currNback+1 Nbacktrials+2 end]);
%             %Remove Skipped trials and save data.
%             VAAV(VAAV(:,end) == 1,:) = [];
%             AVVA(AVVA(:,end) == 1,:) = [];
%             GroupMatrix.N1.VAAV = VAAV(:,1:3);
%             GroupMatrix.N1.AVVA = AVVA(:,1:3);     
%         end
          
        %Removed skipped trials
        %VV(VV(:,end) == 1,:) = [];
        AA(AA(:,end) == 1,:) = [];
        %VXA(VXA(:,end) == 1,:) = [];
        %AXV(AXV(:,end) == 1,:) = [];
        %NA(NA(:,end) == 1,:) = [];
        %NV(NV(:,end) == 1,:) = [];
        %VXV(VXV(:,end) == 1,:) = [];
        %AXA(AXA(:,end) == 1,:) = [];
        %diffModal(diffModal(:,end) == 1,:) = [];
        %Ignoreall(Ignoreall(:,end) == 1,:) = [];
%         if size(NumOfTrials,1) == 102 
%             if currNback == 2 %This is Data2018
%                 VNV(VNV(:,end) == 1,:) = [];
%                 ANA(ANA(:,end) == 1,:) = [];
%                 VNA(VNA(:,end) == 1,:) = [];
%                 ANV(ANV(:,end) == 1,:) = [];
%                 NV(NV(:,end) == 1,:) = [];
%                 NA(NA(:,end) == 1,:) = [];
%             end
%         end        

        %Save data
        %GroupMatrix.(Nback{currNback}).VV = VV(:,1:3);
        GroupMatrix.(Nback{currNback}).AA = AA(:,1:3);
        %GroupMatrix.(Nback{currNback}).VXA = VXA(:,1:3);
        %GroupMatrix.(Nback{currNback}).AXV = AXV(:,1:3);
        %GroupMatrix.(Nback{currNback}).NV = NV(:,1:3);
        %GroupMatrix.(Nback{currNback}).NA = NA(:,1:3);
        %GroupMatrix.(Nback{currNback}).VXV = VXV(:,1:3);
        %GroupMatrix.(Nback{currNback}).AXA = AXA(:,1:3);
        %GroupMatrix.(Nback{currNback}).DiffModality = diffModal(:,1:3);
        %GroupMatrix.(Nback{currNback}).Ignoreall = Ignoreall(:,1:3);
%         if size(NumOfTrials,1) == 102
%             if currNback == 2 %This is Data2018
%                 GroupMatrix.(Nback{currNback}).ANA = ANA(:,1:3);
%                 GroupMatrix.(Nback{currNback}).VNV = VNV(:,1:3);
%                 GroupMatrix.(Nback{currNback}).ANV = ANV(:,1:3);
%                 GroupMatrix.(Nback{currNback}).VNA = VNA(:,1:3);
%                 GroupMatrix.(Nback{currNback}).NA = NA(:,1:3);
%                 GroupMatrix.(Nback{currNback}).NV = NV(:,1:3);
%             end
%         end
        clear VV AA VXA AXV AXA VXV currData diffModal Ignoreall AVA AVVA VAV VAAV AVV VVA VAA ANA VNV ANV VNA NA NV
    end
    clear SkippedTrials
    
    %% Permutate Group matrix
    Perm_GroupRawMatrix = PermResponse(Raw_Matrix,size(Raw_Matrix,1),permnum);
    
    %Build permutation matrix
    currNameP = fieldnames(Perm_GroupRawMatrix); %Getting all the fieldnames within the structure
    for qq = 1:length(fieldnames(Perm_GroupRawMatrix))
        %Current matrix
        currData = Perm_GroupRawMatrix.(currNameP{qq});
        %Response error
        temp_currData = currData(:,var2) - currData(:,var1);
        %Nback trials stored in columns
        for rr = 1:Nbacktrials
            temp_currData(:,rr+1) = circshift(currData(:,var3),rr) - currData(:,var4);
        end
        %Include modality for each trials for A/V condition toward the last column
        if strcmp(keyword,'multi-3') == 1 || size(NumOfTrials,1) == 102
            for rr = 1:Nbacktrials+1
                temp_currData(:,Nbacktrials+1+rr) = circshift(currData(:,7),rr-1);
            end
        end
        
        %Convert Degrees
        temp_currData(:,1:(Nbacktrials+1)) = circ_ang2rad(temp_currData(:,1:(Nbacktrials+1))).*(360/ToneMatrixSize(1));
        
        %% Reviewer suggested analysis (1-back only)
        tempPerm = temp_currData(:,[1 2]);
        %temp(temp(:,1) <= 0 & temp(:,2) <= 0,:) = temp(temp(:,1) <= 0 & temp(:,2) <= 0,:) * -1;
        for t = 1:size(tempPerm,1)
            if tempPerm(t,1) < 0 && tempPerm(t,2) < 0
                tempPerm(t,:) = tempPerm(t,:) * -1;
            elseif tempPerm(t,1) > 0 && tempPerm(t,2) < 0 %Convert this to positive values
                tempPerm(t,:) = tempPerm(t,:) * -1;
            end
        end
         
        %Wrap matrix
        for ff = 1:Nbacktrials+1
            WrapMatrix = temp_currData(:,ff);
            for loopTheWrap = 1:5 %Loop this a few times, more thorough
                WrapMatrix(WrapMatrix > wrapDeg) = WrapMatrix(WrapMatrix > wrapDeg) - maxwrapDeg;
                WrapMatrix(WrapMatrix < -wrapDeg) = WrapMatrix(WrapMatrix < -wrapDeg) + maxwrapDeg;
                tempPerm(tempPerm > wrapDeg) = tempPerm(tempPerm > wrapDeg) - maxwrapDeg;
                tempPerm(tempPerm < -wrapDeg) = tempPerm(tempPerm < -wrapDeg) + maxwrapDeg;
            end
            temp_currData(:,ff) = WrapMatrix;
            clear WrapMatrix
        end 
        
        %Subtract the mean of the new calculations
        temp(:,1) = temp(:,1) - circ_nanmean(temp(:,1));
        Analy.AAperm.Review = temp;
        temp_currData(:,[1 2]) = tempPerm;
        
%         %% Remove permutate outlier 3sd and mean subtraction
%         for eachtrial = 1:size(temp_currData,1) %For each iteration of current permutation
% %            if temp_currData(eachtrial,Nbacktrials+2) == 1 %For visual trials
% %                %Outliers
% %                 if temp_currData(eachtrial,1) >= GroupMatrix.VSD*3 || temp_currData(eachtrial,1) <= GroupMatrix.VSD*-3
% %                      temp_currData(eachtrial,1) = NaN; %Set that value as NaN
% %                 else %Mean error subtraction, subtract 2 mean errors because that is what we did in original analysis
% %                    temp_currData(eachtrial,1) = temp_currData(eachtrial,1) - GroupMatrix.VMeanError - GroupMatrix.VPostMeanError;
% %                 end
% %            else %For auditory trials
% %                %Outliers
% %                 if temp_currData(eachtrial,1) >= GroupMatrix.ASD*3 || temp_currData(eachtrial,1) <= GroupMatrix.ASD*-3
% %                     temp_currData(eachtrial,1) = NaN; %Set that value as NaN
% %                 else %Mean error subtraction, subtract 2 mean errors because that is what we did in original analysis
%                     temp_currData(eachtrial,1) = temp_currData(eachtrial,1) - GroupMatrix.AMeanError - GroupMatrix.APostMeanError;
% %                 end
% %            end
%         end
% 
%         %Saving calculations
%         GroupPermMatrix.(currNameP{qq}) = temp_currData;
        
        %% Permutate a/v
        %Splitting data into smaller chunks if from Multi-3
        for currNback = 1:Nbacktrials
            %matrix with errors, Nbacks, modality
            currData = temp_currData;
            %List of trials to skip
            SkippedTrials(currNback,:) = currNback:size(NumOfTrials,1):size(GroupAnalysisMatrix,1);
            TrialsToNan = SkippedTrials(1:currNback,:); TrialsToNan = TrialsToNan(:);
            %Append trials to skip (Nbacks), 4th column 1 = skip trial
            currData(TrialsToNan,size(currData,2)+currNback) = 1;
            
            %Sorting modalities, always comparing the first modality column (current trial)
            %with modalities in the next column (prev trials)
            %VV = currData(sum(currData(:,(Nbacktrials+2):(Nbacktrials+2+currNback)) == 1,2) == (currNback+1),[1 currNback+1 Nbacktrials+2 end]);
            AA = currData(sum(currData(:,(Nbacktrials+2):(Nbacktrials+2+currNback)) == 2,2) == (currNback+1),[1 currNback+1 Nbacktrials+2 end]);
            %sum of current modality and N-back modality = 2
            %VXA = currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 1) == 2),[1 currNback+1 Nbacktrials+2 end]);
            %AXV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 2) == 2),[1 currNback+1 Nbacktrials+2 end]);
            %VXV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 1) == 2),[1 currNback+1 Nbacktrials+2 end]);
            %AXA = currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 2) == 2),[1 currNback+1 Nbacktrials+2 end]);
            %NV = currData(((currData(:,(Nbacktrials+2)) == 1) + (currData(:,(Nbacktrials+2+currNback)) == 3) == 2),[1 currNback+1 (Nbacktrials+2):end]);
            %NA =  currData(((currData(:,(Nbacktrials+2)) == 2) + (currData(:,(Nbacktrials+2+currNback)) == 3) == 2),[1 currNback+1 (Nbacktrials+2):end]);            
%             if size(NumOfTrials,1) == 102 %This is Data2018
%                 if currNback == 2
%                     VNV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                         (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                         (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                     ANA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                         (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                         (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                     VNA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                         (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                         (currData(:,(Nbacktrials+currNback+2)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                     ANV = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                         (currData(:,(Nbacktrials+currNback+1)) == 3) + ...
%                         (currData(:,(Nbacktrials+currNback+2)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                     NA = [ANA;VNA];
%                     NV = [VNV;ANV];
%                 end
%             end

%             %Diffmodal calculation, up to Nbacktrials = 3
%             if currNback == 1
%                 diffModal = currData(currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+2+currNback)),[1 currNback+1 end]);
%             elseif currNback == 2
%                 diffModal = currData((currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+1+currNback))) + (currData(:,(Nbacktrials+2+currNback)) == currData(:,(Nbacktrials+2))) == 2,[1 currNback+1 Nbacktrials+2 end]);
%             elseif currNback == 3
%                 diffModal = currData((currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+currNback))) ...
%                     + (currData(:,(Nbacktrials+2)) == currData(:,(Nbacktrials+1+currNback))) ...
%                     + (currData(:,(Nbacktrials+2)) ~= currData(:,(Nbacktrials+2+currNback))) == 3, [1 currNback+1 Nbacktrials+2 end]);
%             end
%             
%             %Ignoreall modality
%             Ignoreall = currData(:,[1 currNback+1 Nbacktrials+2 end]);
            
%             %Plot the following and store them at N1, but only when Nback is up to N3
%             %Since this is not a conventional Nback analysis, don't need to remove trials
%             if currNback == 2 && (Nbacktrials >= 3)
%                 VAV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+1+currNback)) == 2) + ...
%                     (currData(:,(Nbacktrials+2+currNback)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 VVA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+1+currNback)) == 1) + ...
%                     (currData(:,(Nbacktrials+2+currNback)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 VAA = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+1+currNback)) == 2) + ...
%                     (currData(:,(Nbacktrials+2+currNback)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 AVA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                     (currData(:,(Nbacktrials+1+currNback)) == 1) + ...
%                     (currData(:,(Nbacktrials+2+currNback)) == 2) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 AVV = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                     (currData(:,(Nbacktrials+1+currNback)) == 1) + ...
%                     (currData(:,(Nbacktrials+2+currNback)) == 1) == 3),[1 currNback+1 Nbacktrials+2 end]);
%                 %Remove skipped trials
%                 VAV(VAV(:,end) == 1,:) = [];
%                 VVA(VVA(:,end) == 1,:) = [];
%                 VAA(VAA(:,end) == 1,:) = [];
%                 AVA(AVA(:,end) == 1,:) = [];
%                 AVV(AVV(:,end) == 1,:) = [];
%                 %Save data
%                 GroupMatrix.N1.Perm.VAV = VAV(:,1:3);
%                 GroupMatrix.N1.Perm.VVA = VVA(:,1:3);
%                 GroupMatrix.N1.Perm.VAA = VAA(:,1:3);
%                 GroupMatrix.N1.Perm.AVA = AVA(:,1:3);
%                 GroupMatrix.N1.Perm.AVV = AVV(:,1:3);
%             elseif currNback == 3 && (Nbacktrials >= 3)
%                 VAAV = currData(((currData(:,(Nbacktrials+2)) == 1) + ...
%                     (currData(:,(Nbacktrials+currNback)) == 2) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 2) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 1) == 4),[1 currNback+1 Nbacktrials+2 end]);
%                 AVVA = currData(((currData(:,(Nbacktrials+2)) == 2) + ...
%                     (currData(:,(Nbacktrials+currNback)) == 1) + ...
%                     (currData(:,(Nbacktrials+currNback+1)) == 1) + ...
%                     (currData(:,(Nbacktrials+currNback+2)) == 2) == 4),[1 currNback+1 Nbacktrials+2 end]);
%                 %Remove skipped trials and save data
%                 VAAV(VAAV(:,end) == 1,:) = [];
%                 AVVA(AVVA(:,end) == 1,:) = [];
%                 GroupMatrix.N1.Perm.VAAV = VAAV(:,1:3);
%                 GroupMatrix.N1.Perm.AVVA = AVVA(:,1:3);
%             end
               
            %Removed skipped trials
            %VV(VV(:,end) == 1,:) = [];
            AA(AA(:,end) == 1,:) = [];
            %VXA(VXA(:,end) == 1,:) = [];
            %AXV(AXV(:,end) == 1,:) = [];
            %NV(NV(:,end) == 1,:) = [];
            %NA(NA(:,end) == 1,:) = [];
            %VXV(VXV(:,end) == 1,:) = [];
            %AXA(AXA(:,end) == 1,:) = [];
            %diffModal(diffModal(:,end) == 1,:) = [];
            %Ignoreall(Ignoreall(:,end) == 1,:) = [];
%             if size(NumOfTrials,1) == 102
%                 if currNback == 2 %This is Data2018
%                     VNV(VNV(:,end) == 1,:) = [];
%                     ANA(ANA(:,end) == 1,:) = [];
%                     VNA(VNA(:,end) == 1,:) = [];
%                     ANV(ANV(:,end) == 1,:) = [];
%                     NV(NV(:,end) == 1,:) = [];
%                     NA(NA(:,end) == 1,:) = [];
%                 end
%             end

            %Save data
            %GroupMatrix.(Nback{currNback}).Perm.VV = VV(:,1:3);
            GroupMatrix.(Nback{currNback}).Perm.AA = AA(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.VXA = VXA(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.AXV = AXV(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.NV = NV(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.NA = NA(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.VXV = VXV(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.AXA = AXA(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.(currName{qq}).DiffModality = diffModal(:,1:3);
            %GroupMatrix.(Nback{currNback}).Perm.(currName{qq}).Ignoreall = Ignoreall(:,1:3);
%             if size(NumOfTrials,1) == 102
%                 if currNback == 2 %This is Data2018
%                     GroupMatrix.(Nback{currNback}).Perm.ANA = ANA(:,1:3);
%                     GroupMatrix.(Nback{currNback}).Perm.VNV = VNV(:,1:3);
%                     GroupMatrix.(Nback{currNback}).Perm.ANV = ANV(:,1:3);
%                     GroupMatrix.(Nback{currNback}).Perm.VNA = VNA(:,1:3);
%                     GroupMatrix.(Nback{currNback}).Perm.NA = NA(:,1:3);
%                     GroupMatrix.(Nback{currNback}).Perm.NV = NV(:,1:3);
%                 end
%             end
            clear VV AA VXA AXV AXA VXV currData diffModal Ignoreall AVA AVVA VAV VAAV AVV VVA VAA ANA VNV ANV VNA NA NV
            
%             if qq == 2500 || qq == 5000 || qq == 7500 || qq == 10000 %Plot the shuffled stuff
%                 figure;
%                 currData = GroupMatrix.(Nback{ii-2}).Perm.AA;
%                 %scatter(currData(:,2),currData(:,1),18,'filled','o','MarkerFaceColor',[1 0 0]);
%                 hold on;
%                 TheBins = runAvg(currData(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
%                 
%                 plot(Xlim,[0 0], 'k', 'LineWidth', 0.25); %Draws a line on Y axis.
%                 RunningAvg = plot(LB:plotTicks:UB, TheBins, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);%Running average plot
%                 
%                 %Axis attribtues.
%                 AnalysisPlot = gca; %this will be used to set attributes of axis
%                 set(AnalysisPlot, 'XTick', XTICK); set(AnalysisPlot, 'Xlim', Xlim); set(AnalysisPlot, 'XTickLabel', XTICKLABEL);
%                 set(AnalysisPlot, 'YTick', YTICK); set(AnalysisPlot, 'Ylim', [YlimLB YlimUB]); set(AnalysisPlot, 'YTickLabel',YTICKLABEL);
%                 xlabel('Previous Trial -- Current Trial (Deg)'); ylabel('Response Error (Deg)');
%                 plottitle = sprintf('Shuffle Number %d %s Trials', qq, CondType); title(plottitle);
%                 %legend([DoG(2) RunningAvg],'Fitted curve','Average','Location','northeast','Orientation','vertical'); legend('boxoff'); %Average is the name of the RunningAvg line
%                 hold off; clear currData RunningAvg TheBins DoGestimates DoGplot
%             end
            
        end
        
        %% Confidence Interval
        currNameCI = fieldnames(GroupMatrix.N1.Perm);
        %Fit data to all the permutations to get estimated alpha
        for ii = 2:Nbacktrials+1
            %For each comparisons
            for currSubset = 1:size(SubsetNames,2)
                %For each permutation
                currField = GroupMatrix.(Nback{ii-1}).Perm.(SubsetNames{currSubset});
                %Get average
                PermBins = runAvg(currField(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
                %Fixing to empirical Values for curve fitting
                %if currSubset == 2 || currSubset == 4 || currSubset == 6 %Fix beta to empirical audio beta
                    [fitresult, gof, fitCoeff] = fitDoGv2(0:plotTicks:UB,PermBins(181:end)',empAbeta); %currField(:,2),currField(:,1)
                %else %Visual and other combined trials
                %    [fitresult, gof, fitCoeff] = fitDoGv2(LB:plotTicks:UB,PermBins',empVbeta);
                %end
                fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNames{currSubset}).Matrix(qq,:) = fitCoeff;
            end
            
%             %For each new comparisons
%             if ii == 2
%                 for currSubset = 1:size(SubsetNamesExtra,2)
%                     %For each permutation
%                     %for currPerm = 1:permnum
%                         %                     if currSubset == 1 %VXVCombine
%                         %                         currField = [GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VV(:,[1 2]); ...
%                         %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VAV(:,[1 2]); ...
%                         %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VAAV(:,[1 2])];
%                         %                     elseif currSubset == 4 %AXACombine
%                         %                         currField = [GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AA(:,[1 2]); ...
%                         %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AVA(:,[1 2]); ...
%                         %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AVVA(:,[1 2])];
%                         %                     else
%                         if size(NumOfTrials,1) == 102
%                             currField = GroupMatrix.N2.Perm.(SubsetNamesExtra{currSubset});
%                         else
%                             currField = GroupMatrix.(Nback{ii-1}).Perm.(SubsetNamesExtra{currSubset});
%                         end
%                         %Fixing to empirical Values for curve fitting
%                         if currSubset == 2 || currSubset == 4 || currSubset == 6 %Fix beta to empirical audio beta
%                             [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empAbeta);
%                         else %Visual and other combined trials
%                             [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empVbeta);
%                         end
%                         fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNamesExtra{currSubset}).Matrix(qq,:) = fitCoeff;
%                     %end
%                 end
%             end
            
        end
        GroupMatrix.(Nback{currNback}).Perm = []; %Reduce file size
        clear temp_currData currData %GroupMatrix
    end

    %% Plots
    for ii = 2:Nbacktrials+1
        for currSubset = 1:size(SubsetNames,2)
            % subplot(2,4,currSubset);
            figure;
            %Plot
            currData = GroupMatrix.(Nback{ii-1}).(SubsetNames{currSubset});
            %scatter(currData(:,2),currData(:,1),18,'filled','o','MarkerFaceColor',[1 0 0]); 
            hold on;
            [TheBins, TheSEM] = runAvg(currData(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
            
            %Fixing to empirical Values for curve fitting
            %if currSubset == 2 || currSubset == 4 || currSubset == 6 %Fix beta to empirical audio beta
                [fitresult, gof, fitCoeff] = fitDoGv2(0:plotTicks:UB,TheBins(181:end)',empAbeta); %currData(:,2),currData(:,1)
            %else %Visual and other combined trials
            %    [fitresult, gof, fitCoeff] = fitDoGv2(LB:plotTicks:UB,TheBins',empVbeta);
            %end
            
            shadedErrorBar(0:plotTicks:UB,TheBins(181:end),TheSEM(181:end)); %Plot SEM shaded region
            DoG = plot(fitresult,'b',currData(:,2),currData(:,1)); DoG(1).LineWidth = 1.5; DoG(1).delete %This removes the blue dots. It's a duplicate of the scatter.
            
            %fitted curve parameters each sbj
            fittedCurve.Paras.(Nback{ii-1}).(SubsetNames{currSubset}).Mat = fitCoeff;
            fittedCurve.Paras.(Nback{ii-1}).(SubsetNames{currSubset}).gof = gof;
            
            plot(Xlim,[0 0], 'k', 'LineWidth', 0.25); %Draws a line on Y axis.
            RunningAvg = plot(0:plotTicks:UB, TheBins(181:end), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);%Running average plot
            
            %Axis attribtues.
            AnalysisPlot = gca; %this will be used to set attributes of axis
            set(AnalysisPlot, 'XTick', XTICK); set(AnalysisPlot, 'Xlim', Xlim); set(AnalysisPlot, 'XTickLabel', XTICKLABEL);
            set(AnalysisPlot, 'YTick', YTICK); set(AnalysisPlot, 'Ylim', [YlimLB YlimUB]); set(AnalysisPlot, 'YTickLabel',YTICKLABEL);
            xlabel('Previous Trial -- Current Trial (Deg)'); ylabel('Response Error (Deg)');
            legend([DoG(2) RunningAvg],'Fitted curve','Average','Location','northeast','Orientation','vertical'); legend('boxoff'); %Average is the name of the RunningAvg line
            hold off; clear currData RunningAvg TheBins DoGestimates DoGplot
        end
        
%         %Plot AVA, AVVA, VAV, VAAV etc
%         if ii == 2
%             
%             for currSubset = 1:size(SubsetNamesExtra,2)
%                 %subplot(2,3,currSubset);
%                 figure;
%                 %Plot
%                 %if currSubset == 1 %Plot VXVCombine
%                 %    currData = [GroupMatrix.(Nback{ii-1}).VV; ...
%                 %        GroupMatrix.(Nback{ii-1}).VAV; ...
%                 %        GroupMatrix.(Nback{ii-1}).VAAV];
%                 %elseif currSubset == 4 %Plot AXACombine
%                 %    currData = [GroupMatrix.(Nback{ii-1}).AA; ...
%                 %        GroupMatrix.(Nback{ii-1}).AVA; ...
%                 %        GroupMatrix.(Nback{ii-1}).AVVA];
%                 %else
%                 if size(NumOfTrials,1) == 102
%                     currData = GroupMatrix.N2.(SubsetNamesExtra{currSubset});
%                 else
%                     currData = GroupMatrix.(Nback{ii-1}).(SubsetNamesExtra{currSubset});
%                 end
%                 
%                 %scatter(currData(:,2),currData(:,1),18,'filled','o','MarkerFaceColor',[1 0 0]);
%                 hold on;
%                 [TheBins, TheSEM] = runAvg(currData(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
%                
%                 %Curve fitting
%                 %Fixing to empirical Values for curve fitting
%                 if currSubset == 2 || currSubset == 3 || currSubset == 5 || currSubset == 7 %Fix beta to empirical audio beta
%                     [fitresult, gof, fitCoeff] = fitDoGv2(currData(:,2),currData(:,1),empAbeta);
%                 else %Visual and other combined trials
%                     [fitresult, gof, fitCoeff] = fitDoGv2(currData(:,2),currData(:,1),empVbeta);
%                 end
%                 %Save estimated parameters
%                 fittedCurve.Paras.(Nback{ii-1}).(SubsetNamesExtra{currSubset}) = fitCoeff;
%                 shadedErrorBar(LB:plotTicks:UB,TheBins,TheSEM);
%                 DoG = plot(fitresult,'b',currData(:,2),currData(:,1)); DoG(1).LineWidth = 1.5; DoG(1).delete %This removes the blue dots. It's a duplicate of the scatter.
%               
%                 plot(Xlim,[0 0], 'k', 'LineWidth', 0.25); %Draws a line on Y axis.
%                 RunningAvg = plot(LB:plotTicks:UB, TheBins, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);%Running average plot
%                
%                 %Axis attribtues.
%                 AnalysisPlot = gca; %this will be used to set attributes of axis
%                 set(AnalysisPlot, 'XTick', XTICK); set(AnalysisPlot, 'Xlim', Xlim); set(AnalysisPlot, 'XTickLabel', XTICKLABEL);
%                 set(AnalysisPlot, 'YTick', YTICK); set(AnalysisPlot, 'Ylim', [YlimLB YlimUB]); set(AnalysisPlot, 'YTickLabel',YTICKLABEL);
%                 xlabel('Previous Trial -- Current Trial (Deg)'); ylabel('Response Error (Deg)');
%                 if median > 0; plottitle = sprintf('Group: N-%d %s Trials, %s\n %sÂ° bins, R-square = %.2f, RMSE = %.2f, Median', (ii-1), CondType, SubsetNamesExtra{currSubset}, UserInput{3}, gof.rsquare, gof.rmse);
%                 else plottitle = sprintf('Group: N-%d %s Trials, %s\n %sÂ° bins, R-square = %.2f, RMSE = %.2f', (ii-1), CondType, SubsetNamesExtra{currSubset}, UserInput{3}, gof.rsquare, gof.rmse); end; title(plottitle);
%                 legend([DoG(2) RunningAvg],'Fitted curve','Average','Location','northeast','Orientation','vertical'); legend('boxoff'); %Average is the name of the RunningAvg line
%                 hold off; clear currData RunningAvg TheBins DoGestimates DoGplot
%             end
%         end
    end
%      
    %% Confidence Interval
%     currName = fieldnames(GroupMatrix.N1.Perm);
%     %Fit data to all the permutations to get estimated alpha
%     for ii = 2:Nbacktrials+1
%         %For each comparisons
%         for currSubset = 1:size(SubsetNames,2)
%             %For each permutation
%             for currPerm = 1:permnum
%                 currField = GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).(SubsetNames{currSubset});
%                 %Fixing to empirical Values for curve fitting
%                 if currSubset == 2 || currSubset == 4 || currSubset == 6 %Fix beta to empirical audio beta
%                     [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empAbeta);
%                 else %Visual and other combined trials
%                     [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empVbeta);
%                 end
%                 fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNames{currSubset}).Matrix(currPerm,:) = fitCoeff;
%             end
%         end
%         
% %         %For each new comparisons
% %         if ii == 2
% %             for currSubset = 1:size(SubsetNamesExtra,2)
% % %                 %For each permutation
% %                 for currPerm = 1:permnum
% %                     if currSubset == 1 %VXVCombine
% %                         currField = [GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VV(:,[1 2]); ...
% %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VAV(:,[1 2]); ...
% %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).VAAV(:,[1 2])];
% %                     elseif currSubset == 4 %AXACombine
% %                         currField = [GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AA(:,[1 2]); ...
% %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AVA(:,[1 2]); ...
% %                             GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).AVVA(:,[1 2])];                        
% %                     else
% %                         currField = GroupMatrix.(Nback{ii-1}).Perm.(currName{currPerm}).(SubsetNamesExtra{currSubset});
% %                     end
% %                     Fixing to empirical Values for curve fitting
% %                     if currSubset == 2 || currSubset == 4 || currSubset == 6 %Fix beta to empirical audio beta
% %                         [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empAbeta);
% %                     else %Visual and other combined trials
% %                         [fitresult, gof, fitCoeff] = fitDoGv2(currField(:,2),currField(:,1),empVbeta);
% %                     end
% %                     fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNamesExtra{currSubset}).Matrix(currPerm,:) = fitCoeff;
% %                 end
% %             end
% %         end
%     end
    
    %% Confidence Interval Plot
    %For N2 and N3, AA and VV are AAA,VVV, AAAA,VVVV. so we can derive the
    %mean and SD, and then do significance testing. 
    %Now build the single SD based on the permutations, N-1 trials testing
    ii = 2;
    for currSubset = 1:size(SubsetNames,2)
        %Get estimated alphas of each subset analysis, all permutations
        EstAlpha.(SubsetNames{currSubset}).Matrix = fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNames{currSubset}).Matrix(:,1);
        temp_Distribution = sort(EstAlpha.(SubsetNames{currSubset}).Matrix,'ascend'); %To calculate C.I.
        %95% C.I.
        TheCI = abs(circ_rad2ang(temp_Distribution([minCI*permnum maxCI*permnum])));
        
        %Calculate the SD of the estimate alpha, convert EstAlpha to vector
        EstAlpha.(SubsetNames{currSubset}).CI = TheCI;
        %Calculate the mean of estimate alpha
        EstAlpha.(SubsetNames{currSubset}).mean = circ_rad2ang(circ_nanmean(EstAlpha.(SubsetNames{currSubset}).Matrix));
    end
    
%     %For the new analysis
%     for currSubset = 1:size(SubsetNamesExtra,2)
%         %Get estimated alphas of each subset analysis, all permutations
%         EstAlpha.(SubsetNamesExtra{currSubset}).Matrix = fittedCurve.Perm.Paras.(Nback{ii-1}).(SubsetNamesExtra{currSubset}).Matrix(:,1);
%         temp_Distribution = sort(EstAlpha.(SubsetNamesExtra{currSubset}).Matrix,'ascend'); %To calculate C.I.
%         %95% C.I.
%         TheCI = abs(circ_rad2ang(temp_Distribution([minCI*permnum maxCI*permnum])));   
%         %Calculate the SD of the estimate alpha, convert EstAlpha to vector
%         EstAlpha.(SubsetNamesExtra{currSubset}).CI = TheCI;
%         %Calculate the mean of estimate alpha
%         EstAlpha.(SubsetNamesExtra{currSubset}).mean = circ_rad2ang(circ_nanmean(EstAlpha.(SubsetNamesExtra{currSubset}).Matrix(:)));
%     end
%     
%     %For VVV, AAA, VVVV, AAAA
%     addStuff = {'VV','AA'}; counterA = 0; moreNames = {'VVV','AAA';'VVVV','AAAA'};
%     for addAnalysis = 1:2 %4 analysis here AAA, AAAA, VVV, VVVV
%         for curraddStuff = 1:2
%             counterA = counterA + 1;
%             currMatrix = fittedCurve.Perm.Paras.(Nback{addAnalysis+1}).(addStuff{curraddStuff}).Matrix(:,1);
%             EstAlpha.(moreNames{addAnalysis,curraddStuff}).mean = circ_rad2ang(circ_nanmean(currMatrix));
%             EstAlpha.(moreNames{addAnalysis,curraddStuff}).Matrix = currMatrix;
%             temp_Distribution = sort(currMatrix,'ascend'); %To calculate C.I.
%             %95% C.I.
%             TheCI = abs(circ_rad2ang(temp_Distribution([minCI*permnum maxCI*permnum])));
%             %Calculate the SD of the estimate alpha, convert EstAlpha to vector
%             EstAlpha.(moreNames{addAnalysis,curraddStuff}).CI = TheCI;
%         end
%     end
% 
%     %We now have all the CIs calculated. Time to plot them.
%     theConds = {{'VV','AA';fittedCurve.Paras.N1.VV(1),fittedCurve.Paras.N1.AA(1)}; ...
%         {'VVV','VVVV','AAA','AAAA'; fittedCurve.Paras.N2.VV(1), fittedCurve.Paras.N3.VV(1), fittedCurve.Paras.N2.AA(1), fittedCurve.Paras.N3.AA(1)}; ... 
%         {'VXA','VAV','VVA','VAA','VAAV','AXV','AVA','AVV','AVVA'; fittedCurve.Paras.N1.VXA(1),fittedCurve.Paras.N1.VAV(1),fittedCurve.Paras.N1.VVA(1),fittedCurve.Paras.N1.VAA(1),fittedCurve.Paras.N1.VAAV(1), ... 
%         fittedCurve.Paras.N1.AXV(1),fittedCurve.Paras.N1.AVA(1),fittedCurve.Paras.N1.AVV(1),fittedCurve.Paras.N1.AVVA(1)}};
%     %Empirical alpha for the above conditions as follows
%     xposition = [.2 .4 .6 .8 1 1.2 1.4 1.6 1.8]; %For organizing the plots along X-axis
%     for toPlot = 1:size(theConds,1);
%         figure; currInfo = theConds{toPlot,1};
%         for ii = 1:size(currInfo,2)
%             %Get the Matrix and the C.I. values previously calculated
%             temp_Distribution = EstAlpha.(currInfo{1,ii}).Matrix;
%             TheCI = EstAlpha.(currInfo{1,ii}).CI;
%             %Plot one SD of the permutated alpha
%             %permDist = plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',1.5); 
%             hold on;
%             %The following line is to make the bootstrap dist thicker, but not interfere with legend
%             %plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',20);
%             oneSD = plot([1 1]+xposition(ii),[-TheCI(1), TheCI(2)],'k');
%             %Plot empirical alpha
%             empAlpha = plot(1+xposition(ii),circ_rad2ang(currInfo{2,ii}),'rx','MarkerSize',6);
%         end
%         
%         %Axis Properties
%         AnalysisPlot = gca;
%         set(AnalysisPlot, 'Ylim', [-20 20]); set(AnalysisPlot, 'Xlim', [.9 3]); set(AnalysisPlot,'XTick',1+xposition(1:size(currInfo,2)));
%         set(AnalysisPlot, 'XTickLabel',currInfo(1,:));
%         ylabel('Previous Trial -- Current Trial (Deg)'); title('Empirical Alpha and 95% Confidence Interval');
%         legend([oneSD empAlpha],'95% C.I.','Empirical Alpha','location','southeast');
%         hold off
%     end

%     %when you run combined analysis, run this
%     %We now have all the CIs calculated. Time to plot them.
%     theConds = {'VV','AA','VXA','AXV','VV','VXV','AA','AXA','VXA','AXV'; ...
%         fittedCurve.Paras.N1.VV(1),fittedCurve.Paras.N1.AA(1), ...
%         fittedCurve.Paras.N1.VXA(1), fittedCurve.Paras.N1.AXV(1), ...
%         fittedCurve.Paras.N2.VV(1), fittedCurve.Paras.N2.VXV(1), ... 
%         fittedCurve.Paras.N2.AA(1), fittedCurve.Paras.N2.AXA(1) ...
%         fittedCurve.Paras.N2.VXA(1), fittedCurve.Paras.N2.AXV(1); ...
%         'N1','N1','N1','N1','N2','N2','N2','N2','N2','N2'};
%     %Empirical alpha for the above conditions as follows
%     xposition = [.2 .4 .6 .8 1 1.2 1.4 1.6 1.8 2.0]; %For organizing the plots along X-axis
%     
%     figure;
%     for ii = 1:size(theConds,2)
%         %Get the Matrix and the C.I. values previously calculated
%         temp_Distribution = sort(fittedCurve.Perm.Paras.(theConds{3,ii}).(theConds{1,ii}).Matrix(:,1),'ascend');
%         TheCI = abs(circ_rad2ang(temp_Distribution([minCI*permnum maxCI*permnum])));
%         %Plot one SD of the permutated alpha
%         %permDist = plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',1.5);
%         hold on;
%         %The following line is to make the bootstrap dist thicker, but not interfere with legend
%         %plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',20);
%         oneSD = plot([1 1]+xposition(ii),[-TheCI(1), TheCI(2)],'k');
%         %Plot empirical alpha
%         empAlpha = plot(1+xposition(ii),circ_rad2ang(theConds{2,ii}),'rx','MarkerSize',6);
%     end
%     
%     %Axis Properties
%     AnalysisPlot = gca;
%     set(AnalysisPlot, 'Ylim', [-4 4]); set(AnalysisPlot, 'Xlim', [.9 4]); set(AnalysisPlot,'XTick',1+xposition(1:size(theConds,2)));
%     set(AnalysisPlot, 'XTickLabel',theConds(1,:));
%     ylabel('Previous Trial -- Current Trial (Deg)'); title('Empirical Alpha and 95% Confidence Interval');
%     legend([oneSD empAlpha],'95% C.I.','Empirical Alpha','location','southeast');
%     hold off
    
    
    %Calculate the p values this way:
    %pval = ((sum((abs(fittedCurve.Perm.Paras.N1.VAV.Matrix(:,1)) >= abs(fittedCurve.Paras.N1.VAV.Mat(1))))+1)/10001);
    
    %Now that we have plotted the results, it's time to head back out the data
    %folder.
    cd ..;
    
    %Turn warnings on again 
    warning on MATLAB:legend:IgnoringExtraEntries
    warning on curvefit:prepareFittingData:removingNaNAndInf
    
catch error %This code will cd into a data folder, even if code fails halfway. This just cd out the data folder, so you can run again with minimal mouse clicks.
    cd ..;  
    %Turn warnings on again 
    warning on MATLAB:legend:IgnoringExtraEntries
    warning on curvefit:prepareFittingData:removingNaNAndInf 
    rethrow(error); 
end
