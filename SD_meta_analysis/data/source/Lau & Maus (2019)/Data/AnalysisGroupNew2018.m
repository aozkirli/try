%Serial dependence multi new group analysis for V-only and A-only
%Pools all the data together, ignores all indv sbj errors, and treats the
%data as 1 yuuuge sbj
try
    %Find data folder in current directory. If you start at the data folder, you'll get an error.
    pathname = cd; DataFolder = fullfile(pathname,'Data'); cd Data;
    %Data folder exist in this directory
    
    %% Variables
    Nbacktrials = 1; AVmatrix = []; NumOfTrials = (1:104)'; VisualMaxResponseTime = 15; AudioMaxResponseTime = 20;
    LB = -0.5*pi; UB = 0.5*pi; Degrees = 2*pi; median = 0; keyWord = {'multi-1','multi-2'}; bindeg = 15/360*Degrees; 
    OutlierRemove = 1; ToneSizeFilenamePosition = 9; wrapDeg = 0.5*pi; maxwrapDeg = 2*wrapDeg; plotTicks = circ_ang2rad(0.5); %This determines size of increment in radians
    SimsPerTrial = 2000; permnum = 5000; var1 = 1; var2 = 6; var3 = 1; var4 = 1; %columns
    SbjStart = 1; conditions = {'V','A'}; Nbacks = {'N1','N2','N3','N4','N5'};
    %Ticks and Tick labels
    XTICK = [-pi/2 -pi/4 0 pi/4 pi/2]; Xlim = [LB UB]; XTICKLABEL = {'-90°','-45°','0°','45°','90°'};
    YlimLB = -5/180 * pi; YlimUB = 5/180 * pi;  YTICK = [YlimLB:(2.5/180*pi):YlimUB]; YTICKLABEL = {'-5°','-2.5°','0°','2.5°','5°'};
   
    %Create a variable with Permnum size Nback names and Permnum names
    for ii = 1:permnum
        Nback{ii,1} = sprintf('N%d',ii);
        permname{ii,1} = sprintf('a%d',ii);
    end

    %% All analysis Lazy mode.
    for lazyMode = 2%1:2
        allnames = dir; x = 1; %counter
        for i = 1:size(allnames,1)
            %If the string matches the filename keyword, save the name in variable
            keywordcheck = strfind(allnames(i).name,keyWord{lazyMode});
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
        for i = 1:size(currname,2); currsession{i} = sessionsize; end
        
        %Join all the names and sessions together
        UserInput{1} = strjoin(currname); UserInput{2} = strjoin(currsession);
        clear x allnames namelist currname currstring currsession sessions
        
        [Raw.(conditions{lazyMode}), SbjNames.(conditions{lazyMode}), SessionInput, ToneMatrixSize] = Loadfiles(UserInput{1},UserInput{2},keyWord{lazyMode},ToneSizeFilenamePosition);
        
        %% Remove long response time trials
        Raw_Matrix = Raw.(conditions{lazyMode});
        ResponseT = [Raw_Matrix(:,4) - Raw_Matrix(:,2), Raw_Matrix(:,7)]; %Show RespT and modality
        if lazyMode == 1 && strcmpi(keyWord{lazyMode},'multi-1') == 1
           ResponseT = [(ResponseT(:,1) > VisualMaxResponseTime),(ResponseT(:,2) > AudioMaxResponseTime)];
        elseif lazyMode == 2 && strcmpi(keyWord{lazyMode},'multi-2') == 1
           ResponseT = [(ResponseT(:,1) > AudioMaxResponseTime),(ResponseT(:,2) > AudioMaxResponseTime)];
        end
        Raw_Matrix(ResponseT(:,1),:) = NaN; Raw_Matrix(ResponseT(:,2),:) = NaN; %Convert the response for that trial to NaN
   
%         %% Mean Subtract in RAW Response
%         temp_Raw = Raw_Matrix(:,[1 6]); %Stim, Resp
%         for currStim = min(temp_Raw(:,1)):max(temp_Raw(:,1)) %all 178 stimulus
%             %Get the matrix for the current stimulus number
%             temp_Mat = temp_Raw(temp_Raw(:,1) == currStim,:);
%             %Mean values
%             temp_Mat(:,2) = circ_ang2rad(temp_Mat(:,2).*2); %Change to rad
%             temp_Mean = circ_mean(temp_Mat(:,2)); %Mean error in response
%             Unity_to_Mean_Distance =  temp_Mean - circ_ang2rad(currStim*2); %Error from Unity 21
%             temp_Mat(:,2) = temp_Mat(:,2) - Unity_to_Mean_Distance; %Mean subtract using distance from unity
%             temp_Mat(:,2) = round(circ_rad2ang(temp_Mat(:,2))./2); %Convert to deg, rounded off
%             %Wrapping
%             temp_Mat(temp_Mat(:,2) < 0,2) = temp_Mat(temp_Mat(:,2) < 0, 2) + 180;
%             temp_Mat(temp_Mat(:,2) > 179,2) = temp_Mat(temp_Mat(:,2) > 179, 2 ) - 180;
%             temp_Raw(temp_Raw(:,1) == currStim,:) = temp_Mat(:,[1 2]);
%         end
%         
%         %Replace the raw response values after cleanup
%         Raw_Matrix(:,6) = temp_Raw(:,2);
        
        Analy.(conditions{lazyMode}).Matrix = Raw_Matrix(:,var2) - Raw_Matrix(:,var1);%Resp - stim
        %Nback
        for aa = 1:Nbacktrials
            Analy.(conditions{lazyMode}).Matrix(:,aa+1) = circshift(Raw_Matrix(:,var3),aa) - Raw_Matrix(:,var4);
            NRemove(aa,:) = [aa, aa+size(NumOfTrials,1):size(NumOfTrials,1):size(Analy.(conditions{lazyMode}).Matrix,1)];
        end
        %Convert to angular info
        Analy.(conditions{lazyMode}).Matrix = circ_ang2rad(Analy.(conditions{lazyMode}).Matrix.*(360/ToneMatrixSize(1)));
        
        %Wrap data
        for ff = 1:Nbacktrials+1
            Temp_currData = Analy.(conditions{lazyMode}).Matrix;
            for loopTheWrap = 1:5 %Loop this a few times, more thorough
                Temp_currData(Temp_currData > wrapDeg) = Temp_currData(Temp_currData > wrapDeg) - maxwrapDeg;
                Temp_currData(Temp_currData < -wrapDeg) = Temp_currData(Temp_currData < -wrapDeg) + maxwrapDeg;
                %%Review comments -- gotta wrap these too
                %temp(temp > wrapDeg) = temp(temp > wrapDeg) - maxwrapDeg;
                %temp(temp < -wrapDeg) = temp(temp < -wrapDeg) + maxwrapDeg;
            end
            Analy.(conditions{lazyMode}).Matrix = Temp_currData;
            clear Temp_currData
        end
        
        %%Subtract the mean of the new calculations
        %temp(:,1) = temp(:,1) - circ_nanmean(temp(:,1));
        %Analy.(conditions{lazyMode}).Review = temp;
        
        %Calculate Group mean error and SD, then subtract mean error
        Gmeanerror = circ_nanmean(Analy.(conditions{lazyMode}).Matrix(:,1));
        Gstderror = circ_nanstd(Analy.(conditions{lazyMode}).Matrix(:,1));
        %Subtract mean error
        currData =  Analy.(conditions{lazyMode}).Matrix(:,1);
        currData = currData - Gmeanerror;
        %Remove outlier
        currData(currData > Gstderror * 3) = NaN;
        currData(currData < Gstderror * -3) = NaN;
        %Calculate mean error and SD post outlier remove
        GmeanerrorPost = circ_nanmean(currData);
        GstderrorPost = circ_nanmean(currData);
%         %Subtract errors again
%         currData = currData - GmeanerrorPost;
%         %Return to matrix
        Analy.(conditions{lazyMode}).Matrix(:,1) = currData;
        clear currData
        
       % %% Reviewer suggested analysis (1-back only)
       % temp = Analy.(conditions{lazyMode}).Matrix;
       % temp(temp(:,2) < 0,:) =  temp(temp(:,2) < 0,:) .* -1;
       % Analy.(conditions{lazyMode}).Matrix = temp;
        
        %% Plots
        %Curve fitting and stuff 
        for i = 1:Nbacktrials
            figure;
            currData =  Analy.(conditions{lazyMode}).Matrix(:,[1 i+1]);
            %Remove trials
            if i == 1
                currData(NRemove(i,:),1) = NaN;
            else
                currData(NRemove(1:i,:),1) = NaN;
            end
            %Remove NaN rows
            currData(isnan(currData(:,1)),:) = [];
            
            %scatter(currData(:,2),currData(:,1),18,'filled','o','MarkerFaceColor',[1 0 0]); 
            hold on;
            [TheBins,TheSEM] = runAvg(currData(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
            
            %Curve fitting
            [fitresult, gof, fitCoeff] = fitDoGv2(LB:circ_ang2rad(0.5):UB,TheBins');%fitDoGv2(currData(:,2),currData(:,1));
            shadedErrorBar(LB:plotTicks:UB,TheBins,TheSEM);
            DoG = plot(fitresult,'b',currData(:,2),currData(:,1)); DoG(1).LineWidth = 1.5; DoG(1).delete %This removes the blue dots. It's a duplicate of the scatter.
            
            %fitted curve parameters each sbj
            fittedCurve.Paras.(Nbacks{i}).(conditions{lazyMode}).Mat = fitCoeff;  fittedCurve.Paras.(Nbacks{i}).(conditions{lazyMode}).gof = gof;
            plot(Xlim,[0 0], 'k', 'LineWidth', 0.25); %Draws a line on Y axis.
            RunningAvg = plot(LB:plotTicks:UB, TheBins, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);%Running average plot
            
            %Axis attribtues.
            AnalysisPlot = gca; %this will be used to set attributes of axis
            set(AnalysisPlot, 'XTick', XTICK); set(AnalysisPlot, 'Xlim', Xlim); set(AnalysisPlot, 'XTickLabel', XTICKLABEL);
            set(AnalysisPlot, 'YTick', YTICK); set(AnalysisPlot, 'Ylim', [YlimLB YlimUB]); set(AnalysisPlot, 'YTickLabel',YTICKLABEL);
            xlabel('Previous Trial -- Current Trial (Deg)'); ylabel('Response Error (Deg)');
            legend([DoG(2) RunningAvg],'Fitted curve','Average','Location','northeast','Orientation','vertical'); legend('boxoff'); %Average is the name of the RunningAvg line
            hold off; clear RunningAvg TheBins DoGestimates DoGplot
            
            %if i == 5 %Build a matrix of permnum alphas and betas for N1 trials, then find the average (for sig testing purposes)
            %    currName = {'VMat','AMat'};
            %    for ttt = 1:permnum
            %    [~,~,fitCoeff] = fitDoGv2(currData(:,2),currData(:,1));
            %        fittedCurve.Paras.(Nbacks{i}).(currName{lazyMode}).Mat(ttt,:) = fitCoeff;
            %    end
            %    clear currData
            %end
        end
    end
    
%     %Plot for reviewer's comment
%     figure; currData = temp;
%     currData(NRemove(1,:),1) = NaN; %Remove 1st trial of each block
%     %Remove NaN rows
%     currData(isnan(currData(:,1)),:) = []; hold on;
%     [TheBins,TheSEM] = runAvg(currData(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
%     %TheBins = TheBins - TheBins(1); %Rescale the runAvg
%     %Curve fitting
%     [fitresult, gof, fitCoeff] = fitDoGv2(LB:circ_ang2rad(0.5):UB,TheBins);%fitDoGv2(currData(:,2),currData(:,1));
%     shadedErrorBar(LB:plotTicks:UB,TheBins,TheSEM);
%     DoG = plot(fitresult,'b',currData(:,2),currData(:,1)); DoG(1).LineWidth = 1.5; DoG(1).delete %This removes the blue dots. It's a duplicate of the scatter.
%     
%     %fitted curve parameters each sbj
%     fittedCurve.Paras.(Nbacks{1}).(conditions{lazyMode}).Mat = fitCoeff;  fittedCurve.Paras.(Nbacks{1}).(conditions{lazyMode}).gof = gof;
%     plot(Xlim,[0 0], 'k', 'LineWidth', 0.25); %Draws a line on Y axis.
%     RunningAvg = plot(LB:plotTicks:UB, TheBins(), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);%Running average plot
%     
%     %Axis attribtues.
%     AnalysisPlot = gca; %this will be used to set attributes of axis
%     set(AnalysisPlot, 'XTick', XTICK); set(AnalysisPlot, 'Xlim', [0 UB]); set(AnalysisPlot, 'XTickLabel', XTICKLABEL);
%     set(AnalysisPlot, 'YTick', YTICK); set(AnalysisPlot, 'Ylim', [YlimLB YlimUB]); set(AnalysisPlot, 'YTickLabel',YTICKLABEL);
%     xlabel('Previous Trial -- Current Trial (Deg)'); ylabel('Response Error (Deg)');
%     legend([DoG(2) RunningAvg],'Fitted curve','Average','Location','northeast','Orientation','vertical'); legend('boxoff'); %Average is the name of the RunningAvg line
%     hold off; clear RunningAvg TheBins DoGestimates DoGplot
    
    %% Permutation
    for lazyMode = 2%1:2
        Perm.Raw.(conditions{lazyMode}) = PermResponse(Raw.(conditions{lazyMode}),size(Raw.(conditions{lazyMode}),1),permnum);
        %Build permutation matrix
        currName = fieldnames(Perm.Raw.(conditions{lazyMode})); %Getting all the fieldnames within the structure
        for qq = 1:length(fieldnames(Perm.Raw.(conditions{lazyMode})))
            %Current matrix
            currData = Perm.Raw.(conditions{lazyMode}).(currName{qq});
            %Response error
            temp_currData = currData(:,var2) - currData(:,var1);
            %Nback trials stored in columns
            for rr = 1:Nbacktrials
                temp_currData(:,rr+1) = circshift(currData(:,var3),rr) - currData(:,var4);
            end            
            %Convert Degrees
            temp_currData(:,1:(Nbacktrials+1)) = circ_ang2rad(temp_currData(:,1:(Nbacktrials+1))).*(360/ToneMatrixSize(1)); 
            
            %% Reviewer perm transformation
            for z = 1:size(temp_currData,1)
                if temp_currData(z,1) <= 0 && temp_currData(z,2) <= 0
                    temp_currData(z,:) = temp_currData(z,:) * -1;
                elseif temp_currData(z,1) > 0 && temp_currData(z,2) < 0
                    temp_currData(z,:) = temp_currData(z,:) * -1;
                end
            end
            
            %Wrap matrix
            for ff = 1:Nbacktrials+1
                WrapMatrix = temp_currData(:,ff);
                for loopTheWrap = 1:5 %Loop this a few times, more thorough
                    WrapMatrix(WrapMatrix > wrapDeg) = WrapMatrix(WrapMatrix > wrapDeg) - maxwrapDeg;
                    WrapMatrix(WrapMatrix < -wrapDeg) = WrapMatrix(WrapMatrix < -wrapDeg) + maxwrapDeg;
                    temp_currData(temp_currData > wrapDeg) = temp_currData(temp_currData > wrapDeg) - maxwrapDeg;
                    temp_currData(temp_currData < -wrapDeg) = temp_currData(temp_currData < -wrapDeg) + maxwrapDeg;
                end
                temp_currData(:,ff) = WrapMatrix;
                clear WrapMatrix
            end
            
            %Review comment -- subtraction
            temp_currData(:,1) = temp_currData(:,1) - circ_nanmean(temp_currData(:,1));
            
            %Remove Outlier, mean subtraction
            currResp = temp_currData(:,1);
            currResp(currResp > Gstderror * 3) = NaN;
            currResp(currResp < Gstderror * -3) = NaN;
            %temp_currData = temp_currData - Gmeanerror - GmeanerrorPost;
            temp_currData = temp_currData - Gmeanerror;
            %Return to matrix
            Perm.Analy.(conditions{lazyMode}).(currName{qq}) = temp_currData;
            
            %For each Nback
            for ii = 1%1:3
                %Running average and fittedCoeff for each permutated set
               % for zz = 1:permnum
                    ttCurr_data = temp_currData; %Perm.Analy.(conditions{lazyMode}).(currName{zz});
                    if ii == 1
                        ttCurr_data(NRemove(ii,:),1) = NaN;
                    else
                        ttCurr_data(NRemove(1:ii,:),1) = NaN;
                    end
                    %Remove NaN rows
                    ttCurr_data(isnan(ttCurr_data(:,1)),:) = [];
                    ttCurr_data = ttCurr_data(:,[1 ii+1]); %select the current Nback column
                    %ttPerm_means(zz,:) = runAvg(ttCurr_data,2,bindeg,LB,UB,2,1,Degrees,median);
                    %PermRunAvg = runAvg(ttCurr_data(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
                    PermRunAvg = runAvg(ttCurr_data(:,[1 2]),2,bindeg,LB,UB,2,1,Degrees,median);
                    PermRunAvg = PermRunAvg - PermRunAvg(1); %Rescale the running avg
                    %Fitting and fixing the beta to the empirical value
                    if lazyMode == 1 %Visual
                        [fitresult, gof, fitCoeff] = fitDoGv2(LB:circ_ang2rad(.5):UB,PermRunAvg',fittedCurve.Paras.N1.V.Mat(2));%fitDoGv2(ttCurr_data(:,2),ttCurr_data(:,1),fittedCurve.Paras.N1.V(2));
                    elseif lazyMode == 2 %Audio
                        %[fitresult, gof, fitCoeff] = fitDoGv2(0:circ_ang2rad(.5):UB,PermRunAvg',fittedCurve.Paras.N1.A.Mat(2));
                        [fitresult, gof, fitCoeff] = fitDoGv2(0:circ_ang2rad(.5):UB,PermRunAvg(181:end)',fittedCurve.Paras.N1.A.Mat(2));%fitDoGv2(ttCurr_data(:,2),ttCurr_data(:,1),fittedCurve.Paras.N1.A(2));
                    end
                    fittedCurve.Perm.Paras.(Nbacks{ii}).(conditions{lazyMode}).Matrix(qq,[1 2]) = fitCoeff;
                %end
            end            
        end
        %Curve fit permutation
        fprintf('Curve fit permutation starting for %s ...\n',(conditions{lazyMode})); clear Perm;
    end
    
%     %Plot Confidence Interval
%     figure;
%     lazyMode = 1; currNback = 0; xposition = [.2 .4 .6 .8 1 1.2]; %For organizing the plots along X-axis
%     for ii = 1:size(xposition,2)
%       
%         %Switch to next modality
%         if ii > 3
%             lazyMode = 2;
%         end
%         
%         %Reset Nback
%         if currNback >= 3
%             currNback = 0;
%         end
%         %Current Nback trial
%         currNback = currNback + 1;
%                 
%         %Find the 2.5% and 97.5% value in the permutation
%         temp_Distribution = sort(fittedCurve.Perm.Paras.(Nbacks{currNback}).(conditions{lazyMode}).Matrix(:,1),'ascend');
%         TheCI = abs(circ_rad2ang(temp_Distribution([.025*permnum .975*permnum])));
%         
%         %Calculate the SD of the permutated Alpha in degrees
%         EstAlpha.(conditions{lazyMode}).(Nbacks{currNback}).SD = TheCI;
%         EstAlpha.(conditions{lazyMode}).(Nbacks{currNback}).mean = circ_rad2ang(circ_nanmean(fittedCurve.Perm.Paras.(Nbacks{currNback}).(conditions{lazyMode}).Matrix(:,1)));
%         
%         %Plot one SD of the permutated alpha
%         %permDist = plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',1.5);
%         hold on;
%         %The following line is to make the bootstrap dist thicker, but not interfere with legend
%         %plot([1 1]+xposition(ii),circ_rad2ang([min(temp_Distribution) max(temp_Distribution)]),'g','Linewidth',20);
%         oneSD = plot([1 1]+xposition(ii),[-TheCI(1), TheCI(2)],'k');
%         %Plot empirical alpha
%         currAlpha = fittedCurve.Paras.(Nbacks{currNback}).(conditions{lazyMode});
%         empAlpha = plot(1+xposition(ii),circ_rad2ang(currAlpha(1,1)),'rx','MarkerSize',6);
%     end
%     
%     %Axis Properties
%     AnalysisPlot = gca;
%     set(AnalysisPlot, 'Ylim', [-7 6]); set(AnalysisPlot, 'Xlim', [.9 2.5]); set(AnalysisPlot,'XTick',[1.2 1.4 1.6 1.8 2.0 2.2]);
%     set(AnalysisPlot, 'XTickLabel',{'V-only N1','V-only N2','V-only N3','A-only N1','A-only N2','A-only N3'});
%     ylabel('Previous Trial -- Current Trial (Deg)'); title('Empirical Alpha and 95% Confidence Interval - Experiment 1');
%     legend([oneSD empAlpha],'95% C.I.','Empirical Alpha','location','southeast');
%     hold off
%     
    %P-value
    %pval = ((sum((abs(ReviewPerm.AudPerms(:,1)) >= abs(ReviewPerm.AudParas(1))))+1)/10001);
        
    cd ..
catch ERR
    cd ..
    rethrow(ERR); 
end