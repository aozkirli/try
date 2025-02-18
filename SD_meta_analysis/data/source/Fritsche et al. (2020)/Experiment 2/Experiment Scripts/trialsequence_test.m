clear variables; close all;

rng('shuffle')

nBlocks     = 8;

numConds            = 2;
cbOrder             = 1;
reps                = 25;
omitSelfAdjacencies = 0;

nSubjects = 24;

%% Generate trialsequences
for iSubject = 1:nSubjects
    
    for iBlock = 1:nBlocks
        
       v_pos_sequence{iSubject}{iBlock} = carryoverCounterbalance(numConds,cbOrder,reps,omitSelfAdjacencies);
       
    end
    
end

%% Evaluate trialsequences
for iBack = 1:10
    
    for iSubject = 1:nSubjects
        
        for iBlock = 1:nBlocks
            
            same_loc{iBack}(iSubject,iBlock) =  sum(v_pos_sequence{iSubject}{iBlock}(iBack+1:end) == v_pos_sequence{iSubject}{iBlock}(1:end-iBack));
            diff_loc{iBack}(iSubject,iBlock) =  sum(v_pos_sequence{iSubject}{iBlock}(iBack+1:end) ~= v_pos_sequence{iSubject}{iBlock}(1:end-iBack));
            
        end
        
    end
    
    mean_same_loc{iBack} = mean(same_loc{iBack},2);
    mean_diff_loc{iBack} = mean(diff_loc{iBack},2);
    
end