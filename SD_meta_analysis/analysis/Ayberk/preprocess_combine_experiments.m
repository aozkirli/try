clc  % Clear command window

% Get the full path of the current script
scriptPath = mfilename('fullpath');

% Extract directory path of the current script
[currentDir, ~, ~] = fileparts(scriptPath);

currentDir = [currentDir filesep '..' filesep '..' filesep 'data']
% Navigate to the current directory
cd(currentDir)

% Remove existing 'studies' folder if it exists, suppress warning ('s')
rmdir('studies', 's')

% Create a new 'studies' folder
mkdir('studies');

% Navigate back to the current directory
cd(currentDir)

% Navigate to the 'experiments' folder
cd('experiments')

% List all .mat files in the 'experiments' folder
matfiles = dir('*.mat'); matfiles = string({matfiles.name});


% Extract the experiment names from the file names
expnames = regexp(matfiles, '^(.*?_20\d{2}\w?)(?=[_.])', 'match', 'once');

% Tabulate occurrences of each unique experiment name
summary = tabulate(expnames);

% Display the summary table
disp(summary)

% Loop over each unique experiment name
for k = 1:size(summary, 1)
    
    % Navigate to the 'experiments' folder
    cd(currentDir)
    cd('experiments')
    
    % Initialize an empty table to store data for this experiment
    tblk = table();
    
    % Get the current experiment name from the summary
    target_name = summary(k, 1);
    
    % Find indices of all files that match this experiment name
    index = find(strcmp(cellstr(expnames), target_name{:}));
    
    % Loop through all matching files
    for i = 1:numel(index)
        
        % Load the data from the .mat file
        load(matfiles(index(i)));
        
        % Add experiment number to the table
        tbl.expnum = i .* ones(size(tbl.obs));
        
        % Concatenate the data into the cumulative table
        tblk = vertcat(tblk, tbl);
    end
    
    % Get the name of the current file to rename later
    rename = matfiles(index(i));
    
    % Remove '_ExpX' tags from the file name
    for jj = 1:summary{k,2}  % Possible experiment indexes
        rename = strrep(rename, ['_Exp' num2str(jj)], '');
        rename = strrep(rename, ' ', '');  % Remove spaces
    end
    
    % Add an index number to the file name
    if k < 10
        num = ['0' num2str(k)];  % Pad single digit with 0
    else
        num = num2str(k);
    end
    
    % Update the file name with the new index
    rename = [num '_' char(rename)];
    
    % Navigate to the 'studies' folder to save the results
    cd(['..' filesep 'studies'])
    
    % Save the cumulative table as a .mat file
    tbl = tblk;
    save(rename, 'tbl')
    
    % Save the table as a .csv file
    writetable(tbl, [rename(1:end-4) '.csv'], 'delimiter', ';')
end
