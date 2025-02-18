%% Script for checking the timing of events in trials

timing_matrix = [];

for iTrial = 1:length(log.trial)
    timing_matrix = [timing_matrix; log.trial{iTrial}.onsets];  
end

% timing within trial
event_timings = timing_matrix(:,2:end) - timing_matrix(:,1:end-1);

response_period_timing = timing_matrix(:,8) - timing_matrix(:,5) + 1/60;

% timing between trials
tmp = timing_matrix(2:end,1) - timing_matrix(1:end-1,8);