function [ CS_spike_train, US_spike_train ] = generate_training_stim(numtrials, ...
                                                                     CS_onset, CS_offset, CS_rate,...
                                                                     US_onset, US_offset, US_rate,...
                                                                     probe_start,probe_freq,slack)
% CS_onset              Can be a vector or a matrix.
%                       Each row represents timings for onsets
%                       When it is a matrix, different rows represent
%                        individual interleaved trials
%                       * leave a row as zeros if you want no stimulation
% CS_offset
% US_onset
% US_offset

if size(CS_onset,1) ~= size(CS_offset,1)
    error('CS : onset timings should be paired with offset timings (#of columns must match)');
end
if size(CS_onset,2) ~= size(CS_offset,2)
    error('CS : onset timings should be paired with offset timings (#of rows must match)');
end
if size(US_onset,1) ~= size(US_offset,1)
    error('US : onset timings should be paired with offset timings (#of columns must match)');
end
if size(US_onset,2) ~= size(US_offset,2)
    error('US : onset timings should be paired with offset timings (#of rows must match)');
end 
if ~isempty(find(CS_offset-CS_onset<0))
    error('CS_onset,offset : onset should not happen later than offset')
end
if ~isempty(find(US_offset-US_onset<0))
    error('US_onset,offset : onset should not happen later than offset')
end
alloffsets = [CS_offset US_offset];
% convert units

%slack = additional time after all offsets
%
    trial_length                        = ceil((max(alloffsets(:))+slack)*1000);          % ms
    CS_onset                            = ceil(CS_onset*1000);              % s
    CS_offset                           = ceil(CS_offset*1000);             % s
    US_onset                            = ceil(US_onset*1000);              % s
    US_offset                           = ceil(US_offset*1000);             % s
    
    
    
    CS_spike_train = zeros(numtrials,trial_length);
    US_spike_train = zeros(numtrials,trial_length);

% generate a 'block' = k trials of CSUS where k = # different structures
% within each trial, there can be single or multiple CS US bursts

isetup = 1;
for itrial = 1:numtrials
    
    isetup = mod(isetup + 1,size(CS_onset,1));
    if isetup == 0
        isetup = size(CS_onset,1);
    end
    CS_spike_train_per_trial = zeros(1,trial_length);
    
    for iburst = 1:size(CS_onset,2)
        CS_spike_train_per_trial(1,CS_onset(isetup,iburst):CS_offset(isetup,iburst)-1) = binornd(1,CS_rate/1000,[1,CS_offset(isetup,iburst)-CS_onset(isetup,iburst)]);
    end
    
    US_spike_train_per_trial = zeros(1,trial_length);
    
    for iburst = 1:size(US_onset,2)
       US_spike_train_per_trial(1,US_onset(isetup,iburst):US_offset(isetup,iburst)-1) = binornd(1,US_rate/1000,[1,US_offset(isetup,iburst)-US_onset(isetup,iburst)]);
    end
    
    CS_spike_train(itrial,:) = CS_spike_train_per_trial;
    if ~isempty(probe_start) && ~isempty(probe_freq)
        if itrial > probe_start && mod(itrial,probe_freq) ==0
            US_spike_train(itrial,:) = zeros(1,length(US_spike_train_per_trial));
        else
            US_spike_train(itrial,:) = US_spike_train_per_trial;
        end
    end
        
end
end
