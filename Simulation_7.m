%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 1: The CSUS/USUS Ratio %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% 12/3/2015 Changes:
% Unit efficacy added lines 20, 47
% change_cs added lines 21, 47


clear
% user-defined parameters ======> edit only this part
    CSUS            = .1:.01:1; %s
    CS_duration     = 0.5; %s
    CS_rate         = 100; %fire/s
    USUS            = 20; %s
    US_duration     = .05; %s
    US_rate         = 500; %fire/s
    numtrials       = 6000;

% create parameters for run_experiment (everything is in seconds)
    if USUS<= max(CSUS) + US_duration
        error(strcat('USUS is too short -- should be greater than ',num2str(max(CSUS) + US_duration)));
    end                      
    trial_length                    = USUS; 
    CS_onset                        = 0.001; %0 seconds returns no learning
    CS_offset                       = CS_duration;
    US_onset                        = CSUS;
    US_offset                       = CSUS + US_duration;
    archive0                        = zeros(1,2500); % default information (=zeros, max_isi-long)
    firing_likelihood               = .5;
    alpha                           = .5; %tolerance to spiking
    halt_at_pause                   = 1;
    probe_start                     = Inf;
    probe_freq                      = 0;
                        
%% RUNy

%%% RUN A : CSUS varies while USUS fixed
for i = 1:length(CSUS)
    disp(strcat('CSUS = ',num2str(CSUS(i))));
    % run experiment, only return inhibition history
    [V,~, inh, ~, ~, ~, ~, ~,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset(i), US_offset(i), US_rate,...
                                           trial_length, 0.5, numtrials,...
                                           archive0,halt_at_pause,...
                                           probe_start, probe_freq);
 
    % determine trials-to-acquisition
    if isempty(pause_trial)
        pause_startA(i) = NaN;
    else
        pause_startA(i) = pause_trial;
    end
end

%
%%% RUN B : CSUS and USUS vary together
USUS = 80*(CSUS+US_duration); 
trial_length = USUS;%?? 100
numtrials = max(pause_startA);
for i = 1:length(CSUS)
    % run experiment, only return inhibition history
     [V,~, inh, ~, ~, ~, ~, ~,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset(i), US_offset(i), US_rate,...
                                           trial_length(i), 0.5, numtrials,...
                                           archive0,halt_at_pause,...
                                           probe_start, probe_freq);
 
    % determine trials-to-acquisition
    if isempty(pause_trial)
        pause_startB(i) = NaN;
    else
        pause_startB(i) = pause_trial;
    end
end

 % PLOT RESULTS
 plot(CSUS, pause_startA,'r',CSUS,pause_startB,'b');
 toc