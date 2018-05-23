function [V, ...
          exc,inh, reserve, archive_tracker,...
          record_activation_energy, ...
          readout_activation_energy, readout_refractory_timer, ...
          pause_trial] = ...
          ...               
          Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                         US_onset, US_offset, US_rate,...
                         trial_length, trail, numtrials,...
                         archive0, halt_at_pause,probe_start,probe_freq)

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%Model of Purkinje Pause%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                %Inspired by F. Johannson and G. Hesslow
                                %Model by J. Kim and M. Ricci
                                %Code by J. Kim and M. Ricci
                                %Last edited by J.Kim, 4/30/16

                                
disp('Starting experiment...');
disp('--Initializing variables');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE STIMULI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Module

    % activation
    record_activation_energy0           = 2;                               % (spikes) temporal summation of spikes needed for write switch
    reserve0                            = 1;                               % Max capacity of reserve
    record_tau_activation               = 70;                              % Time (ms) for record activation to increase by 40%
    record_trigger_on                   = 1;                               % Indicates new CS respects (=1) refractory period
    
    % release
    record_tau_release                  = 100;                              % Time (ms) for reserve to deplete by 40%
    
    % recharge                                     
    record_production_rate              = reserve0/125000;                 % Fraction of reserve plenished every time-step (2.5 mins until full charge)
    record_refractory_timer0            = 1000;                            % Refractory period of write module (ms)
   
    % unit properties
    record_evolution_noise              = 40;                              % Probability of recorder unit evolution delay on each time-step
    record_max_isi                      = size(archive0,2);                % Maximum learnable ISI
   
    
% Read module 

    % activation
    readout_activation_energy0      = 2;                                   % (spikes) temporal summation of spikes needed
    readout_tau_activation          = 200;                                 % Time (ms) for read activation energy to increase by 40%
    
    % release
    readout_consume_rate            = 1/300; %last worked = 1/160K       % Fraction of archive consumed durnig each read
    
    % recharge
    readout_refractory_timer0       = 2000;                                % Read module refractory period (ms)
   
    % unit property
    readout_unit_efficacy           = 45000;  %last worked = 25000K       % Current (mA) delievered per 1 recorder unit
    

% Purkinje cell parameters
    V_rest                          = -70;                                 % Resting potential (mV)
    V_threshold                     = -54;                                 % Spiking threshold(mV)
    V_spike                         = 10;                                  % Spiking potential(mV)
    V_hyperpol                      = -85;                                 % Hyperpolarization potential(mV)
    CS_weight                       = 10;                                  % Current (mA) per CS spike
    tau_membrane                    = 5;                                   % PC membrane time constant
    pacemaker_mag                   = 10;                                  % Current (mA) from each pacemaker event
    pacemaker_p                     = 0.3;                                 % Rate of pacemaker Poisson process
    
  
    
% Experimental criterion (pause detection)
    pause_criterion             = .25;                              % Minimum value of archive for pause detection
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE STIMULI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--Generating stim'); 

[CS_spike_train, US_spike_train]  = generate_training_stim(numtrials, ...
                                                         CS_onset, CS_offset, CS_rate,...
                                                         US_onset, US_offset, US_rate,probe_start,probe_freq,trail);


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE STATE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Voltage, Current and Switch Energies
    V                                   = zeros(numtrials,size(CS_spike_train,2));                 % Membrane Potential (mV)
    exc                                 = zeros(numtrials,size(CS_spike_train,2));
    inh                                 = zeros(numtrials,size(CS_spike_train,2));
    reserve                             = zeros(numtrials,size(CS_spike_train,2));                 % record of reserve every timestep
    archive_tracker                     = zeros(numtrials,size(archive0,2));                       % record of archive every trial
    record_activation_energy            = zeros(numtrials,size(CS_spike_train,2));                 % record of AE every timestep
    record_refractory_timer             = zeros(numtrials,size(CS_spike_train,2)); 
    readout_activation_energy           = zeros(numtrials,size(CS_spike_train,2));                 % record of AE every timestep
    readout_refractory_timer            = zeros(numtrials,size(CS_spike_train,2));                 % record of refractory timer every timestep
    pause_trial                         = [];
    
    
% cell
    V(1,1)                            = V_rest;                              % Initial membrane potential (mV)
    exc(1,1)                          = 0;
    inh(1,1)                          = 0;
    
% READ
    batch_reading                     = zeros(1,round(record_max_isi));      % beginning batch (=none)
    archive                           = archive0;                            % beginning reserve (=none)
    trigger_on                        = 0;                                   % beginning state (=off)
    readout_activation_energy(1,1)    = readout_activation_energy0;          % beginning activation energy
    readout_refractory_timer(1,1)     = 0;                                   % beginning refectory timer (=nonrefractory)
    
% WRITE
    batch_writing                     = zeros(1,round(record_max_isi));      % beginning batch (=none)
    reserve(1,1)                      = 0;                                   % beginning reserve (=none)
    record_activation_energy(1,1)     = record_activation_energy0;           % beginning activation energy
    record_refractory_timer(1,1)      = 0;      
    
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE CELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--Running trials'); 
for itrial = 1:numtrials
    disp(strcat('Trial',num2str(itrial)));
    for t = 2:size(CS_spike_train,2)





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE ONE TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % WRITE MODULE
        [batch_final, batch_writing, reserve(itrial,t), record_activation_energy(itrial,t),...
            record_refractory_timer(itrial,t),record_trigger_on] = ...
            ...
            WRITE(CS_spike_train(itrial,t),  US_spike_train(itrial,t),...
                   ...
                   batch_writing, ...
                   reserve(itrial,t-1), ...
                   record_activation_energy(itrial,t-1), ...
                   record_refractory_timer(itrial,t-1),...
                   ...
                   record_activation_energy0, ...
                   record_refractory_timer0,...
                   record_production_rate, ...
                   record_tau_release,...
                   record_tau_activation,...
                   record_evolution_noise,...
                   record_max_isi,...
                   reserve0,...
                   record_trigger_on);

        % Update archive with WRITE output
        archive = archive + batch_final;
        % READ MODULE
        [hyperpolarize, batch_reading, archive, trigger_on, readout_activation_energy(itrial,t), readout_refractory_timer(itrial,t)] = ...     
            ... 
            READ( CS_spike_train(itrial,t), ...
                     ...
                     batch_reading, ...
                     archive, ...
                     trigger_on, ...
                     readout_activation_energy(itrial,t-1), ...
                     readout_refractory_timer(itrial,t-1), ...    
                     ...
                     readout_activation_energy0,...
                     readout_refractory_timer0,...
                     readout_unit_efficacy,...
                     readout_consume_rate,...
                     readout_tau_activation);  


        % Non-pacemaker currents
        exc(itrial,t) = CS_weight*CS_spike_train(itrial,t); % CS excitation
        inh(itrial,t) = hyperpolarize; % Inhibition from archive

        %Calculate membrane potential
        if V(itrial,t-1) == V_spike;

            V(itrial,t) = V_hyperpol;

        elseif V(itrial,t-1) >= V_threshold

            V(itrial,t) = V_spike;

        else
            pacemaker =  pacemaker_mag*binornd(1,pacemaker_p);
            V(itrial,t) = V(itrial,t-1) - (V(itrial,t-1)-V_rest)/tau_membrane + exc(itrial,t) - inh(itrial,t)+ pacemaker;
        end





    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KEEP TRACK OF PER-TRIAL SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    archive_tracker(itrial,:) = archive; % to see the evolution of archive

    

    %detection of pause trial
    if halt_at_pause && itrial>3
        
        % compute the indices of isi(s)
        current_setup = mod(itrial,size(CS_onset,1));
        if current_setup == 0
            current_setup = size(CS_onset,1);
        end
        isi_idx = [];
        for iburst = 1:size(CS_onset,2)
            isi_idx = [isi_idx (US_onset(current_setup,iburst)-0.05)*1000:1:US_onset(current_setup,iburst)*1000];
            isi_idx = round(isi_idx);
        end
        % aggregate 4 trials, if the total number of spikes from 5ms before US until US is less than 1 , it is considered pause 
%             current_ISI_rate = length([find(V(itrial-3,isi_idx) ==  V_spike) ...
%                                        find(V(itrial-2,isi_idx) ==  V_spike) ...
%                                        find(V(itrial-1,isi_idx) ==  V_spike) ...
%                                        find(V(itrial,isi_idx) ==  V_spike)]);
%             fprintf(strcat('------',num2str(sum(archive)),' units in archive\n'));
%             fprintf(strcat('------',num2str(1000*current_ISI_rate),'Hz\n'));
%             if current_ISI_rate < 1
%                 pause_trial = itrial;
%                 break
%             end
        disp(['inh sum =', num2str(sum(inh(itrial,isi_idx)))]);
        if sum(inh(itrial,isi_idx)) > 200*size(CS_onset,2);
            pause_trial = itrial;
            break
        end
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPLY INTER-TRIAL CHANGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if itrial < numtrials
        vacation                = trial_length*1000 - size(CS_spike_train,2);
        V(itrial+1,1)           = V_rest;
        exc(itrial+1,1)         = 0;
        inh(itrial+1,1)         = 0;
        
        % READ
        batch_reading                            = zeros(1,round(record_max_isi));      % APPROXIMATED : because iti is much longer than the maximum encodable time
        trigger_on                               = 0;                                   
        readout_activation_energy(itrial+1,1)    = readout_activation_energy0;          % APPROXIMATED : because iti is much longer than the time it takes for replanishment
        readout_refractory_timer(itrial+1,1)     = max(readout_refractory_timer(itrial,end) - vacation,0);                                        
    
        % WRITE
        batch_writing                            = zeros(1,round(record_max_isi));      % APPROXIMATED : because iti is much longer than the maximum encodable time
        reserve(itrial+1,1)                      = min(reserve(itrial,end) + record_production_rate*vacation, reserve0);                                  % 
        record_activation_energy(itrial+1,1)     = record_activation_energy0;           % APPROXIMATED : because iti is much longer than the time it takes for replanishment
        record_refractory_timer(itrial+1,1)      = max(record_refractory_timer(itrial,end) - vacation,0);      
    end
end

if halt_at_pause
    disp(['Pause at trial ', num2str(pause_trial)]);
end
disp('done')
end