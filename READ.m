                                    % FINAL OUTPUT
function [ hyperpolarize, ...          % amount of instantaneous hyperpolarization
           ...
           ...                      % RECURRENT STATE VALUES
           batch_reading,...           %
           archive, ...                %
           trigger_on, ...             %
           activation_energy, ...      %
           refractory_timer] = ...     %
  ...           
  ...
  ...                               % EXTERNAL STATE VALUES
  READ( CS_spike, ...               % current presynaptic spike (1 or 0)
           ...
           ...                      % RECURRENT STATE VALUES
           batch_reading, ...          % distribution vector of recorded timer units undergoing readout (it's length = remaining readout)
           archive, ...                % distribution vector of total recorded timer units
           trigger_on, ...             % is readout taking place? (are timer units primed?)
           activation_energy, ...      % remaining activation energy
           refractory_timer, ...       % remaining refractory period
           ...
           ...                      % INITIAL STATE VALUES
           activation_energy0, ...      % activation energy (= close temporal sum of presynaptic spikes) required to start readout (needs to be small)
           refractory_timer0, ...      % refractory period (has to be >> readout_timer0)
           ...
           ...                      % STRUCTURAL PARAMETERS
           unit_efficacy, ...          % inhibitory power per timer unit
           consume_rate, ...           % portion of recorded timer unit reserve consumed per readout
           tau_activation)             % how dense should the presynaptic spikes be to trigger the switch?

                                       % ** all time-dependent values should be normalized for step size 
                                       
% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SWITCH : a 'pulse' thresholded activation function with refractoriness
%          triggered by PF spikes, not postsynaptic voltage
%               => predicts no CR upon direct pkj cell stimulation during test
%               => predicts no repetitive CR upon CS-CS-US test or CS-US-CS 
%                  if the second CS immediately follows
%
% RELEASE (of recorded units) : impulse-like exponential
%               => only bio-algorithmically plausible alternative (= sustained 
%                  constant) doesn't make sense
%
% EFFECT (of recorded units) : impulse or sustained
%               => predicts different inhibition profiles
%

% Edited 2/28/2016 by M. G. Ricci: now require archive to encode > 100 ms
% in order to initiate reading. See line 85. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise = 0; % (standard dev) inhibitory effect has random sampling-like property (binomial distribution)
hyperpolarize = 0;

if trigger_on == 1 %If recorder units are primed
    
    %Countdown for readout refractory period
    refractory_timer = refractory_timer - 1;
       
    if length(batch_reading) == 1 %If all units have been read
        
        hyperpolarize = batch_reading(1);
        trigger_on = 0; %deprime recorder units
        batch_reading = 0;
        
    elseif length(batch_reading) > 1
        
        % apply hyperpolarization
        hyperpolarize = batch_reading(1)*unit_efficacy*min(max(normrnd(1,noise),0),2);
        batch_reading = batch_reading(2:end);
        
    else
        error('length of batch cannot be less than 1');
    end
    
elseif trigger_on == 0
    if refractory_timer > 0
        % refractory period
        refractory_timer = refractory_timer - 1;

    elseif refractory_timer <= 0
        if activation_energy <= 0 && length(archive) > 100 %Added to account for 100 ms threshold
            % start trigger
            refractory_timer = refractory_timer0;
            activation_energy = activation_energy0;
            trigger_on = 1;
            % apply consumption
            batch_reading = archive*consume_rate;
            archive = archive*(1-consume_rate);
            
        elseif activation_energy > 0
            
            % adjust activation energy based on recuperation and presynaptic spikes 
            activation_energy = activation_energy + (activation_energy0-activation_energy)/tau_activation;
            activation_energy = activation_energy - CS_spike;
            
        else
            error('looks like activation_energy is not a number');
        end
    else
        error('looks like refractory_timer is not a number');
    end
else
    error('trigger_on is either 1 or 0');
end

end
    
