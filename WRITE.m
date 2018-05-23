function [batch_final,...              % distribution of the current batch of recorded timer unit (upon US)
          ...
          ...                       % RECURRENT STATE VALUES
          batch_writing,...            % distribution of the current batch of recorded timer unit (before US)
          reserve,...                  % 
          activation_energy,...% Activation energy is the quantity that detects CS onset. When the enough
          ...                              % CS spikes happen quickly, activation energy decreases below a threshold
          ...                              % (here 0). 
          record_refractory_timer,...
          trigger_on] = ...     % 
...  
...
...                                 % EXTERNAL STATE VALUES
   WRITE(CS_spike, ...                % current CS presynaptic spike (1 or 0)
          US_spike, ...                % current US presynaptic spike (1 or 0)
         ...
          ...                       % RECURRENT STATE VALUES
          batch_writing,...            % distribution of the current batch of recorded timer unit (before US)
          reserve,...                  % remaining quantity of recordable timer units
          activation_energy, ...       % remaining activation energy
          record_refractory_timer, ...
          ...
          ...                       % INITIAL STATE VALUES
          activation_energy0, ...      % total activation energy.
          record_refractory_timer0,...
          ...
          ...                       % STRUCTURAL PARAMETERS
          production_rate,...          % units per step. constant production of recordable timer units with saturation (=reserve0)
          tau_release, ...             % release rate
          tau_activation, ...          % AE recharge time constant. how dense should the presynaptic spikes be to trigger the switch?
          evolution_noise, ...         % amount of noise in evolution (unit of timesteps/timestep) (=variance of Gaussian filter)
          max_isi, ...                 % maximum recordable time (in timesteps)
          reserve0,...
          trigger_on)                
                                      

%%%%%%%%%%%%%%% REPLENISH ACTIVATION ENERGY, UPDATE REFRACTORY TIMER %%%%%%%%%%%%%%%

activation_energy = activation_energy ...
                    + (activation_energy0-activation_energy)/tau_activation ... 
                    - CS_spike; 
record_refractory_timer = max(record_refractory_timer -1,0);
                
                
%%%%%%%%%%%%%%% EVOLVE CURRENT BATCH %%%%%%%%%%%%%%%           
% evolve batch (before or after saving??)
batch_writing = [0 batch_writing(1:end-1)]; %slide evolving recorder units to the 
%right  and add 0 in first slot indicating a recorder unit encoding
%time t
center = ceil(evolution_noise*2)+1; %mean of blurring kernel depends on ISI, but what is 1.5 for?
kernel = exp(-(((1:2*ceil(evolution_noise*2)+1)-center).^2)./evolution_noise); %Gaussian kernel on 1:noise*1.5
kernel = kernel/sum(kernel); %normalize
batch_writing = conv(batch_writing,kernel,'same'); %blur recorded units. is there a simpler way to do this?
                

%%%%%%%%%%%%%%% REPLENISH BATCH %%%%%%%%%%%%%%%      
% Fraction of recorder units left. The reserve is incremented by
% production_rate each time step until the maximum value, reserve0.
reserve = min(reserve + production_rate, reserve0); 


%%%%%%%%%%%%%%% SAVE %%%%%%%%%%%%%%%  
if US_spike == 1 && trigger_on
    % consolidate evolving timer units into recorded timer units
    batch_final = batch_writing; %stop recorder units for readout
    batch_writing = zeros(1,round(max_isi)); %reset recorder units
    trigger_on = 0;
    record_refractory_timer = record_refractory_timer0;
else
    batch_final = zeros(1,round(max_isi)); %set final batch to all zeros
end


%%%%%%%%%%%%%%% SWITCH ON/OFF %%%%%%%%%%%%%%%    
% release new timer units
if activation_energy < 0 && record_refractory_timer == 0
    trigger_on = 1;
    record_refractory_timer  = record_refractory_timer0;
end
if reserve/tau_release <= production_rate
    trigger_on = 0;
end



%%%%%%%%%%%%%%% RELEASE %%%%%%%%%%%%%%% 
if trigger_on
    newbatch = reserve/tau_release;  
    %The reserve ejects some fraction of its contents every ms

    batch_writing(1) = batch_writing(1) + newbatch;
    %The component of batch_writing indicating that the CS has
    %triggered recorder unit release is rewritten to indicate the
    %number of recorder units released at time t;

    reserve = reserve - newbatch; 
    %The reserve is depleted by the number of released recorder units
    
end
    
        
        

end
    
