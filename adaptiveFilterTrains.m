function [newTrains, newPSTH, newP] = adaptiveFilterTrains(sr, param)
% mimic a dynamic syapse as in Houghton20xx
%
% USAGE
%  [newTrains, newPSTH] = adaptiveFilterTrains(sr, param)
% PARAMS
%  sr    - SpikeResp object
%  param - vector holding [tau1, tau2, phi]
% RETURNS
%  newTrains - filtered trains
%  newPSTH   - filtered PSTH (average over trials)
%  newP      - hidden variable "p" governing synaptic state 

% extract params
tau1 = param(1);
tau2 = param(2);
phi  = param(3);

dt = sr.binSize;     % get resolution of trains
buffer = ceil(4*max(tau1,tau2)/dt);
sr.getRawTrains();   % reset trains to binary
% pre-allocate for speed

newTrains = zeros(sr.cats, sr.trials, size(sr.trains,3) + buffer);
newPSTH = zeros(sr.cats, size(sr.trains,3) + buffer);
f = zeros(sr.trials, size(sr.trains,3) + buffer);
prePendTrain = zeros(sr.trials, buffer);
p = f+1;

% iterate over cats and times (loop is vectorized over trials)
for cat = 1:sr.cats
   catTrain = [squeeze(sr.trains(cat,:,:))];
   for t = 1:size(newTrains,3)-1
      f(:,t+1) = f(:,t) - f(:,t) * dt/tau1;                 % exponential decay of membrane voltage f 
      p(:,t+1) = p(:,t) + (1 - p(:,t)) * dt/tau2;             % p is the synapse's hidden state
      if t<size(catTrain,2) && any(catTrain(:,t)==1)        % SPIKE!!
         spikeIdx = catTrain(:,t)==1;                       % find trials where the spike occurs
         f(spikeIdx,t+1) = f(spikeIdx,t) + p(spikeIdx,t);   % increment f by p
         p(spikeIdx,t+1) = phi * p(spikeIdx,t+1);           % scale p by phi
         %f(spikeIdx,t+1) = b * f(spikeIdx,t) + 1;
      end
   end
   newTrains(cat,:,:) = f;                                  % put results into return vars
   newP(cat,:,:) = p;                                  % put results into return vars
   newPSTH(cat,:) = squeeze(mean(newTrains(cat,:,:),2));
end


