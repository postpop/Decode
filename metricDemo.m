clear all
cats = 10; % different stimuli
trials = 10; % trials

% generates random spike trains - each response contains 20 spikes with
% stimulus specific offset
for cat = 1:cats
   spikeTimes(cat,:,:) = 10 + 5*cat + 2*cumsum(2+randn(1,trials,20),3);
end
resolution = 1;%ms
sr = SpikeResp(spikeTimes, resolution);
clf
subplot(221)
sr.plotDot();
title('spike times')
ylabel('stimulus')
subplot(222)
sr.filterTrains(2)
for cat = 1:cats
   hold on
   plot(sr.PSTH(cat,:)+cat*.4,'k');
end
xlabel('time [ms]')
title('PSTH')
%%
% inPart = ceil((1:(sr.cats*sr.trials))/sr.trials);
dm = DecodeMetric2(sr);
tauAll = 2.^([-2:2:16 ])% ms - 
thetaAll = 1;
dm.optimize(tauAll, thetaAll);
subplot(223)
imagesc(dm.clustMat)
colorbar
xlabel('decoded stimulus')
ylabel('actual stimulus')
title('confusion matrix')
axis('square')

subplot(224)
plot(dm.MITau,'.-k','MarkerSize',18)
xlabel('\tau [ms]')
ylabel('MI [bits]')