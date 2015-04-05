function [accuracy, predlabel,H,fitScape,decisionValues] = svmGrid(K0, labels, sAll, cAll)
%[accuracy, predlabel,H,fitScape,decisionValues] = svmGrid(K0, labels, sAll, cAll)

%% scale data to between 0 and 1 for RBF as this is optimal
K0 = (K0-min(K0(:)))/range(K0(:));
%%

pars = allcomb(sAll, cAll);
%acc = zeros(size(pars,1),1);
runs = 20;

cats = length(unique(labels));
trials = length(labels)/cats;
trialIdx = 1:trials;

allTest = zeros(runs,trials*cats);
for run = 1:runs%2 trials ea cat, test set
   [idx] = not(ismember(trialIdx,randsample(trialIdx,min(floor(trials*.9),trials-1))));
   allTest(run,:) = repmat(idx',[],cats);
end
allTrain = not(allTest);
pars = allcomb(sAll, cAll);

predlabel = zeros(size(allTrain,1),sum(allTest(run,:)));
accuracy = zeros(size(allTrain,1),3);
H = zeros(size(allTrain,1),cats,cats);
fitScape = zeros(runs,length(cAll),length(sAll));
decisionValues = zeros(runs,size(K0,1),10);
parfor run = 1:runs
   %determine optimal model parameters by 2fold xval
   acc = zeros(size(pars,1),1);

   for par = 1:length(pars)
      tmp = svmtrain(labels(allTrain(run,:)),K0(allTrain(run,:),:),['-s 0 -t 2 -b 0 -v 2 -m 1024 -g ' ...
         num2str(pars(par,1)) ' -c ' num2str(pars(par,2))  ]);
      if ~isempty(tmp)
         acc(par) = tmp;
      end
   end
   [tmp,par] = max(acc);
   fitScape(run,:,:) = reshape(acc,length(cAll),[]);   
   par = par(1);
   %train optimal model on full training set
   modelOpt = svmtrain(labels(allTrain(run,:)),K0(allTrain(run,:),:),['-s 0 -t 2 -b 0 -m 1024 -g ' ...
      num2str(pars(par,1)) ' -c ' num2str(pars(par,2))  ]);
   %optPar = pars(par,:);
   %test model on test data
   [predlabel(run,:), accuracy(run,:)] = svmpredict(labels(find(allTest(run,:))), K0(find(allTest(run,:)),:), modelOpt, '-b 0');
   %[dummy1, dummy2, decisionValues(run,:,1:length(modelOpt.rho))] = svmpredict(labels, K0, modelOpt, '-b 0');
   H(run,:,:) = confusionmat(predlabel(run,:),labels(find(allTest(run,:))))/sum(allTest(run,:))*length(unique(labels(find(allTest(run,:)))));
end


