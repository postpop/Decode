classdef DecodeBayes < handle
   
   properties % Leave public, as get methods in other langauges are public
      
      sr, inPart
      kernel
      optTau,  tauAll,tau
      H, perf, mi
      HCnt, perfCnt, miCnt
      HTau, perfTau, miTau,
      
   end
   
   
   methods (Access='public')
      function self = DecodeBayes(sr, inPart)
         %DecodeBayes(sr, inPart)
         if isa(sr,'SpikeResp')
            self.sr = sr;
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
         
         if nargin==1
            inPart = ceil((1:(sr.cats*sr.trials))/sr.trials);
         end
         self.inPart = inPart(:);
      end
      
      function getKernel(self,tau,whichFilter)
         self.tau=tau;
         if nargin==2 || whichFilter==0
            x =0:self.sr.binSize:3*tau;
            filter = x.*exp(-(x./tau).^2);%alpha
         elseif whichFilter==1
            x =0:self.sr.binSize:3*tau;
            filter = exp(-(x./tau));%exp
         else
            x =-3*tau:self.sr.binSize:3*tau;
            filter = exp(-(x./tau).^2);%gauss
         end
         filter = filter./norm(filter);
         self.kernel = filter;
      end
      
         
      function optimizeTau(self, tauAll,whichOpt,targetTrains)
         %optimizeTau(tauAll,[whichOpt=0 correct/ 1 is MI],[targetTrains])
         
         if nargin<=2
            whichOpt=0;
         end
         
         if nargin<4%cross validation
            targetTrains= self.sr.trains;
         end
         targetRates = sum(targetTrains,2);
         dataRates = sum(reshape(permute(self.sr.trains,[2,1,3]),self.sr.cats*self.sr.trials,[]),2);         
         self.HCnt = clusterNearestMean(dataRates, targetRates, self.inPart,1);
         self.perfCnt = mean(diag(self.HCnt));
         self.miCnt = MI(self.HCnt);
         
         
         self.tauAll = tauAll;
         for tau = 1:length(tauAll)            
            getKernel(self,self.tauAll(tau),3);
            [cnvTrains, rawTrains] = self.sr.filterTrains(self.kernel);
            if nargin==4%cross validation
               targetTrainsTmp = [zeros(size(targetTrains,1),floor(length(self.kernel)/2)),targetTrains,zeros(size(targetTrains,1),ceil(length(self.kernel)/2)-1)];
               if length(targetTrainsTmp)>length(cnvTrains)
                  lenDiff = length(targetTrainsTmp)- length(cnvTrains);
                  cnvTrains = [cnvTrains,zeros(size(targetTrainsTmp,1),lenDiff)];
               end
               Htmp = clusterNearestMean(cnvTrains, targetTrainsTmp, self.inPart);
            else%normal clustering
               Htmp = clusterNearestMean(cnvTrains, rawTrains, self.inPart);
            end
            self.HTau(tau,:,:) = Htmp;
            self.miTau(tau) = MI(Htmp);
            self.perfTau(tau) = mean(diag(Htmp));
         end
         [self.perf, optIdx] = max(self.perfTau);
         self.mi = max(self.miTau);
         if whichOpt==1
            [self.mi, optIdx] = max(self.miTau);
         end
         self.optTau = self.tauAll(optIdx);
         self.H = squeeze(self.HTau(optIdx,:,:));
      end
      
      
      function plotOptimizeTau(self)
         maxMI = log2(length(unique(self.inPart)));
         subplot(122)
         plot(log2(self.tauAll),self.miTau,'k','LineWidth',1.5)
         hold on
         plot(log2(self.tauAll),maxMI*self.perfTau,'r','LineWidth',1.5)
         xlabel('log_2(\tau)')
         axis('tight')
         set(gca,'YLim',[0 maxMI])
         legend({'MI [bit]','3*correct'}),legend('boxoff')
         
         subplot(121)
         myPcolor(self.H);
         title([self.perf,self.mi])
         set(gca,'CLim',[0 1])
         axis('square')
         colorbar
         drawnow
      end
      
   end
   
end
