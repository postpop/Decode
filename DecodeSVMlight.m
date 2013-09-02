classdef DecodeSVMlight < handle
   
   properties % Leave public, as get methods in other langauges are public
      sr
      inPart
      kernel, binSize, tau
      clustMat, correct
      HTau, diagHTau, perfTau, MITau
      tauAll, cAll, sAll
      optTau, optPerf, optMI
      whichFilter
      
   end
   
   
   methods (Access='public')
      function self = DecodeSVMlight(sr, inPart)
         % Metric(sr, inPart)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.binSize = self.sr.binSize;
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
            x =0:self.binSize:3*tau;
            filter = x.*exp(-(x./tau).^2);%alpha
         elseif whichFilter==1
            x =0:self.binSize:3*tau;
            filter = exp(-(x./tau));%exp
         else
            x =-3*tau:self.binSize:3*tau;
            filter = exp(-(x./tau).^2);%gauss
         end
         filter = filter./norm(filter);
         self.kernel = filter;
      end
      
      
      function optimizeTau(self, tauAll, cAll, sAll)
         
         cls = length(unique(self.inPart));
         self.tauAll = tauAll;
         self.cAll = cAll;
         self.sAll = sAll;
         self.HTau = zeros(length(tauAll),cls,cls );
         self.diagHTau = zeros(length(tauAll), cls);
         fprintf('\n')
         for tau = 1:length(self.tauAll)
            % get distance matrix
            getKernel(self,self.tauAll(tau));
            self.sr.filterTrains(self.kernel)
            feats = reshape(permute(self.sr.trains,[2,1,3]),...
               self.sr.cats*self.sr.trials,[]);
            [accuracy, predLabel,H] = ...
               svmGridLIGHT(feats, self.inPart, sAll, cAll);
            
            self.HTau(tau,:,:) = squeeze(mean(H,1));
            self.MITau(tau) = MI(squeeze(mean(H,1)));
            self.diagHTau(tau,:) = diag(squeeze(mean(H,1)));
         end
         fprintf('!\n')
         %% FIND OPTIMAL TAU and associated performance
         self.perfTau = mean(self.diagHTau,2);
         [self.optMI optIdx] = max(self.MITau);
         self.optTau = self.tauAll(optIdx);
         self.optPerf = self.perfTau(optIdx);
         self.correct = self.optPerf;
         self.clustMat = squeeze(self.HTau(optIdx,:,:));
         self.tauAll = tauAll;
      end
      
      function plotOptimizeTau(self)
                  
         subplot(121)
         myPcolor(self.clustMat);
         axis('square')
         set(gca,'YTick' ,1.5:1:self.sr.cats+.5,'YTickLabel' ,self.sr.catNames);
         set(gca,'XTick' ,1.5:1:self.sr.cats+.5,'XTickLabel' ,self.sr.catNames);
         set(gca,'CLim',[0 1]);
         colorbar()
         title('%correct')
         subplot(122)
         
         [AX] = plotyy(log2(self.tauAll), self.perfTau,log2(self.tauAll), self.MITau);
         set(AX,'XLim',[min(self.tauAll) max(self.tauAll)])
         set(get(AX(1),'Ylabel'),'String','%correct')
         set(AX(1),'YLim',[0 1])
         set(get(AX(2),'Ylabel'),'String','MI')
         set(AX(2),   'YLim',[0 log2(length(unique(self.inPart)))])
         xlabel('log_2(\tau)')
         
      end
      
      
   end
   
   
end
