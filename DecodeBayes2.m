classdef DecodeBayes2 < handle
   
   properties % Leave public, as get methods in other langauges are public
      
      sr, inPart
      kernel
      optTau,  tauAll,tau
      distMat, H, perf, mi
      distMatTau, HTau, perfTau, miTau,
      
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
      
           
      function optimize(self, tauAll, thetaAll)
         self.optimizeTauTheta(tauAll, thetaAll)      
      end
      
      function optimizeTauTheta(self, tauAll, thetaAll)      
         %optimizeTauTheta(tauAll,thetaAll)
         self.tauAll = tauAll;
         for tau = 1:length(tauAll)
            getKernel(self,self.tauAll(tau));
            self.sr.filterTrains(self.kernel);            
            convTrains = reshape(permute(self.sr.trains,[2,1,3]),...
               self.sr.cats*self.sr.trials,[]);
            self.sr.getRawTrains();
            rawTrains = self.srTrains();
            for theta = 1:length(thetaAll)
               embedVectors(convTrains);
               embedVectors(rawTrains);
               Htmp = clusterNearestMean(convTrains, rawTrains, self.inPart);
               self.HTau(tau,theta,:,:) = Htmp;
               self.miTau(tau,theta) = MI(Htmp);
               self.perfTau(tau,theta) = mean(diag(Htmp));
            end
         end
         [self.mi, optIdx] = max(self.miTau);
         self.optTau = self.tauAll(optIdx);
         self.perf = max(self.perfTau);
         self.H = squeeze(self.HTau(optIdx,:,:));
      end
      
      
      function plotOptimizeTau(self)
         subplot(122)
         plot(log2(self.tauAll),self.miTau,'k','LineWidth',1.5)
         hold on
         plot(log2(self.tauAll),3*self.perfTau,'r','LineWidth',1.5)
         xlabel('log_2(\tau)')
         axis('tight')
         set(gca,'YLim',[0 3])
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
