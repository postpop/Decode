classdef Decode2 < handle
   
   properties % Leave public, as get methods in other langauges are public
      
      srN
      inPart, nuisancePart
      kernel, binSize, tau
      distMat, clustMat, correct
      entropy
      distMatTau, HTau, diagHTau
      optTau, optThe, optMI, optPerf, perfTau, MITau
      tauAll, thetaAll
      whichFilter, whichDist
      cnvTrains
      
      sameDist, diffDist, sameDistDiff
   end
   
   methods (Abstract)
      optimizeTau(self);
   end
   
   methods (Access='public')
      
      function self = Decode2(srN, inPart)
         % Metric(sr1, inPart)
         if isa(srN,'SpikeResp')
            self.srN = srN;
            self.binSize = self.srN(1).binSize;
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
         
         if nargin==1
            inPart = ceil((1:(srN(1).cats*srN(1).trials))/srN(1).trials);
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
      
      function [sameDist, diffDist, sameDistDiff] = getDistDist(self)
         %[sameDist, diffDist, sameDistDiff] = getDistDist()
         
         
         cats = self.srN(1).cats;
         trials = self.srN(1).trials;
         for cat1 = 1:cats
            tmp = self.distMat( (cat1-1)*trials + (1:trials), (cat1-1)*trials + (1:trials));
            self.sameDist(cat1,:) = squareform(tmp);
         end
         
         cnt = 0;
         for cat1 = 1:cats
            for cat2 = cat1+1:cats
               tmp = self.distMat( (cat1-1)*trials + (1:trials), (cat2-1)*trials + (1:trials));
               cnt = cnt + 1;
               self.sameDistDiff(cnt,:) = [self.sameDist(cat1,:) self.sameDist(cat2,:)];
               self.diffDist(cnt,:) = tmp(:);
            end
         end
         sameDist = self.sameDist;
         diffDist = self.diffDist;
         sameDistDiff = self.sameDistDiff;
         
      end
      
      function plotOptimizeTau(self)
         
         subplot(221)
         myPcolor(num2cellstr(self.thetaAll/pi*180,2),self.tauAll,self.perfTau)
         set(gca,'CLim',[0 1])
         title('%correct')
         
         subplot(222)
         myPcolor(num2cellstr(self.thetaAll/pi*180,2),self.tauAll,self.MITau)
         hold on, plot(find(self.thetaAll==self.optThe)+.5,find(self.tauAll==self.optTau)+.5,'ok');
         set(gca,'CLim',[0 log2(length(unique(self.inPart)))])
         title('MI')
         
         subplot(223)
         myPcolor(self.clustMat)
         set(gca,'CLim',[0 1])
         title('optimal H')
         
         subplot(224)
         [AX] = plotyy(1:self.srN(1).cats, diag(self.clustMat),1:self.srN(1).cats,log2(self.srN(1).cats)-self.entropy);
         set(AX,'XLim',[1 self.srN(1).cats])
         set(get(AX(1),'Ylabel'),'String','%correct')
         set(AX(1),'YLim',[0 1])
         set(get(AX(2),'Ylabel'),'String','MI')
         set(AX(2),   'YLim',[0 log2(length(unique(self.inPart)))])
         xlabel('stim#')
      end
      
   end
   
end
