classdef DecodeMetric < handle
   
   properties % Leave public, as get methods in other langauges are public
      
      sr
      inPart, nuisancePart
      kernel, binSize, tau
      distMat, clustMat, correct
      distMatTau, HTau, diagHTau
      optTau, perfTau, tauAll
      whichFilter, whichDist
      distMatLag
      tempProf
      
   end
   
   
   methods (Access='public')
      function self = DecodeMetric(sr, inPart)
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
%           
      function getDistMat(self, whichDist)
         if whichDist~=3
            self.sr.getRawTrains();
            self.sr.filterTrains(self.kernel)
            cnvTrains = reshape(permute(self.sr.trains,[2,1,3]),...
               self.sr.cats*self.sr.trials,[]);
         end
         
         if nargin==1 || whichDist==0
            % CALC distance matrix after vR - Euclidean Dist
            self.distMat = squareform(pdist(cnvTrains,'euclidean'));
         elseif nargin==2 && whichDist==1
            % CALC uncentered xcorr after Schreiber
            self.distMat = squareform(pdist(cnvTrains,'cosine'));
         elseif nargin==2 && whichDist==2
            % CALC xcorr after Schreiber
            self.distMat = squareform(pdist(cnvTrains,'correlation'));
         elseif nargin==2 && whichDist==3
            % CALC distance after VictorPurpura
            sStart =cumsum(size(self.sr.spikeTimes,3)*ones(1,self.sr.cats*self.sr.trials))-size(self.sr.spikeTimes,3);
            sEnd=sum(self.sr.spikeTimes>0,3)';sEnd=sEnd(:)';
            sEnd = sStart+sEnd;
            
            count = 0;
            for cat = 1:self.sr.cats
               for trials = 1:self.sr.trials
                  count = count+1;
                  spt(count,:) = self.sr.spikeTimes(cat,trials,:);
               end
            end
            spt = spt';
            self.distMat = reshape(spkdl(spt(:),sStart,sEnd,self.tau/25.4),self.sr.trials*self.sr.cats,[]);
         elseif nargin==2 && whichDist==4
            % proposal from Jan Hendrik, supposed to be phase invariant!!!
            self.distMat = zeros(self.sr.cats*self.sr.trials);
            hilbertTrains = hilbert(cnvTrains')';
            for tr1 = 1:size(hilbertTrains,1)
               for tr2 = tr1+1:size(hilbertTrains,1)
                  x = hilbertTrains(tr1,:);
                  y = hilbertTrains(tr2,:);
                  self.distMat(tr1,tr2) = abs((x./abs(x)) * (y./abs(y))');
                  self.distMat(tr2,tr1) = self.distMat(tr1,tr2);
               end
            end
            self.distMat = max(self.distMat(:))-self.distMat;
            for tr1 = 1:size(hilbertTrains,1)
               self.distMat(tr1,tr1) = 0;
            end
         end
         
      end
      
      function cluster(self)
         [self.clustMat] = clusterHeur(self.distMat, self.inPart);
      end
      
      function clusterCross(self,nuisancePart)
         [self.clustMat] = clusterHeurCross(self.distMat, self.inPart, nuisancePart);
      end
      
      
      function optimizeTau(self, tauAll, varargin)
         self.nuisancePart = [];
         if nargin==3
            self.nuisancePart = varargin{1};
         end
         
         whichDist=0;
         if length(varargin)==2
            whichDist= varargin{2};
         end
         
         cls = length(unique(self.inPart));
         self.tauAll = tauAll;
         self.distMatTau = zeros(length(tauAll), self.sr.cats*self.sr.trials,self.sr.cats*self.sr.trials);
         self.HTau = zeros(length(tauAll),cls,cls );
         self.diagHTau = zeros(length(tauAll), cls);
         fprintf('\n')
         for tau = 1:length(self.tauAll)
            % get distance matrix
            getKernel(self,self.tauAll(tau));
            self.getDistMat(whichDist);
            self.distMatTau(tau,:,:) = self.distMat;
            % cluster
            if isempty(self.nuisancePart)
               fprintf('.');
               self.cluster();
            else
               fprintf('x');
               self.clusterCross(self.nuisancePart);
            end
            self.HTau(tau,:,:) = self.clustMat;
            self.diagHTau(tau,:,:) = diag(self.clustMat);
         end
         fprintf('!\n')
         %% FIND OPTIMAL TAU and associated performance
         self.perfTau = mean(self.diagHTau,2);
         [optPerf optIdx] = max(self.perfTau);
         self.optTau = self.tauAll(optIdx);
         self.distMat = squeeze(self.distMatTau(optIdx,:,:));
         self.clustMat = squeeze(self.HTau(optIdx,:,:));
         self.tauAll = tauAll;
         
      end
      
      function plotOptimizeTau(self)
         
         trials = self.sr.trials;
         cats = self.sr.cats;
         
         subplot(3,6,[1 2 3 7 8 9])
         hold on
         myPcolor(self.distMat);
         vline(trials+1:trials:trials*cats,'k')
         hline(trials+1:trials:trials*cats,'k')
         hold off
         axis('square');
         axis('tight');
         set(gca,'YLim',[1 cats*trials+1],'XLim',[1 cats*trials+1])
         set(gca,'YTick' ,trials/2:trials:size(self.distMat,1)-trials/2,'YTickLabel' ,self.sr.catNames);
         set(gca,'XTick' ,trials/2:trials:size(self.distMat,1)-trials/2,'XTickLabel' ,self.sr.catNames);
         
         subplot(3,6,[4 5 6 10 11 12])
         myPcolor(self.clustMat);
         axis('square')
         set(gca,'YTick' ,1.5:1:cats+.5,'YTickLabel' ,self.sr.catNames);
         set(gca,'XTick' ,1.5:1:cats+.5,'XTickLabel' ,self.sr.catNames);
         set(gca,'CLim',[0 1]);
         colorbar()
         
         subplot(3,6,[13 14])
         semilogx(2.45*self.tauAll,self.perfTau);
         ylabel('\tau [ms]')
         title(['\tau_{opt}=' num2str(2.45*self.optTau,2)]);
         set(gca,'YLim',[0 1])
         axis('tight')
         subplot(3,6,[15 16])
         Z = linkage(squareform(self.distMat));
         dendrogram(Z);
         subplot(3,6,[17 18])
         try
            [Y] = mdscale(squareform(self.distMat),2);
            gscatter(Y(:,1),Y(:,2),self.inPart);
            legend('off')
         catch ME
            disp(ME)
         end
      end
      
      function tempProfile(self,D)
         %multi-dimensional scaling
         self.tempProf=[];
         try
            [C] = mdscale(squareform(self.distMat),D);
            %create Response matrix
            self.getKernel(self.optTau);
            R = [ones(self.sr.cats*self.sr.trials,1),self.convTrains];
            self.tempProf = C\R;
         catch ME
            disp(ME)
         end
      end
      
   end
   
   
end
