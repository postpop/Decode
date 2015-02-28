classdef DecodeMetric2 < Decode2
   
   properties % Leave public, as get methods in other langauges are public
      biasMI
      biasClustMat
      biasDistMat
   end
   
   
   methods (Access='public')
      
      function self = DecodeMetric2(srN, inPart)
         % population decoding using vector space embedding after Houghton et al. 2008
         % with a nearest neighbour classifier
         %
         %ARGS:
         %  srN    - array of SpikeResp objects; each representing the
         %           responses of a different cell in a population to the
         %           same stimulus set
         %  inPart - class membership of each trial
         if nargin==1
            inPart = ceil((1:(srN(1).cats*srN(1).trials))/srN(1).trials);
         end
         
         self = self@Decode2(srN, inPart);
      end
      
      
      function getDistMat(self, whichDist)
         
         if nargin==1 || whichDist==0
            % CALC distance matrix after vR - Euclidean Dist
            self.distMat = squareform(pdist(self.cnvTrains,'euclidean'));
         elseif nargin==2 && whichDist==1
            % CALC uncentered xcorr after Schreiber
            disp('schreiber')
            self.distMat = squareform(pdist(self.cnvTrains,'cosine'));
         elseif nargin==2 && whichDist==2
            % CALC xcorr 
            disp('pearson')
            self.distMat = squareform(pdist(self.cnvTrains,'correlation'));
         elseif nargin==2 && whichDist==3
            % CALC distance after VictorPurpura
            disp('not implemented yet...')
         end
         
      end
      
      function cluster(self)
         [self.clustMat] = clusterHeur(self.distMat, self.inPart);
      end
      
      function clusterCross(self,nuisancePart)
         [self.clustMat] = clusterHeurCross(self.distMat, self.inPart, nuisancePart);
      end
      
      function optimizeTau(self, tauAll, thetaAll, nuisancePart, whichDist)
         self.optimize(tauAll, thetaAll, nuisancePart, whichDist)
      end
      
      function optimize(self, tauAll, thetaAll, nuisancePart, whichDist)
         % optimize parameters of the embedding and the classifier
         %
         %ARGS
         % tauAll   - array or matrix of values to try for the alpha-filter
         %            width. if matrix than the second dim must equal the
         %            number of cells. this makes the algorithm chose tau
         %            for each cell independently.
         % thetaAll - array of values to try for the embedding angle
         %            between 0 and pi
         % whichDist - 0 - euclidean, 1 - schreiber, 2 - Pearson's coef of
         %              corr
         % nuisancePart - 
         self.nuisancePart = [];         
         if nargin>3
            self.nuisancePart = nuisancePart;
         end
         self.whichDist = 0;
         if nargin==5
            self.whichDist = whichDist;           
         end
         % get all tau combinations if needed
         if numel(tauAll)>length(tauAll) %each cell its own tau
            % get all combinations of taus
            cmdString = ['allcomb(' mat2str(tauAll(1,:))];
            for s = 2:length(self.srN)
               cmdString = [cmdString ', ' mat2str(tauAll(s,:)')];
            end
            cmdString = [cmdString ')'];
            tauAllCmb = eval(cmdString);
         else %one tau for all
            tauAllCmb = tauAll(:);
         end
         
         cls = length(unique(self.inPart));
         self.tauAll = tauAllCmb;
         self.thetaAll = thetaAll;
         self.distMatTau = zeros(length(tauAll),length(thetaAll), self.srN(1).cats*self.srN(1).trials,self.srN(1).cats*self.srN(1).trials);
         self.HTau = zeros(length(tauAll),length(thetaAll), cls,cls );
         self.diagHTau = zeros(length(tauAll), length(thetaAll), cls);
         srNtmp = self.srN;
         for tau = 1:size(self.tauAll,1)
            % embed spike trains
            cnv = zeros(self.srN(1).cats*self.srN(1).trials,length(self.srN),length(self.srN(1).trains));%preallocate
                       
            for s = 1:length(self.srN)
               if size(self.tauAll,2)>1 %each cell its own tau
                  thisTau = self.tauAll(tau,s);
               else %on tau for all
               	thisTau = self.tauAll(tau);
               end    
%                % ensure that there are no more than 100bins per tau
%                if 100*self.srN(s).binSize<min(thisTau)
%                   newBinSize = 100*self.srN(s).binSize;      
%                   self.srN(s) = SpikeResp(self.srN(s).spikeTimes, newBinSize);               
%                   self.binSize = newBinSize;                  
%                end
               getKernel(self,thisTau);
               self.srN(s).filterTrains(self.kernel,'full')
               cnv(:,s,1:length(self.srN(s).trains)) = reshape(permute(self.srN(s).trains,[2,1,3]),...
                                                         self.srN(s).cats*self.srN(s).trials,[]);
            end
            
            for the = 1:length(thetaAll)
               self.cnvTrains = embedVectors(cnv,thetaAll(the));
               self.getDistMat(self.whichDist);
               self.distMatTau(tau,the,:,:) = self.distMat;
               % cluster
               if isempty(self.nuisancePart)
                  %fprintf('.');
                  self.cluster();
               else
                  %fprintf('x');
                  self.clusterCross(self.nuisancePart);
               end
               self.HTau(tau,the,:,:) = self.clustMat;
               self.diagHTau(tau,the,:,:) = diag(self.clustMat);
               self.perfTau(tau,the) = mean(diag(self.clustMat));
               self.MITau(tau,the) = MI(self.clustMat);
            end
         end
         %% FIND OPTIMAL TAU and associated performance
         [self.optMI,idxMI] = max(self.MITau(:));
         [i,j] = ind2sub(size(self.MITau),idxMI);
         self.optPerf = self.perfTau(i,j);
         self.optTau = self.tauAll(i);
         self.optThe = self.thetaAll(j);
         self.clustMat = squeeze(self.HTau(i,j,:,:));
         self.distMat = squeeze(self.distMatTau(i,j,:,:));
         
         optH = self.clustMat;
         optH(optH==0) = optH(optH==0)+eps;
         self.entropy = -nansum(optH.*log2(optH),2);
         %% reset
         self.kernel = [];
         self.cnvTrains = [];
         self.srN = srNtmp;
         for s = 1:length(self.srN)
            self.srN(s).getRawTrains();
         end
         self.binSize = self.srN(s).binSize;
      end
      
      function [biasMI] = getBias(self, varargin)         

         if nargin>1
            thisTau = varargin{1};
            thisThe = varargin{2};
         else
            thisTau = self.optTau;
            thisThe = self.optThe;
         end
         
         runs = 10;
         cls = length(unique(self.inPart));
         self.biasMI = zeros(runs,1);
         self.biasClustMat = zeros(runs,cls,cls);
         self.biasDistMat = zeros(runs,self.srN(1).cats*self.srN(1).trials,self.srN(1).cats*self.srN(1).trials);
         for run = 1:runs
            dmBias = DecodeMetric2(self.srN, self.inPart(randperm(length(self.inPart))));
            dmBias.optimize(thisTau, thisThe);
            self.biasMI(run) = dmBias.optMI;
            self.biasClustMat(run,:,:) = dmBias.clustMat;
            self.biasDistMat(run,:,:) = dmBias.distMat;
         end
         biasMI = mean(self.biasMI);
      end
   end
   
end
