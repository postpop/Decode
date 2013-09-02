classdef DecodeSVM2 < Decode2
   
   properties % Leave public, as get methods in other langauges are public
      fitScape      
   end
   
   
   methods (Access='public')
      
      function self = DecodeSVM2(srN, inPart)
         % population decoding using vector space embedding after Houghton et al. 2008
         % with an RBF SVM as classifier       
         %
         %ARGS:
         %  srN    - array of SpikeResp objects; each representing the
         %           responses of a different cell in a population to the
         %           same stimulus set
         %  inPart - [OPTIONAL] class membership of each trial, if omitted,
         %           infers inPart from cats given in srN
         if nargin==1
            inPart = ceil((1:(srN(1).cats*srN(1).trials))/srN(1).trials);
         end
         
         self = self@Decode2(srN, inPart);
      end
      
      function optimizeTau(self, tauAll, thetaAll, cAll, sAll)
         self.optimize(tauAll, thetaAll, cAll, sAll)
      end
      
      function optimize(self, tauAll, thetaAll, cAll, sAll)
         % optimize parameters of the embedding and the classifier
         %
         %ARGS
         % tauAll   - array of values to try for the alpha-filter width 
         % thetaAll - array of values to try for the embedding angle
         %            between 0 and pi
         % cAll     - C parameter of the C-svm used 
         % sAll     - sigma parameter of the RBF kernel
         
         cls = length(unique(self.inPart));
         self.tauAll = tauAll;
         self.thetaAll = thetaAll;
         self.distMatTau = zeros(length(tauAll),length(thetaAll), self.srN(1).cats*self.srN(1).trials,self.srN(1).cats*self.srN(1).trials);
         self.HTau = zeros(length(tauAll),length(thetaAll), cls,cls );
         self.diagHTau = zeros(length(tauAll), length(thetaAll), cls);
         for tau = 1:length(self.tauAll)
            % get distance matrix
            getKernel(self,self.tauAll(tau));
            cnv = zeros(self.srN(1).cats*self.srN(1).trials,length(self.srN),length(self.srN(1).trains));
            for s = 1:length(self.srN)
               self.srN(s).getRawTrains();
               self.srN(s).filterTrains(self.kernel)
               cnv(:,s,1:length(self.srN(s).trains)) = reshape(permute(self.srN(s).trains,[2,1,3]),...
                  self.srN(s).cats*self.srN(s).trials,[]);
            end
            
            for the = 1:length(thetaAll)
               self.cnvTrains = embedVectors(cnv,thetaAll(the));
               [accuracy, predLabel,H,fs] = ...
                  svmGrid(self.cnvTrains, self.inPart, sAll, cAll);             
                  %svmGridLinear(self.cnvTrains, self.inPart, cAll);
                  
               
               self.clustMat = squeeze(mean(H,1));
               self.HTau(tau,the,:,:) = self.clustMat;
               self.diagHTau(tau,the,:,:) = diag(self.clustMat);
               self.perfTau(tau,the) = mean(diag(self.clustMat));
               self.MITau(tau,the) = MI(self.clustMat);
               self.fitScape(tau,the,:,:) = squeeze(mean(fs))';
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
         
      end
      
   end
   
end
