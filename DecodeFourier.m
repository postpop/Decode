classdef DecodeFourier < Decode2
    
    properties % Leave public, as get methods in other langauges are public
        spectraTrains;
        fitScape
        
        
    end
    
    
    methods (Access='public')
        
        function self = DecodeFourier(srN, inPart)
            % population decoding using vector space embedding after Houghton et al. 2008
            % with an nearest neighbour classifier
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
                self.distMat = squareform(pdist(self.cnvTrains,'cosine'));
            elseif nargin==2 && whichDist==2
                % CALC xcorr after Schreiber
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
        
        function optimizeTau(self, whichFreqs, thetaAll, cAll, sAll, whichDist, cutoffs)
            self.optimize(whichFreqs, thetaAll, cAll, sAll, whichDist, cutoffs);
        end
        
        function optimize(self, whichFreqs, thetaAll, cAll, sAll, whichDist, cutoffs)
            % optimize parameters of the embedding and the classifier
            %
            %ARGS
            % tauAll   - array of values to try for the alpha-filter width
            % thetaAll - array of values to try for the embedding angle
            %            between 0 and pi
            % cAll     - C parameter of the C-svm used
            % sAll     - sigma parameter of the RBF kernel
            % cutoffs  - cutoff frequency for spectrum as an equivalent to tau
            %            in the time domain metric
            cls = length(unique(self.inPart));
            
            tauAll = 1:length(cutoffs);
            
            self.tauAll = tauAll;
            self.thetaAll = thetaAll;
            self.distMatTau = zeros(length(tauAll),length(thetaAll), self.srN(1).cats*self.srN(1).trials,self.srN(1).cats*self.srN(1).trials);
            self.HTau = zeros(length(tauAll),length(thetaAll), cls,cls );
            self.diagHTau = zeros(length(tauAll), length(thetaAll), cls);
            self.srN(1).getSpec();
            
            for tau = 1:length(self.tauAll)
                % get distance matrix
                freqIdx = ismember(self.srN(1).specFreq, whichFreqs);
                freqIdx(freqIdx>cutoffs(tau)) = [];
                spectra = reshape(permute(self.srN(1).spec(:,:,freqIdx),[2,1,3]),...
                    self.srN(1).cats*self.srN(1).trials,[]);
                if whichDist==1%amplitude, set all phases
                    spectra = pol2cart(1, abs(spectra));
                end
                if whichDist==2%phase, scramble amplitudes
                    spectra = pol2cart(angle(spectra),1);
                end
                spectra = [real(spectra) imag(spectra)];
                for the = 1:length(thetaAll)
                    self.spectraTrains = embedVectors(spectra,thetaAll(the));
                    [accuracy, predLabel,H,fs] = ...
                        svmGrid(self.spectraTrains, self.inPart, sAll, cAll);
                    self.clustMat = squeeze(mean(H,1));
                    %svmGridLinear(self.spectraTrains, self.inPart, cAll);
                    self.cnvTrains = spectra;
                    %self.getDistMat(0);
                    %self.distMatTau(tau,the,:,:) = self.distMat;
                    %self.cluster();
                    
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
            
        end
        
        
        
    end
    
end
