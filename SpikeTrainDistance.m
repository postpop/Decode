classdef SpikeTrainDistance < handle

    properties (Constant=true)
        VR=1;
        VP=2;
        RC=3;
        FF=4;
        FA=5;
        FP=6;

    end
    
    properties (GetAccess='private',SetAccess='private')

    end
    properties

        sr
        metric
        distMat

    end

    methods (Access='public')
        function self = SpikeTrainDistance(varargin)
            if nargin==2 && findstr(class(varargin{1}),'SpikeResp')
                self.sr = varargin{1};
                self.metric = varargin{2};
            else
                self.sr = SpikeResp(varargin{1},varargin{2});
                self.metric = varargin{3};
            end

        end

        function distMat = getDistMatrix(self,param)
            if self.metric>3
                sr.getSpec();
                harmonicsFull = squeeze(reshape(...
                    sr.spec(:,:,ismember(sr.specFreq, param))...
                    ,cats*trials,1,[]));
            end


            switch self.metric
                case self.VR% van Rossum
                    distMat = squareform(pdist(...
                        convTrain(spikeTrains, getFilter(self,param,0) ),'euclidean'...
                        ));
                    
                %case self.VP% VictorPurpura

                case self.RC% RCORR
                    distMat = squareform(pdist(...
                        convTrain(spikeTrains, getFilter(self,param,0) ),'cosine'...
                        ));

                case self.FF% FOURIER full: amp+phase
                    distMat = pdist([real(harmonicsFull) imag(harmonicsFull)]);
                    
                case self.FA% FOURIER amp
                    distMat = pdist(abs([real(harmonicsFull) imag(harmonicsFull)]));
                    
                case self.FP% FOURIER phase
                    harmonicsFull = harmonicsFull./abs(harmonicsFull);
                    distMat = pdist(([real(harmonicsFull) imag(harmonicsFull)]));
            end

        end
        function filter = getFilter(self,tau,whichFilter)
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
        end

        function plotDistMatrix(self)
        end

    end

end