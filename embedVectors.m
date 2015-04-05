function [embeddedVectors] = embedVectors(spikeTrains,theta)
%embeds vectors for computing the multi cell distance after Houghton 2007
%  embeddedVectors = embedVectors(spikeTrains,theta)
%ARGS:
%  spikeTrains - (cats*trials) x cels x time [ms]
%  theta       - "interaction" angle (usually between 0 and pi/2)
%RETURNS:
%  embeddedVectors         - embedded spiketrains

%%
trains = size(spikeTrains,1);
cels = size(spikeTrains,2);
trainLen = size(spikeTrains,3);

%% form rotation vectors for each cell
vec = zeros(cels);
vec(1,1) = 1;
cosTheta = cos(theta);
sinTheta = sin(theta);

if cels>=2
   for cel = 2:cels
      vec(cel,:) = vec(cel-1,:);
      vec(cel,cel-1) = vec(cel-1,cel-1) * cosTheta;
      vec(cel,cel) = vec(cel-1,cel-1) * sinTheta;
      [cosTheta sinTheta]  = getNewTheta(cosTheta);
      if isnan(sinTheta)
         sinTheta = 0;
      end
      if isnan(cosTheta)
         cosTheta = 0;
      end
   end
end


%% DEBUG
% ang = zeros(cels);
% for cel1 = 1:cels
%    for cel2 = 1:cels
%       x = vec(cel1,:);      y = vec(cel2,:);
%       ang(cel1,cel2) = acos(dot(x, y)/(norm(x)*norm(y)));
%    end
% end
% disp(theta/pi*180);
% disp(vec);disp('angles:');disp(unique(ang(:)));
%% EMBED vectors
embeddedVectors = zeros(trains,trainLen*cels);
for train1 = 1:trains
   t1 = reshape(spikeTrains(train1,:,:),cels,trainLen)'*vec;
   embeddedVectors(train1,:) = t1(:)';         
end
%clf;plot(embeddedVectors(1,:));title(theta),drawnow

%%
function [cosTheta, sinTheta] = getNewTheta(cosTheta)
cosTheta = (cosTheta - cosTheta^2) / (1 - cosTheta^2);
sinTheta = sqrt(1 - cosTheta^2);
%% EOF