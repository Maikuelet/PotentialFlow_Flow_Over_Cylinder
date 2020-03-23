
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  MESH GENERATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nodeX,faceX,nodeY,faceY] = UniformMesh(domainP, meshSizes)

xLength = domainP([1],[2])-domainP([1],[1]);  
yLength = domainP([2],[2])-domainP([2],[1]);                              
domainLengths=[xLength, yLength]; 

%% X AXIS
dim=1;
[nodeX,faceX]=facesZVB(domainLengths(dim),...
          meshSizes(dim),domainP([dim],[1]));



%% Y AXIS
dim=2;
[nodeY,faceY]=facesZVB(domainLengths(dim),...
  meshSizes(dim),domainP([dim],[1]));

%nodeY = flip(nodeY);
%faceY = flip(faceY);

end

function [nx,fx]=facesZVB(length,numCV,initPoint)

fx=linspace(initPoint,initPoint+length,numCV+1);
nx(1,1)=initPoint;
nx(1,2:numCV+1)=(fx(2:end)+fx(1:end-1))*0.5;
nx(1,numCV+2)=initPoint+length;

end
