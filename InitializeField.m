%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  FIELD INITIAILIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stream,p,T,rho,v] = InitializeField(p0, T0, istream,rho0, nodeX,nodeY,v0,mat);
sizeX = numel(nodeX);
sizeY = numel(nodeY);

stream = zeros(sizeY,sizeX);
p      = zeros(sizeY,sizeX);
T      = zeros(sizeY,sizeX);
rho    = zeros(sizeY,sizeX);

stream(:,:) = istream;
p(:,:)      = p0;
T(:,:)      = T0;
rho(:,:)    = rho0;


%% Velocity Field
v.vp     = zeros(sizeY,sizeX);

%Face Velocity  Nfaces = Nnodes-1
v.vxn    = zeros(sizeY-1,sizeX-1);
v.vye    = zeros(sizeY-1,sizeX-1); %No velocity expected in Y direction

% Velocity Field
v.vp(:,:)     = v0;
v.vxn(:,:)    = v0;
%vye = 0

%% Obstacle filling 

for i = 2:sizeX-1
    for j =2:sizeY-1
        
        if mat(j,i) == 1
            rho(j,i)     = 0;
            stream(j,i)  = v0*(nodeY(sizeY)/2); 
            v.vp(j,i)    = 0;
        end
        
    end
end

end