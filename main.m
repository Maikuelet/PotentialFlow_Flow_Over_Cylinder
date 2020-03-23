%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  ISENTROPIC POTENTIAL FLOW   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Miquel Altadill Llasat  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all
close all


InputData;

%% MESH GENERATION
[nodeX,faceX,nodeY,faceY] = UniformMesh(domainP, meshSizes);

%% MATERIAL
[mat] = material(nodeX,nodeY,radi);

%% FIELD INITIALIZATION
[stream,p,T,rho,v] = InitializeField(p0, T0, istream,rho0, nodeX,nodeY,v0,mat);

%% COEFFICIENTS
[coeff] = interiorcoefficients(rho,nodeX,faceX,nodeY,faceY,rho0);

%% BOUNDARY CONDITIONS
[stream,coeff,v] = BoundaryConditions(stream,v0,H,nodeY,meshSizes,coeff,v);

%% SOLVER
% - GS Based Solver
[stream,p,T,rho,v] = GSsolver(stream,coeff,nodeX,nodeY,sigma,meshSizes,p,T,rho,T0,p0,v0,fluidC,v,mat);

%% PostProcess
PostProcess(stream,p,T,rho,v);