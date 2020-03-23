%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  INTERIOR COEFFICIENT   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coeff,rhoh] = interiorcoefficients(rho,nodeX,faceX,nodeY,faceY,rho0)
    
    sizeX = numel(nodeX);
    sizeY = numel(nodeY);
    
    coeff.ap = zeros(sizeY,sizeX);
    coeff.ae = zeros(sizeY,sizeX);
    coeff.aw = zeros(sizeY,sizeX);
    coeff.an = zeros(sizeY,sizeX);
    coeff.as = zeros(sizeY,sizeX);
    coeff.bp = zeros(sizeY,sizeX);
    
    
    % Interior nodes Filling
    for i = 2:sizeX-1       
        for j = 2:sizeY-1
            
            Dy = faceY(j)-faceY(j-1); 
            Dx = faceX(i)-faceX(i-1);
            
            dPE = nodeX(i+1)-nodeX(i);
            dPW = nodeX(i) - nodeX(i-1);
            dPN = nodeY(j+1) - nodeY(j);
            dPS = nodeY(j)  - nodeY(j-1);
            
            N = nodeY(j+1);
            S = nodeY(j-1);
            E = nodeX(i+1);
            W = nodeX(i-1);

            fn = faceY(j);
            fs = faceY(j-1);
            fe = faceX(i);
            fw = faceX(i-1);
            
            Px = nodeX(i);
            Py = nodeY(j);
            
            if rho(j,i) ~= 0
            
            rhohe = (dPE) / ( ((fe-Px)/(rho0/rho(j,i))) + ((E-fe)/(rho0/rho(j,i+1))) );
            rhohw = (dPW) / ( ((Px-fw)/(rho0/rho(j,i))) + ((fw-W)/(rho0/rho(j,i-1))) );
            rhohn = (dPN) / ( ((fn-Py)/(rho0/rho(j,i))) + ((N-fn)/(rho0/rho(j+1,i))) );
            rhohs = (dPS) / ( ((Py-fs)/(rho0/rho(j,i))) + ((fs-S)/(rho0/rho(j-1,i))) );
                
            coeff.ae(j,i) = rhohe * (Dy/dPE);
            coeff.aw(j,i) = rhohw * (Dy/dPW);
            coeff.an(j,i) = rhohn * (Dx/dPN);
            coeff.as(j,i) = rhohs * (Dx/dPS);
            coeff.ap(j,i) = coeff.ae(j,i)+coeff.aw(j,i)+coeff.an(j,i)+coeff.as(j,i);
            coeff.bp(j,i)= 0;
            
            end
                            
        end        
    end
    
    
end
