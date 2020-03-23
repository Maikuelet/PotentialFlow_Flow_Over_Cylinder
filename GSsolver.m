%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  GAUSS-SEIDEL SOLVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stream,p,T,rho,v] = GSsolver(stream,coeff,nodeX,nodeY,sigma,meshSizes,p,T,rho,T0,p0,v0,fluidC,v,mat)

    

    N = meshSizes(1);
    M = meshSizes(2);
    
    eT   = 1; 
    ep   = 1;
    erho = 1;
    n = 0;
    
    
    while eT> sigma  || ep> sigma  || erho > sigma 
    
    % Initialize stream parameter;
    
    n = n+1    
    % Discretized Stream Func Solve
    % Solver Core
    for i = 2:N+1
        for j = 2:M+1
            
            if mat(j,i) == 0
            stream(j,i) = (coeff.ae(j,i)*stream(j,i+1) + ...
                           coeff.aw(j,i)*stream(j,i-1) + ...
                           coeff.an(j,i)*stream(j+1,i) + ...
                           coeff.as(j,i)*stream(j-1,i)) / coeff.ap(j,i);                     
            end
            
        end        
    end
    
    [stream] = neumannbc(stream,N);
    
    %Density Temperature and Pressure Solve
    % Field Solver
    Tini   = T;
    pini   = p;
    rhoini = rho;
    
    for i = 2:N+1
        for j = 2:M+1
            
            % Only compute at fluid media mat=0
            if mat(j,i) ==0
                
                v.vxn(j,i) = (stream(j+1,i)-stream(j,i) ) / (nodeY(j+1) - nodeY(j));            
                v.vye(j,i) = (stream(j,i+1)-stream(j,i) ) / (nodeX(i+1) - nodeX(i));

                % Due to the for dimensions it is necessary to keep face
                % velocity to zero, this if ensures that
                if j == M+1
                   v.vxn(j,i) = 0;
                   v.vye(j,i) = 0;
                end            

                vxP = (v.vxn(j,i)+v.vxn(j-1,i))/2;
                vyP = (v.vye(j,i)+v.vye(j,i-1))/2;

                % Due to the for dimensions the mean value is the point value
                if j == 2
                   vxP = v.vxn(j,i);
                   vyP = v.vye(j,i);
                end            
                if j == M+1
                   vxP = v.vxn(j-1,i);
                   vyP = v.vye(j,i-1);
                end

                v.vp(j,i) = sqrt(vxP^2 + vyP^2);                  

                % If -> Total energy conserved
                % Then:            
                T(j,i) = T0 + (v0^2-v.vp(j,i)^2)/(2*fluidC.cp);

                % Isentropic relation applyed
                p(j,i) = p0 * (T(j,i)/T0)^(fluidC.gam/(fluidC.gam-1));                
                
                % Density from gas relation
                % Air considered as an ideal gas
                rho(j,i) = p(j,i)/(fluidC.R*T(j,i));             
                
            end
                        
        end        
    end
    
    [p] = neumannbc(p,N);
    [rho] = neumannbc(rho,N);
    [T] = neumannbc(T,N);
    [v.vp] = neumannbc(v.vp,N);
    
    eT   = max(max(abs(T-Tini)));
    ep   = max(max(abs(p-pini)));
    erho = max(max(abs(rho-rhoini)));
    
    %[coeff] = interiorcoefficients(rho,nodeX,faceX,nodeY,faceY,rho0);
    
    
    end
    
    

end