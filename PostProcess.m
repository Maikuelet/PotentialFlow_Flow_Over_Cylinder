%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  POSTPROCESS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PostProcess(stream,p,T,rho,v)

figure(1);
h = heatmap(stream);  %heatmap
colormap(jet);
colorbar;
title('Stream Function');
xlabel('x'), ylabel('y');



figure(2);
contour(stream,20);
colormap(jet);
colorbar;
title('Stream Function Streamlines');
xlabel('x'), ylabel('y');


figure(3);
h = heatmap(p);  %heatmap
colormap(jet);
colorbar;
title('Pressure Field');
xlabel('x'), ylabel('y');


figure(4);
h = heatmap(T);  %heatmap
colormap(jet);
colorbar;
title('Temperature Field');
xlabel('x'), ylabel('y');


figure(5);
h = heatmap(rho);  %heatmap
colormap(jet);
colorbar;
title('Density Field');
xlabel('x'), ylabel('y');

figure(6)
h = heatmap(v.vp);
colormap(jet);
colorbar;
title('Velocity Field');
xlabel('x'), ylabel('y');

end