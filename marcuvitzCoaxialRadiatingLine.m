function [G] = marcuvitzCoaxialRadiatingLine(f,a,b)
%Marcuvitz Equação 1b 4-16 
close all
lambda = 3*10^8./f;
c = (2/3)*(1/log(a/b));
d = (pi^2*(b^2-a^2)./(lambda.^2)).^2;
G= c*d;
eixo = (a-b)./lambda;
figure('Name','','NumberTitle','off')
plot (eixo,G)
title('Figura 4.16-2 Conductance of coaxial guide radiating into half space')
