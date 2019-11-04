function [SWR] = CalcSWR(f,a,bWindow,c,G,B,L)

BWindow = CapacitiveWindowMarcuvitz(f,a,bWindow,c);
lambda = physconst('LightSpeed')./f;
Y = G+i*B;
Y_l = (Y+j.*tan(2*pi./lambda*L))./(1+j*Y.*tan(2*pi./lambda*L));
Y_l = Y_l+j*BWindow;
r = (1 - Y_l)./(1+Y_l);
SWR = (1+abs(r))./(1-abs(r));


plot(f/10^9,SWR)
xlabel('Frequência GHz')
grid on
grid minor
xlim ([1 10])
ylim([1 5])
%set(gca, 'YScale', 'log')
ylabel('SWR','interpreter','latex','Rotation',0)
legend({['Raio do Disco Capacitivo = ' num2str(bWindow) ' Distando ' num2str(L) ]},'Interpreter','latex','Location','best')
title('Cellflex 7/8','interpreter','latex')



