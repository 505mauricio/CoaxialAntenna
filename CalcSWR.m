function [SWR] = CalcSWR(f,G,B,L,BWindow)
f = [1*10^9:10^7:f];
lambda = physconst('LightSpeed')./f;
Y = G+i*B;
Y_l = (Y+j*tan(2*pi./lambda*L))./(1+j*Y*tan(2*pi./lambda*L));
Y_l = Y_l+j*BWindow;
r = (1 - Y_l)./(1+Y_l);
SWR = (1+abs(r))./(1-abs(r));


plot(f/10^9,SWR)
xlabel('Frequência GHz')
grid on
grid minor
ylim([1 5])
ylabel('SWR','interpreter','latex','Rotation',0)




