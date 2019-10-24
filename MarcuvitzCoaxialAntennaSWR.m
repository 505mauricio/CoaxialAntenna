function [SWR,G,B,fLimite] = MarcuvitzCoaxialAntennaSWR(fLimite,a,b)
%Função que plota SWR de antenas coaxiais segundo a formula do marcuvitz,
%para o intervalo de 1GHz até uma fLimite

[G,FreqConductance] = OneAConductanceMarcuvitz(fLimite,a,b);
[B,FreqSusceptance] = TwoASusceptanceMarcuvitz(fLimite,a,b);
Y = G+i*B;
r = (1 - Y)./(1+Y);
step = pi/500;
z = step:step:pi/2;
temp = (fLimite-10^9)/(length(z)-1);
f = 10^9:temp:fLimite;
lambda = physconst('LightSpeed')./f;
SWR = (1+abs(r))./(1-abs(r));

plot(f./10^9,SWR)
xlim([1 80])
ylim([1 3])
grid on
grid minor
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','best')

xlabel('Frequencia GHz')
ylabel('SWR','interpreter','latex','Rotation',0)




