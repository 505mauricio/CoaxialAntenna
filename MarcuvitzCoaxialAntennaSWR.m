function [SWR] = MarcuvitzCoaxialAntennaSWR(f,a,b)
G = OneAConductanceMarcuvitz(f,a,b)';
B = TwoASusceptanceMarcuvitz(f,a,b);
Y = G+i*B;
r = (1 - Y)./(1+Y);
step = pi/500;
z = step:step:pi/2;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;
lambda = physconst('LightSpeed')./f;
SWR = (1+abs(r))./(1-abs(r));

plot((a-b)./lambda,G)
grid on
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
title('Formula 1a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)

plot((a-b)./lambda,B)
grid on
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
title('Formula 1a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)

plot(f./10^9,SWR)
xlim([10 80])
grid on
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
xlabel('Frequência GHz','interpreter','latex')
ylabel('SWR','interpreter','latex','Rotation',0)


