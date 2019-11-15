function [resp,f] = TwoASusceptanceMarcuvitz(f,a,b)
%Implementação da formula 2a referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook 
%Valido enquanto lambda > 2(a-b)/gamma_1
%f = 8*10^10;
%a = 12.5*10^-3/2;
%b = 4.47*10^-3/2;
%a/b = 2.7964
step = pi/250;
z = step:step:pi;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;
lambda = physconst('LightSpeed')./f;
k = (2*pi./lambda)';
integral = sum(step*(2*sinint(k.*(a^2+b^2-2*a*b.*cos(z)).^0.5) - sinint(2*k*a.*sin(z./2)) -sinint(2*k*b.*sin(z./2))),2);    
resp = (1/(pi*log(a/b))).*integral;


plot((a-b)./lambda,resp)
grid on
grid minor
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
title('Formula 2a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{B}{Y_0}$','interpreter','latex','Rotation',0)

