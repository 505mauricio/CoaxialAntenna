function [resp] = OneAConductanceMarcuvitz(f,a,b)
%Implementação da formula 1a referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook 
%Valido enquanto lambda > 2(a-b)/3.142
%f = 8*10^10;
%a = 12.5*10^-3/2;
%b = 4.47*10^-3/2;
%a/b = 2.7964
step = pi/500;
z = step:step:pi/2;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;
lambda = physconst('LightSpeed')./f;
k = 2*pi./lambda;
resp = zeros(1,length(z-1));
for j = 1:length(f)    
    integral = sum((step./sin(z)).*(besselj(0,(2*pi/lambda(j))*a.*sin(z))- besselj(0,(2*pi/lambda(j))*b.*sin(z)) ).^2);    
    resp(j) = (1/log(a/b))*integral;
end

plot((a-b)./lambda,resp)
grid on
grid minor
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
title('Formula 1a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)

