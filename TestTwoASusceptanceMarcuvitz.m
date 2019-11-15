function [resp,f] = TestTwoASusceptanceMarcuvitz(f,a,b)
%Implementação do codigo TwoASusceptanceMarcuvitz(f,a,b) usando função
%simplificada de sinint para ganhar desempenho
%Valido enquanto lambda > 2(a-b)/3.142
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
for i = 1:length(f)   
    integral = 0;
    for j = 1:length(z)
        integral(i) =integral + step*(2*SinIntegral(k(i).*(a^2+b^2-2*a*b.*cos(j)).^0.5) - SinIntegral(2*k(i)*a.*sin(j/2)) -SinIntegral(2*k(i)*b.*sin(j/2)));
    end
    resp(i) = (1/log(a/b))*integral(i);
end
%plot((a-b)./lambda,resp)
%grid on
%legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
%title('Formula 1a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
%xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
%ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)

