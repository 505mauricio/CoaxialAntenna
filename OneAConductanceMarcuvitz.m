function [resp] = OneAConductanceMarcuvitz(f,a,b)
%Implementação da formula 1a referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook 
%Valido enquanto lambda > 2(a-b)/3.142
%f = 1*10^9:5800000:3*10^10;
%a = 12.5*10^-3/2;
%b = 4.47*10^-3/2;
%a/b = 2.7964
step = pi/500;
z = step:step:pi/2;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;
lambda = physconst('LightSpeed')./f;
k = 2*pi./lambda;
% Ja = besselj(0,k*a.*sin(z));
% Jb = besselj(0,k*b.*sin(z));
% resp = zeros(1,length(z-1));
% for i=1:length(z)
%     integral = sum((step./sin(z)).*(Ja(i)-Jb(i))^2);
%     resp(i) = (1/log(a/b))*integral;
% end
integral =0;
resp = zeros(1,length(z-1));
for j = 1:length(f)    
    for i=1:length(z)
        integral = integral+(step./sin(i)).*(besselj(0,(2*pi/lambda(j))*a.*sin(i))- besselj(0,(2*pi/lambda(j))*b.*sin(i)) ).^2;    
    end
    resp(j) = (1/log(a/b))*integral;
    integral = 0;
end

plot((a-b)./lambda,resp)
grid on
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
title('Formula 1a) referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook','interpreter','latex')
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)

% plot(z,Ja)
% grid on
% legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
% title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
% xlabel('z','interpreter','latex')
% ylabel('$J_\nu(z)$','interpreter','latex')

% plot(z,Jb)
% grid on
% legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
% title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
% xlabel('z','interpreter','latex')
% ylabel('$J_\nu(z)$','interpreter','latex')