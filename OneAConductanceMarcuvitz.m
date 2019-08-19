function [resp] = OneAConductanceMarcuvitz(f,a,b)
%Implementação da formula 1a referente a parte 4-16 do livro Marcuvitz-Waveguide Handbook 
%f = 1*10^9:5800000:3*10^10;
%a = 8.13*10^-3/2;
%b = 2.74*10^-3/2;
step = pi/10000;
z = step:step:pi/2;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;
lambda = physconst('LightSpeed')./f;
k = 2*pi./lambda;
Ja = besselj(0,k.*a.*sin(z));
Jb = besselj(0,k.*b.*sin(z));
for i=1:length(z)
    integral = sum((step./sin(z)).*(Ja(i)-Jb(i)).^2);
    resp(i) = (1/log(a/b)).*integral;
end

plot((a-b)./lambda,resp)
grid on
legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')