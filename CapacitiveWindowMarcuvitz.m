function [B] = CapacitiveWindowMarcuvitz(f,a,b,c)
% Janela capacitiva, Marcuvitz Waveguide Handbook pag 229
% a = Inner Conductor: 9.32/2*10^-3 Cellflex 7/8"
% c = Outer Conductor: 25.2/2*10^-3
%b = a+(c-a)/2
%x_1 = 3.097;%valor de 3.412 retirado da tabela na pag 74, x_01
%x_1 = 3.097/(c/a-1);
x_1 = besscrosszero(0,c/a,1);
gamma_1 = x_1*(c/a-1)/pi; 
lambda = physconst('LightSpeed')./f;
b_0 = c-a;
d = c-b;
A = 1./(1-(2*b_0./lambda).^2).^0.5-1;
A_1 = b./a*log(c./a)./(c./a-1).*((((c./b)-1)./(log(c./b))).^2);

A_2 = pi^2*(a./b)./(gamma_1.*(1-(2*b_0./(lambda.*gamma_1)).^2).^0.5)...
      .*((c./a)-1)./((besselj(0,x_1).^2./besselj(0,x_1*c/a).^2)-1)...
      .*(((besselj(0,x_1).*bessely(0,x_1*b./a)-bessely(0,x_1)*besselj(0,x_1*b./a))./((c./b)-1)).^2) ...
      -(1./(1-(2.*b_0./lambda).^2).^0.5 * (2./pi*b_0./d.*sin(pi.*d./b_0)).^2);


% A_2 = pi^2*a./b./(gamma_1*(1-(2*b_0./(gamma_1.*lambda)).^2).^0.5)/(gamma_1.*(1-(2.*b_0./(gamma_1.*lambda)).^2).^0.5).*(c/a-1)/(besselj(0,x_1)^2/besselj(0,x_1*c/a)^2-1)...
% *((besselj(0,x_1)*bessely(0,x_1*b/a)-bessely(0,x_1)*besselj(0,x_1*b/a))/(c/b-1))^2 ...
% -(2/pi*b_0/d*sin(pi*d/b_0))^2./(1-(2.*b_0./lambda).^2).^0.5;

B = 2*b_0./lambda.*A_1.*(4*log(csc(pi*d/(2*b_0)))+((4*A*cos(pi*d/(2*b_0)).^4))./((1+A*sin(pi*d/(2*b_0)).^4))...
+(((b_0./lambda).^2).*(1-3*sin(pi*d/(2*b_0)).^2).^2.*cos(pi*d/(2*b_0)).^4)+A_2);

