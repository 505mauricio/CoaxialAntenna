a = 9.32/2*10^-3;
c = 25.2/2*10^-3;
b = [a:a/100:c/1.2];
f = 10*10^9;
BMarcuvitz = zeros(length(b),length(f));
for i=1:length(f)
    BMarcuvitz(:,i) = CapacitiveWindowMarcuvitz(f(i),a,b,c);
end

plot(c./b,BMarcuvitz)
grid on
grid minor
title('Valores de B normalizado para 10GHz segundo a equação do Marcuvitz','interpreter','latex')
xlabel('$\displaystyle\frac{c}{b}$','interpreter','latex')
ylabel('$\displaystyle\frac{B}{Y_0}$','interpreter','latex','Rotation',0)
