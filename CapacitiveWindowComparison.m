a = 9.32/2*10^-3;
c = 25.2/2*10^-3;
b = [c/2.5,c/2,c/1.5];
f = [1*10^9:1*10^8:10*10^9];
BMarcuvitz = zeros(length(b),length(f));
fLimiteCollin = physconst('LightSpeed')*0.1/(c-a);
fLimiteCollin = [1*10^9:1*10^8:fLimiteCollin];
BCollin = zeros(length(b),length(fLimiteCollin));
for i=1:length(f)
    BMarcuvitz(:,i) = CapacitiveWindowMarcuvitz(f(i),a,b,c);
end
for i=1:length(fLimiteCollin)
    BCollin(:,i) = CapacitiveWindowCollin(fLimiteCollin(i),a,b,c);
end

plot(f/10^9,BMarcuvitz(1,:))
hold on
plot(f/10^9,BMarcuvitz(2,:))
plot(f/10^9,BMarcuvitz(3,:))
plot(fLimiteCollin/10^9,BCollin(1,:))
plot(fLimiteCollin/10^9,BCollin(2,:))
plot(fLimiteCollin/10^9,BCollin(3,:))
hold off

legend({['M $\displaystyle\frac{c}{b}= $' num2str(c/b(1))],['M $\displaystyle\frac{c}{b}= $' num2str(c/b(2))]...
,['M $\displaystyle\frac{c}{b}= $' num2str(c/b(3))]...
,['C $\displaystyle\frac{c}{b}= $' num2str(c/b(1))],['C $\displaystyle\frac{c}{b}= $' num2str(c/b(2))]...
,['C $\displaystyle\frac{c}{b}= $' num2str(c/b(3))]},'Interpreter','latex','Location','Best')



title('Comparação entre janelas capacitivas Marcuvitz vs Collin para diversas Frequências','interpreter','latex')
%xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
xlabel('Frequência GHz')
ylabel('$\displaystyle\frac{B}{Y_0}$','interpreter','latex','Rotation',0)
