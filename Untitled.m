function SWRPlot(fLimite,a,b)
%Função que plota SWR de antenas coaxiais segundo a formula do marcuvitz,
%para o intervalo de 1GHz até uma fLimite
SWR1 = MarcuvitzCoaxialAntennaSWR(fLimite,a,b);

step = pi/500;
z = step:step:pi/2;
temp = (fLimite-10^9)/(length(z)-1);
fLimite = 10^9:temp:fLimite;

plot(fLimite./10^9,SWR1)
xlim([10 80])
ylim([1 5])
grid on
grid minor
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
xlabel('Frequência GHz','interpreter','latex')
ylabel('SWR','interpreter','latex','Rotation',0)
