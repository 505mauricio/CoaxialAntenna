f = 8*10^10;
a = 12.5*10^-3/2;
b = 4.47*10^-3/2;

SWR1 = MarcuvitzCoaxialAntennaSWR(f,a,b);
SWR2 = MarcuvitzCoaxialAntennaSWR(f,a,a/2);

step = pi/500;
z = step:step:pi/2;
temp = (f-10^9)/(length(z)-1);
f = 10^9:temp:f;

plot(f./10^9,SWR1)
xlim([10 80])
ylim([1 5])
grid on
grid minor
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')
xlabel('Frequência GHz','interpreter','latex')
ylabel('SWR','interpreter','latex','Rotation',0)

hold on

plot(f./10^9,SWR2)
legend({['$\displaystyle\frac{a}{b}= $' num2str(a/b)]},'Interpreter','latex','Location','Best')