lambda = physconst('LightSpeed')./fCA400Susceptance;
plot((8.13/2*10^-3-2.74/2*10^-3)./lambda,CA400Susceptance)
hold on
grid on
grid minor
lambda = physconst('LightSpeed')./fSCF1250Susceptance;
plot((12.3/2*10^-3-3.56/2*10^-3)./lambda,SCF1250Susceptance)
lambda = physconst('LightSpeed')./fCA600Susceptance;
plot((12.5/2*10^-3-4.47/2*10^-3)./lambda,CA600Susceptance)
lambda = physconst('LightSpeed')./fCA900DBSusceptance;
plot((18.54/2*10^-3-6.6/2*10^-3)./lambda,CA900DBSusceptance)
lambda = physconst('LightSpeed')./fLCF7850DBSusceptance;
plot((25.2/2*10^-3-9.32/2*10^-3)./lambda,LCF7850DBSusceptance)
lambda = physconst('LightSpeed')./fLCF15850Susceptance;
plot((46.5/2*10^-3-17.6/2*10^-3)./lambda,LCF15850Susceptance)
hold off

legend({['CA400'],['SCF12-50'],['CA600'],['CA900DB'],['LCF78-50'],['LCF158-50']},'Interpreter','latex','Location','best')
xlim([0 0.9])
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{B}{Y_0}$','interpreter','latex','Rotation',0)