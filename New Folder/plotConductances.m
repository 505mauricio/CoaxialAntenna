lambda = physconst('LightSpeed')./fCA400Conductance;
[val,idx]=min(abs(5.484387197112775e+10-fSCF1250Conductance));
plot((8.13/2*10^-3-2.74/2*10^-3)./lambda(1:idx),CA400Conductance(1:idx))
hold on
grid on
grid minor

lambda = physconst('LightSpeed')./fSCF1250Conductance;
[val,idx]=min(abs(3.369152829627711e+10-fSCF1250Conductance));
plot((12.3/2*10^-3-3.56/2*10^-3)./lambda(1:idx),SCF1250Conductance(1:idx))


lambda = physconst('LightSpeed')./fCA600Conductance;
[val,idx]=min(abs(3.686485234919332e+10-fSCF1250Conductance));
plot((12.5/2*10^-3-4.47/2*10^-3)./lambda(1:idx),CA600Conductance(1:idx))


lambda = physconst('LightSpeed')./fCA900DBConductance;
[val,idx]=min(abs(2.479008562169565e+10-fSCF1250Conductance));
plot((18.54/2*10^-3-6.6/2*10^-3)./lambda(1:idx),CA900DBConductance(1:idx))


lambda = physconst('LightSpeed')./fLCF7850DBConductance;
[val,idx]=min(abs(1.865575338351915e+10-fSCF1250Conductance));
plot((25.2/2*10^-3-9.32/2*10^-3)./lambda(1:idx),LCF7850DBConductance(1:idx))


lambda = physconst('LightSpeed')./fLCF15850Conductance;
[val,idx]=min(abs(1.025629907711073e+10-fSCF1250Conductance));
plot((46.5/2*10^-3-17.6/2*10^-3)./lambda(1:idx),LCF15850Conductance(1:idx))
hold off

legend({['CA400'],['SCF12-50'],['CA600'],['CA900DB'],['LCF78-50'],['LCF158-50']},'Interpreter','latex','Location','best')
xlim([0 0.9])
xlabel('$\displaystyle\frac{a-b}{\lambda} $','interpreter','latex')
ylabel('$\displaystyle\frac{G}{Y_0}$','interpreter','latex','Rotation',0)