Y = CA400Conductance+i*CA400Susceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(5.484387197112775e+10-fSCF1250Conductance));
plot(fCA400Conductance(1:idx)./10^9,SWR(1:idx))
hold on
grid on
grid minor

Y = SCF1250Conductance+i*SCF1250Susceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(3.369152829627711e+10-fSCF1250Conductance));
plot(fSCF1250Conductance(1:idx)./10^9,SWR(1:idx))


Y = CA600Conductance+i*CA600Susceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(3.686485234919332e+10-fSCF1250Conductance));
plot(fCA600Conductance(1:idx)./10^9,SWR(1:idx))

Y = CA900DBConductance+i*CA900DBSusceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(2.479008562169565e+10-fSCF1250Conductance));
plot(fCA900DBSusceptance(1:idx)./10^9,SWR(1:idx))

Y = LCF7850DBConductance+i*LCF7850DBSusceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(1.865575338351915e+10-fSCF1250Conductance));
plot(fLCF7850DBConductance(1:idx)./10^9,SWR(1:idx))

Y = LCF15850Conductance+i*LCF15850Susceptance;
r = (1 - Y)./(1+Y);
SWR = (1+abs(r))./(1-abs(r));
[val,idx]=min(abs(1.025629907711073e+10-fSCF1250Conductance));
plot(fLCF15850Conductance(1:idx)./10^9,SWR(1:idx))

hold off

ylim([1 5])
xlabel('Frequencia GHz')
ylabel('SWR','interpreter','latex','Rotation',0)
legend({['CA400'],['SCF12-50'],['CA600'],['CA900DB'],['LCF78-50'],['LCF158-50']},'Interpreter','latex','Location','best')