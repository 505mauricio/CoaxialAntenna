x_1 = besscrosszero(0,(8.13/2*10^-3)/(2.74/2*10^-3),1);
gamma_1 = x_1*((8.13/2*10^-3)/(2.74/2*10^-3)-1)/pi; 
validadeCA400 = physconst('LightSpeed')*gamma_1/(2*(8.13/2*10^-3-2.74/2*10^-3));

x_1 = besscrosszero(0,(12.3/2*10^-3)/(3.56/2*10^-3),1);
gamma_1 = x_1*((12.3/2*10^-3)/(3.56/2*10^-3)-1)/pi; 
validadeSCF1250 = physconst('LightSpeed')*gamma_1/(2*(12.3/2*10^-3-3.56/2*10^-3));

x_1 = besscrosszero(0,(12.5/2*10^-3)/(4.47/2*10^-3),1);
gamma_1 = x_1*((12.5/2*10^-3)/(4.47/2*10^-3)-1)/pi; 
validadeCA600 = physconst('LightSpeed')*gamma_1/(2*(12.5/2*10^-3-4.47/2*10^-3));

x_1 = besscrosszero(0,(18.54/2*10^-3)/(6.6/2*10^-3),1);
gamma_1 = x_1*((18.54/2*10^-3)/(6.6/2*10^-3)-1)/pi; 
validadeCA900 = physconst('LightSpeed')*gamma_1/(2*(18.54/2*10^-3-6.6/2*10^-3));

x_1 = besscrosszero(0,(25.2/2*10^-3)/(9.32/2*10^-3),1);
gamma_1 = x_1*((25.2/2*10^-3)/(9.32/2*10^-3)-1)/pi; 
validadeLCF7850DB = physconst('LightSpeed')*gamma_1/(2*(25.2/2*10^-3-9.32/2*10^-3));

x_1 = besscrosszero(0,(46.5/2*10^-3)/(17.6/2*10^-3),1);
gamma_1 = x_1*((46.5/2*10^-3)/(17.6/2*10^-3)-1)/pi; 
validadeLCF15850 = physconst('LightSpeed')*gamma_1/(2*(46.5/2*10^-3-17.6/2*10^-3));