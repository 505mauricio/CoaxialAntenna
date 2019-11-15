freq = 20*10^9;
lambda = 3*10^8/freq;
beta = 2*pi/lambda;
L = 0.02*lambda;
Ya = 0.8123 + j*0.9096;
Ya_L = (Ya+j*tan(beta*L))/(1+j*Ya*tan(beta*L))