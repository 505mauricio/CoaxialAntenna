function [B] = CapacitiveWindowCollin(f,a,b,c)
%Valido para c-a<=0.1*lambda_0
%Referente a equação 8.3 da pagina 552 do livro "Foundations For Microwave
%Engineering
%No livro do collin "b" é o raio do condutor externo, "a" o raio do interno
%e "c" o raio da estrutura metalica adicionada. Para manter a consistencia
%os nomes das variáveis foram trocados para serem coerentes com os nomes do
%codigo do modelo do Marcuvitz

lambda_0 = physconst('LightSpeed')./f;

B = 8*(c-b).^2./(lambda_0*b)...
    .*log(c/a)./(log(c./b).^2)...
    .*log(csc(pi/2*((c-b)./(c-a))));