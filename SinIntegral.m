function [resp] = SinIntegral(x)

step = x/10^7;
range = step:step:x;
resp = sum(sin(range)./range)*step;

