function [phase, modulus2] = besselphase(nu, x)
    % BESSELPHASE Phase function for the Bessel function derivatives
    % This private function does not contain argument checks. Please use
    % production code found in specphase package. 
    
    J = besselj(nu,x);
    Y = bessely(nu,x);

    % Fix matlab bug which results in incorrect overflow  near 0
    Y(isinf(Y)) = -Inf;

    phase = atan2(Y, J);

    if x > nu
        approx = sqrt(x^2 - nu^2)-nu*asec(x/nu) - pi/4;
        phase = phase - round((phase-approx)/(2*pi))*2*pi;
    end

    if nargout > 1
        modulus2 = (J.^2+Y.^2); 
    end
end

