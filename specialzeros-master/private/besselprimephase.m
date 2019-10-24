function [phase, modulus2] =  besselprimephase(nu,x)
    % BESSELPRIMEPHASE Phase function for the Besse function derivatives
    % This private function does not contain argument checks. Please use
    % production code found in specphase package. 
   
    dJ = 0.5*(besselj(nu-1,x)-besselj(nu+1,x));
    dY = 0.5*(bessely(nu-1,x)-bessely(nu+1,x));

    % Fix matlab bug which results in incorrect overflow 
    % and problems with our evaluation
    dY(isnan(dY) || isinf(dY)) = Inf;

    phase = atan2(dY,dJ);
    if x > nu
        approx = sqrt(x^2 - nu^2)-nu*asec(x/nu) + pi/4;
        phase = phase - round((phase-approx)/(2*pi))*2*pi;
    end

    if nargout > 1
        modulus2 = (dJ.^2+dY.^2); 
    end
end
