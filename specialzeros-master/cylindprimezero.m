function x = cylindprimezero(nu, m, t)
    % x = CYLINDPRIMEZERO(nu, m, t)  The mth positive real root of the function
    %     J'(nu,x) cos(pi t) + Y'(nu,x) sin(pi t) == 0 
    % where J' and Y' are the derivatives of Bessel functions of the first and
    % second find.
    % 
    % m can be an array. nu and t must be scalars.
    %
    % The order nu must be postive.
    
    
    if length(nu) > 1
        error('cylindprimezero: the order must be a scalar')
    end
    if nu < 0
        error('cylindprimezero: the order must nonnegative')
    end

    if nargin < 3
        t = 0; 
    end

    % Which multiple of pi the first zero occurs at
    offset = 0;

    
    [phaseNu, Modulus2Nu] = besselprimephase(nu,nu);
    phaseTurningPoint = phaseNu/pi - 1/2;
    %Root less than nu
    
    if t < 0
        offset = -1;
        if t == phaseTurningPoint %Double root
            offset = 0;
        end
        if t < phaseTurningPoint
            offset = 1;
        end
    end 
    
    if nu == 0 && t == 0
        offset = 1;
    end
    
    x = zeros(size(m));

    for i=1:numel(m)
        targetPhase = (m(i)+offset+t-1/2)*pi;

        if m(i) == 1 
            if t == phaseTurningPoint && nu ~= 0
                x(i) = nu;
                continue
            end
            if offset == -1
                f =  @(x) besselprimephase(nu,x) - (t + 1/2)*pi 
                x(i)  = ridders(f, 0, nu, 2*eps);
                continue
            end
        end

        % First term in McMahon's expansion for root
        x(i) = targetPhase + (2*nu - 1)*pi/4;

        if m(i) == 2 && offset == -1
            x(i) = nu*(1 + pi*sqrt(Modulus2Nu/2*(t - phaseNu/pi+1/2)));
        end

        relerr = 1;
        n = 0;
        while abs(relerr) > eps && n < 20;
            [th,M2] = besselprimephase(nu,x(i));
            % Newton's method
            relerr = (th-targetPhase)*pi/2*M2/(1-(nu/x(i))^2);
            x(i) = x(i)*(1 - relerr);
            n = n + 1;
        end
        if n == 20
            warning(['cylindprimezero: Failed to attain tolerence ' ... 
                     'on root m(%d) = %d. Instead root has' ...
                     'relative error = %5g.'], ...
                     i, m(i), relerr)
        end
    end
end
