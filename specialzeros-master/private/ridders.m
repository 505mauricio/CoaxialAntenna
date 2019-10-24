function [d,n] = ridders(f, a, b, tol)
    % RIDDERS' root finding method. A quadratically convergent
    % bracketing method
    %
    %   d = RIDDERS(f, a, b, tol) solve f(d) == 0 for d in (a,b)
    %   to a relavtive tolerence of tol
    
    fa = f(a);
    fb = f(b);
    if fa*fb > 0
        error(['ridders: function values at the end points are the same sign;'  ...
               'the interval can not be guaranteed to bound a root'])
    end

    sgn = sign(fa);
    n = 0; 
    while abs((b-a)./max(abs([a,b]))) > tol
        c = (a+b)/2;
        fc = f(c);
        
        d = c + (b-c)*sgn*fc/sqrt(fc^2 - fa*fb);
        fd = f(d);
        if fd*fc < 0  % root between fd and fc
            if d < c
                a = d;
                fa = fd;
                b = c;
                fb = fc;
            else
                a = c;
                fa = fc;
                b = d;
                fb = fd;
            end
        elseif fa*fd < 0
            b = d;
            fb = fd;
        elseif fb*fd < 0
            a = d;
            fa = fd;
        else % fd==0
            return
        end
        n = n + 1;
    end
end
