function [k, itter] = besscrosszero(nu, l, N, varargin)
    % BESSCROSSZERO the roots of the Bessel function cross products
    %
    %     k = BESSCROSSZERO(m, labmda, N) the Nth root of the Bessel
    %     function cross products J(m,x)*Y(m,labmda*x) - Y(m,x)*J(m,labmda*x)
    %
    %     k = BESSCROSSZERO(m, labmda, N, BC) the Nth root of the function,
    %     when, BC = 'DD',
    %             J(m,x)*Y(m,labmda*x) - Y(m,x)*J(m,labmda*x)
    %           BC ='DN':
    %             Y'(m,x)*J(m,labmda*x) - J'(m,x)*Y(m,labmda*x)
    %           BC ='ND':
    %             J(m,x)*Y'(m,labmda*x) - Y(m,x)*J'(m,labmda*x)
    %           BC ='NN':
    %             J'(m,x)*Y'(m,labmda*x) - Y'(m,x)*J'(m,labmda*x)
    %
    %     N can be an array

    numvarargs = length(varargin);
    optargs = {'DD', eps}; %Defaults
    optargs(1:numvarargs) = varargin;
    [T,tol] = optargs{:};

    if numel(T) == 2 && ~regexp(T,'[ND]{2}')
        error('Type must be either ''D'' or ''N''.')
    end

    %Transform roots if necessary
    if l < 1
        if strcmp(T, 'DN')
            T = 'ND';
        elseif strcmp(T,'ND')
            T = 'DN';
        end
        k = zerobesscross(nu,1/l,N,T,tol)/l;
        return
    end

    if l == 1
        k = nan(size(N));
        return
    end

    i = (T(1) == 'N');
    j = (T(2) == 'N');

    delta = ((i==1).*(nu==0) - 1).*(j==1);
    l1 = l-1;
    mu = 4*nu^2;

    tp = theta(i,j,nu,nu,l);

    k = zeros(size(N));
    
    if nargout > 1
        itter = zeros(size(N));
    end
    
    % Initial guess for roots after the tp using McMahon's expansion
    mcm = ((N+delta)*pi >= tp);
    switch T
        case 'DD'
            k(mcm) = N(mcm)*pi/(l1);
            p = (mu-1)/(8*l);
            % q = (mu-1)*(mu-25)*(l^3-1)/(6*(4*l)^3*l1);
            % r = (mu-1)*(mu^2-114*mu+1073)*(l^5-1)/(5*(4*l)^5*l1);
        case 'DN'
            k(mcm) = (N(mcm)-0.5)*pi/(l1);
            p =-((mu+3)-(mu-1)*l)/(8*l*l1);
            % q =-((mu^2+46*mu-63)-(mu-1)*(mu-25)*l^3)./(6*(4*l)^3*l1);
            % r =-((mu^3+185*mu^2-2053*mu+1899) - (mu-1)*(mu^2-114*mu+1073)*l^5)/(5*(4*l)^5*l1);
        case 'ND'
            k(mcm) = (N(mcm)-0.5)*pi/(l1);
            p = ((mu+3)*l-(mu-1))/(8*l*l1);
            % q = ((mu^2+46*mu-63)*l^3-(mu-1)*(mu-25))./(6*(4*l)^3*l1);
            % r = ((mu^3+185*mu^2-2053*mu+1899)*l^5 - (mu-1)*(mu^2-114*mu+1073))/(5*(4*l)^5*l1);
        case 'NN'
            k(mcm) = (N(mcm)-1+(nu==0))*pi/(l1);
            p = (mu+3)/(8*l);
            % q = (mu^2+46*mu-63)*(l^3-1)./(6*(4*l)^3*l1);
            % r = (mu^3+185*mu^2-2053*mu+1899)*(l^5-1)/(5*(4*l)^5*l1);
    end
    if nu > 0 && l > 2
         k(mcm) = k(mcm) + p./k(mcm);
    end
    %McMahon's expansion for the roots/stationary points of the J_nu(x)
    k(~mcm) = (N(~mcm)+nu/2-(2*i-1)./4)*pi/l;

    for m = 1:length(N)
        n = N(m);
        x = k(m);
        err = 1;
        itt = 0;
        while err > tol && itt < 25
            %Newton's method
            [t,xdt] = theta(i,j,nu,x,l);
            update = (t - (n+delta)*pi)./(xdt);
            x = x*(1-update);
            err = abs(update);
            itt=itt+1;
            if nargout > 1
               itter(m) = itt;
            end
        end
        if itt >= 25
            warning('unable to reach tol at,nu = %f n = %d',nu,n)
        end
        k(m) = x;
    end
end

function [t,xdt] = theta(i,j,n,x,l)
    if i == 1
        [t1,M2] = besselprimephase(n,l*x);
        xdt1 = 2*(1-(n./(l*x)).^2)./(pi*M2);
    else
        [t1,M2] = besselphase(n,l*x);
        xdt1 = 2./(pi*M2);
    end

    if j == 1
        [t2,M2] = besselprimephase(n,x);
        xdt2 = 2*(1-(n./x).^2)./(pi*M2);
    else
        [t2,M2] = besselphase(n,x);
        xdt2 = 2./(pi*M2);
    end

    t = t1-t2;
    xdt = xdt1 - xdt2;
end
