% called by diffusionReactionCentroDriver.m


function u = diffusionReactionCentro(xc,yc,bi,visc,dt,finalT,CENTRO,s)

    GAMMA = 1/visc;
    a = sqrt(GAMMA/(4*visc));
    b = sqrt(GAMMA*visc);
    c = a*(b-1);

    function v = rk4(V,t,k,F)               % 4th order Runge-Kutta
        s1 = feval(F,V,t);                                 
        s2 = feval(F,V + k*s1/2,t+k/2);                     
        s3 = feval(F,V + k*s2/2,t+k/2);                        
        s4 = feval(F,V + k*s3,t+k);                       
        v = V + k*(s1 + 2*s2 + 2*s3 + s4)/6;     
    end

    function ex = exact(x,y,t) 
        ex = 1./( 1 + exp( a*(x + y - b*t ) + c ) );
    end

    function fp = fStandard(V,t,dt)        % u_t = F(u), standard
        V(bi) = exact(xc(bi), yc(bi), t);
           fp = visc*DS*V + GAMMA*(V.^2).*(1 - V);              
    end

    function fp = fCentro(V,t,dt)    % u_t = F(u), centrosymmetry
        V(bi) = exact(xc(bi), yc(bi), t);
           fp = visc*rbfCentro.centroMult(V,L,M,2) + GAMMA*(V.^2).*(1 - V);              
    end

    safe = false;
    mu = 5e-15;
    phi = iqx();
    t = 0;
    U = exact(xc,yc,0);  % initial condition


    if CENTRO
        N = length(xc);
        tic
        [r, rx, ry] = phi.distanceMatrix2d(xc(1:N/2),yc(1:N/2),xc,yc);
        B = phi.rbf(r,s);
        H = phi.L(r, s);
        [kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu)

        D = rbfCentro.centroDM(B,H,N,2,mu,safe); 
        [L,M] = rbfCentro.centroDecomposeMatrix(D,2);
        
        while t<finalT
            u = rk4(U,t,dt,@fCentro);
            t = t + dt; 
            u(bi) = exact(xc(bi), yc(bi), t);  
            U = u;       
        end
        

        errorCentro = norm( u - exact(xc,yc,t), inf)
        toc

    else

        tic
        [r, rx, ry] = phi.distanceMatrix2d(xc,yc);
        B = phi.rbf(r,s);
        H = phi.L(r, s);
        kappaB = cond(B)
        DS = rbfx.dm(B,H,mu,safe);
        
        while t<finalT
            u = rk4(U,t,dt,@fStandard);
            t = t + dt; 
            u(bi) = exact(xc(bi), yc(bi), t);  
            U = u;       
        end
        errorStandard = norm( u - exact(xc,yc,t), inf)
        toc 
       
    end



end




