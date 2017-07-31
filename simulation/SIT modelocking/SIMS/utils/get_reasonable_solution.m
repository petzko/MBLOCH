function X = get_reasonable_solution(options,s,T,g0,q0,kappa,gamma,tol,MAXREP)

X0 = [rand(1,2), rand(1,2), rand(1)];
[X,fval,exitflag,output] = fsolve('background_stability_III',X0,options,s,T,g0,q0,kappa,gamma);
ctr = 1;
while ((sum(abs(imag(X))) > tol) || (X(end))<tol || exitflag <= 0) && (ctr < MAXREP)
    X0 = [rand(1,2), rand(1,2), rand(1)];
    [X,fval,exitflag,output] = fsolve('background_stability_III',X0,options,s,T,g0,q0,kappa,gamma);
    if mod(ctr,50) == 0
        display(['repsol: ' num2str(ctr) ]);
    end
    ctr = ctr + 1;
end
X = real(X);

end