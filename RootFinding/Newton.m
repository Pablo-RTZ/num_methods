function [sol, dif, iter, ACOC]= Newton(f, df, x0, tol, maxiter)
% Add docstrings, nargout and translate code TODO

arguments

    f (1,1) function_handle
    df (1,1) function_handle
    x0 (1,1) double = 0
    tol (1,1) double = 1e-10
    maxiter (1,1) double = 50

end

iter = 0;
dif = tol + 1; %Garantiza entrada en el bucle
I = []; %Inicializa I

while and(dif > tol, iter < maxiter)
    x1= x0 - f(x0)/df(x0);
    dif = abs(x1 - x0);
    I = [I dif];
    x0 = x1;
    iter = iter + 1;
end

if iter >= maxiter %No ha convergido
    sol = [];
    ACOC = [];
    disp ("No se ha convergido a la solución")

else
    sol = x1;
    %ACOC = fACOC(I);
end

end