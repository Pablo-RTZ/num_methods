function Cn = Cheby(n)
% Cn = Cheby(n)
% Returns coefficients for the n-th Chebyshev polynomial

tT0=1; Cn=tT0; degree=1;

while degree<=n
    if degree==1
        tT1=[1 0];
        Cn=tT1;
        degree=degree+1;
    else
        if degree>=2
            tT2=2*conv([1 0],tT1)-[0 0 tT0];
            Cn=tT2;
            tT0=tT1;
            tT1=tT2;
            degree=degree+1;
        end
    end
end
end