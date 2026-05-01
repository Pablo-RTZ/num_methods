function sol = IS(A,b)
% sol = IS(A,b)
% Implements inverse substitution to solve diagonal systems

b = b(:);
n = length(b);
sol = zeros(n,1);
sol(n) = b(n)./A(n,n);

for k=n-1:-1:1
    sol(k) = (b(k)-A(k,k+1:n)*sol(k+1:n))/A(k,k);
end

end