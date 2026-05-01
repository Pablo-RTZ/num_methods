function sol = DS(A,b)
% sol = DS(A,b)
% Implements direct substitution to solve triangular systems

b = b(:);
n = length(b);
sol = zeros(n,1);
sol(1) = b(1)./A(1,1);

for k = 2:n
    sol(k) = (b(k)-A(k,1:k-1)*sol(1:k-1))/A(k,k);
end

end