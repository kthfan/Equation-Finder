function J = Df(F, x, dx)
    len = numel(x);
    x1 = x - dx*eye(len);
    x2 = x + dx*eye(len);
    f1 = zeros(1, len);
    f2 = zeros(1, len);
    for i=1:len
        f1(i) = F(x1(i, :));
        f2(i) = F(x2(i, :));
    end
    J = ((f2 - f1) / (2*dx));
end