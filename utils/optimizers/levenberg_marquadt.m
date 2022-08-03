function p = levenberg_marquadt(F, narg, tol, max_iter)
    p0 = rand(1, narg);
    f_last = F(p0);
    p = p0;
    m = norm(f_last);
    I = eye(narg);
    s = 1;
    mu = 1;
    for iter=1:max_iter
        if any(isnan(f_last))
            p0 = (rand(1, narg)-0.5)*100;
            f_last = F(p0);
            p = p0;
            m = norm(f_last);
            continue
        end
        mu = norm(f_last)*mu / m;
        J = Df(F, p0, 1e-6);
        inv_J = J * pinv(mu*I + (J')*J);
        p = p0 - s * inv_J * f_last;
        f = F(p);
        if norm(f-f_last) < tol
            break
        end
        p0 = p;
        f_last = f;
    end
end

