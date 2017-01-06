function J = crit_J(p, D)
    nu = p(3:end);
    theta = transpose([p(1) p(2)]);
    [m n] = size(D);
    [Col,Lin] = meshgrid(1:m,1:n);
    Gal = Sersic(nu,Lin,Col);

    H = [ones(m*n,1), Gal(:)];
    J = transpose(D(:) - H*theta)*(D(:) - H*theta);
end