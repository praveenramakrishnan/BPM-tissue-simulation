function [E] = efield_focused_lambda_1300nm(X, Y, Z)
    lambda = 1300e-9;
    E = efield_plane(X, Y, Z, lambda);
end
