function [coord_list, I_list, bg] = Coordinate_Descent_2D(...
    coord_list, I_list, bg, im, height, width, sigma)
% PURPOSE:
% Gradient descent on 2D coordinate and weight refinement.
%---------------------------------------------------
% USAGE:
% [coord_list, I_list, bg] = Coordinate_Descent_2D(...
%   coord_list, I_list, bg, height, width, sigma)
%---------------------------------------------------
% INPUTS:
% coord_list:   list of input coordinate
% I_list:       list of input intensity
% bg:           background photon
% im:           acquired image
% height:       image height
% width:        image width
% sigma:        standard deviation of system PSF (pixel)
%---------------------------------------------------
% OUTPUTS:
% coord_list:   list of output coordinate
% I_list:       list of output intensity
% bg:           background photon
%---------------------------------------------------

%%
iter = 1;
max_iter = 200;

[xx, yy] = meshgrid(1 : width, 1 : height);
xx_r = xx(:); yy_r = yy(:);

A = Gen_Meas_From_2D_Coord(coord_list, height, width, sigma);

res_list = [];

%%
while iter < max_iter
    mu = A * I_list + bg;
    res_list(iter) = sum(mu - im .* log(mu + 1e-5));
    %% background descent
    grad_bg = sum((mu - im) ./ (mu + 1e-5));
    % backtracking
    alpha = 0.5;
    epsilon = 1;
    if grad_bg > 0
        epsilon = min(epsilon, bg / grad_bg);
    end
    df = -1e-6 * norm(grad_bg) ^ 2;
    mu = A * I_list + bg;
    fc = sum(mu - im .* log(mu + 1e-5));

    bg_temp = bg - epsilon * grad_bg;
    mu_base = A * I_list;
    mu = mu_base + bg_temp;

    while sum(mu - im .* log(mu + 1e-5)) > fc + epsilon * df
        epsilon = epsilon * max(alpha, 0.1);
        alpha = alpha * 0.8;
        if epsilon < 1e-10
            epsilon = 0;
            break;
        end
        bg_temp = bg - epsilon * grad_bg;
        mu = mu_base + bg_temp;
    end
    bg = bg - epsilon * grad_bg;
    
    %%
    for i = 1 : size(coord_list, 1)
        %% intensity descent
        temp_xx = erf((xx + 0.5 - coord_list(i, 1)) / sqrt(2) / sigma) -...
            erf((xx - 0.5 - coord_list(i, 1)) / sqrt(2) / sigma);
        temp_yy = erf((yy + 0.5 - coord_list(i, 2)) / sqrt(2) / sigma) -...
            erf((yy - 0.5 - coord_list(i, 2)) / sqrt(2) / sigma);
        
        grad_I = sum((mu - im) ./ (mu + 1e-5) .* temp_xx(:) .* temp_yy(:));
        
        % backtracking
        alpha = 0.5;
        epsilon = 1000;
        if grad_I > 0
            epsilon = min(epsilon, I_list(i) / grad_I);
        end
        df = -1e-6 * norm(grad_I) ^ 2;
        mu = A * I_list + bg;
        fc = sum(mu - im .* log(mu + 1e-5));
        
        I = I_list(i) - epsilon * grad_I;
        mu_base = A(:, [1 : i - 1, i + 1 : end]) * I_list([1 : i - 1, i + 1 : end], :);
        mu = mu_base + A(:, i) * I + bg;
        
        while sum(mu - im .* log(mu + 1e-5)) > fc + epsilon * df
            epsilon = epsilon * max(alpha, 0.1);
            alpha = alpha * 0.8;
            if epsilon < 1e-10
                epsilon = 0;
                break;
            end
            I = I_list(i) - epsilon * grad_I;
            mu = mu_base + A(:, i) * I + bg;
        end
        I_list(i) = I_list(i) - epsilon * grad_I;
        %% x y coord descent
        mu = A * I_list + bg;
        
        d_mu_x = I_list(i) / (2 * sqrt(2 * pi) * sigma) *...
            (exp(-(xx - 0.5 - coord_list(i, 1)) .^ 2 / (2 * sigma ^ 2)) -...
            exp(-(xx + 0.5 - coord_list(i, 1)) .^ 2 / (2 * sigma ^ 2))) .* temp_yy;
        d_mu_y = I_list(i) / (2 * sqrt(2 * pi) * sigma) * temp_xx .*...
            (exp(-(yy - 0.5 - coord_list(i, 2)) .^ 2 / (2 * sigma ^ 2)) -...
            exp(-(yy + 0.5 - coord_list(i, 2)) .^ 2 / (2 * sigma ^ 2)));
        
        grad_xy = [sum((mu - im) ./ (mu + 1e-5) .* d_mu_x(:)),...
            sum((mu - im) ./ (mu + 1e-5) .* d_mu_y(:))];
        
        % backtracking
        alpha = 0.5;
        epsilon = 1e-1;
        df = -1e-6 * norm(grad_xy) ^ 2;
        mu = A * I_list + bg;
        fc = sum(mu - im .* log(mu + 1e-5));
        
        coord = coord_list(i, :) - epsilon * grad_xy;
        mu_base = A(:, [1 : i - 1, i + 1 : end]) * I_list([1 : i - 1, i + 1 : end], :);
        meas_col = 1 / 4 *...
            (erf((xx_r + 0.5 - coord(1)) / sqrt(2) / sigma) -...
            erf((xx_r - 0.5 - coord(1)) / sqrt(2) / sigma)) .*...
            (erf((yy_r + 0.5 - coord(2)) / sqrt(2) / sigma) -...
            erf((yy_r - 0.5 - coord(2)) / sqrt(2) / sigma));
        mu = mu_base + meas_col * I_list(i) + bg;
        
        while sum(mu - im .* log(mu + 1e-5)) > fc + epsilon * df
            epsilon = epsilon * max(alpha, 0.1);
            alpha = alpha * 0.8;
            if epsilon < 1e-10
                epsilon = 0;
                meas_col = A(:, i);
                break;
            end
            coord = coord_list(i, :) - epsilon * grad_xy;
            meas_col = 1 / 4 *...
                (erf((xx_r + 0.5 - coord(1)) / sqrt(2) / sigma) -...
                erf((xx_r - 0.5 - coord(1)) / sqrt(2) / sigma)) .*...
                (erf((yy_r + 0.5 - coord(2)) / sqrt(2) / sigma) -...
                erf((yy_r - 0.5 - coord(2)) / sqrt(2) / sigma));
            mu = mu_base + meas_col * I_list(i) + bg;
        end
        A(:, i) = meas_col;
        coord_list(i, :) = coord_list(i, :) - epsilon * grad_xy;
    end
    
    if iter > 1 && res_list(iter - 1) - res_list(iter) < 1e-2
        break;
    end
    iter = iter + 1;
end