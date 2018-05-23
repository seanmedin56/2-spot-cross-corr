% I think I took this from some version of the pipeline
function [intensity, fit] = gauss_integral(snip)
%This function takes a 2d pixelated spot and calculates the intensity by
%fitting the snippet to a 2d gaussian and subtracting out the background
%noise from the sum of the intesnities

    [mesh_x, mesh_y] = meshgrid(1:size(snip,1), 1:size(snip,2));
            
    %Single gaussian function
    singleGaussian = @(params) (params(1).*...
        exp(-(...
        (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / 2*params(5)^2))  .* (mesh_x-params(2)).^2 ...
        - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / 4*params(5)^2)) .* (mesh_x-params(2)).*(mesh_y-params(4))...
        + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / 2*params(5)^2)).* (mesh_y-params(4)).^2 ...
            )))...
        + params(6) - double(snip);

    x0 = [max(max(snip)), round(length(snip)/2), 1, round(length(snip)/2), ...
        1,min(min(snip)), 0];

    %Perform fitting
    lsqOptions=optimset('Display','none',...
    'maxiter',10000); 

    [fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        x0,zeros(1,7),inf(1,7), lsqOptions);
    
    
    intensity = sum(sum(double(snip) - fit(6)));

end

