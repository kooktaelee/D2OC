function [PDF, array_x1, array_x2] = gauss2d(array_x1,array_x2,mean,covar_matrix)
% Generate Gaussian Probability Distribution Function
% Input
% array_x1 = discretisized array of x1 (1xn array)
% array_x2 = discretisized array of x2 (1xn array)
% mean = mean value of PDF along the x1 axis (1x2 array)
% covar_matrix = convariance matrix of PDF (2x2 matrix)
% Output
% PDF = Gaussian Probability Distribution Function (nxn matrix)

sigma_1 = sqrt(covar_matrix(1,1));
sigma_2 = sqrt(covar_matrix(2,2));
const_1A = 1/(sqrt(2*pi)*sigma_1);
const_1B = 1/(2*(sigma_1^2));
const_2A = 1/(sqrt(2*pi)*sigma_2);
const_2B = 1/(2*(sigma_2^2));
for cnt_1 = 1:length(array_x1)
    for cnt_2 = 1:length(array_x2)
    PDF(cnt_2,cnt_1) = const_1A*exp(-const_1B*((array_x1(cnt_1)-mean(1))^2))...
        *const_2A*exp(-const_2B*((array_x2(cnt_2)-mean(2))^2));
    end
end

end