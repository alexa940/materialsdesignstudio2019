function [ K, cholL, s] = GP_Kernel_2D( X, L, sf, sn )
%GP_Kernel Code to calculate the GP Kernel (Squared Exponential)
%   Uses the X, L and sf values to calculate the Squared Exponential Kernel
%   matrix and its inversion for use in a GP regression model.

[kd, dim] = size(X);

if length(sn) > 1;
    s_noise = repmat((sn.^2),1,kd);
    sn_delta = s_noise.*(eye(kd));
else
    sn_delta = sn*(eye(kd));
end

c = 0;
for i = 1:dim;
    m = repmat(X(:,i),1,kd);
    c = c + ((m-m.')./L(i)).^2;
end


K = (sf^2)*exp(-c/2) + sn_delta;
display('#### Kernel Calculated ####')
% complete a Cholesky Decomposition of the Kernel Matrix

[cholL,s] = chol(K, 'lower');

% check that the Cholesky Decomposition was successful (the Kernel Matrix
% is positive definite and calculate the inverse if it has been.
if s ~= 0
    cholL = [];
end
    display('#### Cholesky Decomposition Complete ####')
end

