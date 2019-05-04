function [ Ks, Kss ] = GP_testKernel_2D( testX, X, L, sf )
%GP_testKernel Calculation of the kernel between the test points
%   Calculation of the two vectors (test and Conditioning) to get the
%   additional rows and columns for the prediction of additional points by
%   the GP regressor

[ksd,n] = size(testX);
[kd,~] = size(X);

Ks = zeros(ksd,kd);

for p=1:ksd;
    for q=1:kd;
        c = 0;
        for i = 1:n;
            c = c + ((testX(p,i)-X(q,i))./L(i)).^2;
        end
        Ks(p,q) = (sf^2)*exp( -(1/2)*c);
    end;
end

c = 0;
for i = 1:n;
    m = repmat(testX(:,i),1,ksd);
    c = c + ((m-m.')./L(i)).^2;
end

Kss = (sf^2)*exp(-c/2);

end

