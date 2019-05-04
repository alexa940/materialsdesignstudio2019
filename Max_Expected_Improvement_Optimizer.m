vidname=strcat('Filename'); %Determine file name that will be saved later
% Read in the Original Dataset 
MasterMatrix = readtable('Filename.xlsx');
MasterMatrix = MasterMatrix{:,:};

MasterMatrix = abs(MasterMatrix);   
X=MasterMatrix(:,1:vfresolution);       %X Training dataset from file
Y=-1*MasterMatrix(:,vfresolution+1);    %Y Training dataset from file - Negated for optimizer

maxiter = 100;      %Number of new points to generate
h = waitbar(0,'Iterating'); %Displays how far into the iterations we are

%FEA basis
time=1000;          %Time at which we are optimizing
dT=1;               %Boundary condition
vfresolution=20;    %Degrees of freedom from radial volume fractions

%% Define the Parameters
L = repelem([0.5],vfresolution);   %Lengthscale hyperparameter
sf = 1;         %Noise hyperparameter

ii=1;           %Used to track iterations
while ii < maxiter
    waitbar(ii/maxiter,h);                          %Displays how far into the iterations we are
    Ybest=repmat(min(Y),10000,1);                   %Current highest performance metric output
    [N,VFres] = size(X);                            %N=number of samples VFres=number of volume fractions
    GAUSSIANpdf = @(x)(1/2)*(1+erf(x/sqrt(2)));     %Gaussian probability density function
    GAUSSIANcdf = @(x)1/sqrt(2*pi)*exp(-x.^2/2);    %Gaussian cumulitive distribution function

    Variances=0.1;                                  %Variance
    notestpoints=10000;                             %Number of test points
    X_test = lhsdesign(notestpoints,vfresolution);  %Latin hypercube sampling to generate test space
    X_test = X_test(all( X_test,2),:);              %Delete Rows wtih 0s - check

    %Compute Kernels
    [~, cholL, s] = GP_Kernel(X, L, sf, Variances);
    if s~= 0
        Disp('Could not complete Cholesky Decomposition, Matrix must be positive definite. Consider changing HyperParameters')
    else
        alpha1 = transpose(cholL)\(cholL\Y);
        [ks, kss] = GP_testKernel(X_test, X, L, sf);
        % predict the mean
        Mean = ks*alpha1;
        v = cholL\(ks.');
        % Predict the variance
        Var = kss - (v.')*(v);
    end
    s = sqrt(abs(Var));
    s = s+1e-16;
    % Calculate Expected Improvement:
    test1 = (Ybest-Mean);
    test2 = GAUSSIANcdf(((Ybest-Mean).')/s);
    test3 = GAUSSIANpdf(((Ybest-Mean).')/s);
    EI=-1*(((Ybest-Mean).')*(GAUSSIANcdf(((Ybest-Mean).')/s).')+s*(GAUSSIANpdf(((Ybest-Mean).')/s).'));
     
    %Get the best predicted performance of the test space:
    minimum=min(min(EI));
    [I,y]=find(EI==minimum);
    x_new = X_test(I,:);
    
    %Put new sample point to the list of evaluation points:
    X(end+1,:) = x_new;                                         % Add new point to X training set
    Y(end+1,:) = -1*MDS_functionQonly_t(dT,time,X(end,:));      % Calculate the function value at query point

    ii=ii+1;        %Iteration complete
    disp(Y(end,:)); %Display new added sample performance
    
    %Save every loop:
    MasterMatrix = [X Y];
    Savefile=abs(MasterMatrix);
    save(vidname,'Savefile');
    xlswrite('OptimizingCompKernel.xlsx',Savefile);
    
    disp(N-10000);          %Display size of training set

end
