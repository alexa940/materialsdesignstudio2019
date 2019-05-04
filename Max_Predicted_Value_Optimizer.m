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

ii=1;               %Used to track iterations
while ii < maxiter
    waitbar(ii/maxiter,h);          %Displays how far into the iterations we are
    Ybest=repmat(min(Y),10000,1);   %Current highest performance metric output

    notestpoints=10000;             %Number of test points
    X_test0 = lhsdesign(notestpoints,vfresolution); %Latin hypercube sampling to generate test space
    X_test(:,1)=X_test0(:,1);       %First value of volume fractions at test points
    randomval=randi(5,1);
    
    %Used to constrain the volume fractions monotonically:
    for column = 2:vfresolution
        X_test(:,column)=X_test(:,column-1)-X_test0(:,column).*X_test(:,column-1)/randomval;
    end
    
    gprMdl = fitrgp(X,Y,'Basis','linear','FitMethod','exact','PredictMethod','exact'); %Fits the Gaussian Proceess Rgressor from the training data
    ypred = predict(gprMdl,X_test); %Predict values in the test space
    
    %Get the best predicted performance of the test space:
    minimum=min(min(ypred));        
    [I,y]=find(ypred==minimum);
    x_new = X_test(I,:);
    
    X(end+1,:) = x_new;% Add new point to X training set
    Y(end+1,:) = -1*MDS_functionQonly_t500(dT,time,X(end,:));     % Calculate the function value at query point and add to Y training set

    disp(Y(end,:));     %Display new added sample performance
    ii=ii+1;            %Iteration complete
    
    %Save every loop
    MasterMatrix = [X Y];
    Savefile=abs(MasterMatrix);
    save(vidname,'Savefile');
    xlswrite('OptimizingBI500contin.xlsx',Savefile);
    
    disp(size(X,1));    %Display size of training set

end
