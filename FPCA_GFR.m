% Use the Matlab functions used in the PACE package

addpath('PACE')

%  Import GFR data

data = cell2mat(slope);
sub = unique(data(:,1));
time = data(:,3);
ncohort = length(sub);
y = cell(1,ncohort);
t = cell(1,ncohort);
lint=10;
xi=zeros(ncohort,2);

for i=1:ncohort
    index = (data(:,1) ==sub(i));
    t{i} = data(index,3);
    t{i}=t{i}';
    y{i} = data(index,2);
    y{i}=y{i}';
end
regular=0;
%========================================================================================================
%New style of setting the input arguments for function FPCA() that calls PCA().
p = setOptions('yname','x', 'regular', regular, 'selection_k', 'FVE','FVE_threshold', 0.99,'screePlot',1, 'designPlot',1,'corrPlot',1,'numBins',0, 'verbose','on');  
%Here, the options "p" is set as the follows: yname = 'x', no. of PC is selected by FVE and FVE_threshold = 0.9,
%create scree plot and design plot (by default, they are not displayed), 
%do not perform prebinning for the input data (by default, the program try to bin the data)
%and the rest uses default values                           
%Note the default value of selection_k is BIC1, which is recommended.
%By default, we display the diagnostic messages. You can suppress the messages by setting 'verbose' to be 'off'.

%Use FPCA() to recover the functional object for y and the results is a cell array
%that can be assessible by getVal(), e.g., to get eigenfunctions "phi", use
%phi = getVal(yy, 'phi');
%This method is recommended!

fprintf(1,'Recover the individual trajectories using FPCA():\n');
time=cputime;
[yy] = FPCA(y,t,p);
time = cputime-time;
display(['Total time in FPCA fitting is : ' num2str(time) ' sec.']);   %extract some of the results for the plots below:
out1 = getVal(yy,'out1');      %vector of time points for mu, phi and ypred
mu = getVal(yy,'mu');          %estimated mean function
out21 = getVal(yy,'out21');    %vector of time points for xcov
xcovfit = getVal(yy,'xcovfit');%fitted covariance surface
xcov = getVal(yy,'xcov');      %estimated smooth covariance evaluated at out21
phi = getVal(yy,'phi');        %estimated eigenfunctions
no_opt = getVal(yy,'no_opt');  %number of optimal FPC scores
xcorr = getVal(yy,'xcorr');    %fitted correlation surface
lambda = getVal(yy,'lambda');
sigma = getVal(yy,'sigma');
sigma_new = getVal(yy,'sigmanew');
error = getVal(p,'error');
method = getVal(p,'method');
shrink = getVal(p,'shrink');
rho = getVal(yy,'rho_opt');
n = length(getVal(yy,'y'));
noeig=getVal(yy,'noeig');

%========================================================================================================

%create Kth mode variation plots for k = 1:no_opt
KModeVariationPlot(yy)  %kth mode variation plots for all k from 1 to no_opt

%KModeVariationPlotJD(yy) change figure  
%Uncomment this if you just want k=1, use
%KModeVariationPlot(yy,1);

%Uncomment this to create scree plot after FPCA()
%createScreePlot(yy);
xi_new=FPCApred(yy,y,t,regular)
xi=getXI(y, t, mu, phi, lambda, sigma, sigma_new, noeig, error, method, shrink, out1, regular, rho);

   
%===========Plot observed, true and predicted curves===================================
ypred=FPCAeval(yy,[],out1);         %obtain predicted curves for all existing subjects
                                    %also allows each subject evaluated at
                                    %different time points
%or simply use ypred = getVal(yy,'ypred'); %all subjects evaluated at the
%same out1

% For prediction of existing or new subjects, use  
% where yy have been obtained from FPCA(). newy and newt are cell arrays
% for new subjects, which have the same structure like y and t.  
[yprednew, xi_new, xi_var] =  FPCApred(yy, newy, newt);
%if all new subjects are evaluated at the same time, newt can be a row vector for 
%time points for one subject

xi(i,:)=[4*randn(1) 1.9*randn(1)];     %generate 2 Principal components
                                       %1st PC score: N(0,9), 2nd PC score: N(0,4)
xtrue=cell(1,ncohort);
for i=1:ncohort
    xtrue{i}=mu_true(out1,lint)+xi(i,:)*xeig(out1,lint,2);
end
%plot 4 randomly selected curves
figure
set(gcf,'Position',[1 29 1247 705]);
k = 1;

for i=mysample(1:ncohort,9,0)        %randomly sample 9 of the curves (without replacement) for plotting
    subplot(3,3,k)
    plot(t{i},y{i},'*',out1,xtrue{i},'b-',out1,ypred{i},'r-');
    title(['Subject ' num2str(i)]);
    xlabel('t');
    ylabel('X(t)');
    k= k+1;
    legend('obs','true','pred','Location','Best')
end
 