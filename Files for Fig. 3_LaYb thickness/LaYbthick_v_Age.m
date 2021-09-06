%% Analysis of La/Yb based thickness values v time
clear; clc; clf; format shorte;
% Date: July, 2021
% Author: Priydarshi Chowdhury (Monash University)
% >> This code is a part of PNAS manuscript (Chowdhury et al.)
% >> Contact: pchowdhury59@gmail.com / priyadarshi.chowdhury@monash.edu

%% User defined details
%name = input('Enter the file name without extension: ','s');
name = 'Thick LaYb Data'; % name of the data file without extension

% Median plotting (1=yes, 0=no) and bin size if yes
median_plot = 1; binsize = 10;
topX = 3600; botX = 3000; % up. and lr. age limit in Ma for binning

%% Reading the input data
fname  = [name,'.txt']; fid = fopen(fname,'rt');
header = fscanf(fid,'%s',2); data = fscanf(fid,'%f',[2 inf]);
fclose(fid); data = data'; 
%grouping each column into different variables
x_val0 = data(:,1); y_val0 = data(:,2);
fprintf('Total no. of input data = %d\n',length(y_val0));

% Remove NaN from 'y_val0' dataset
A = rmmissing([x_val0,y_val0]);
x_val0= A(:,1); y_val0= A(:,2);
fprintf('Data excluding missing values/NaNs = %d\n',length(x_val0));

%% Remove outliers from 'y_val0' dataset using MEAN method
[y_val1,tf] = rmoutliers(y_val0,'mean');
for i=1:length(y_val0)
    if tf(i)==0
    y_val(i,1) = y_val0(i); x_val(i,1) = x_val0(i);
    else
    y_val(i,1) = NaN; x_val(i,1) = NaN;
    end
end
x_val = x_val(~isnan(x_val)); y_val = y_val(~isnan(y_val));
fprintf('Outliers = %d',length(y_val0)-length(y_val));

% Fitting curve to the LaYb thickness values
myFit = fitlm(x_val,y_val,'quadratic');
display(myFit)

%% Grouping y_val by x_val bins and statistical calc for each bin
% creating the bin edges and bining the x_val
binEdges = botX:binsize:topX;
[N,edges,ind] = histcounts(x_val,binEdges);
for i = 1:length(N)
    flag = (ind==i);
    binElements = y_val(flag); 
    % remove outliers from each bin using quartile method
    binElements = rmoutliers(binElements,'quartiles');
    binMedian (i) = median(binElements);
end
% Mid-point of 'x_val' bins for plotting the median of y_val dataset
for i=1:length(binEdges)-1
    edgeMidpoint(i) = (binEdges(i)+ binEdges(i+1))/2;
end

%% Plotting
figure(1)
set(gcf,'renderer','Painters')

% plotting all thickness values
plot(x_val0,y_val0,'d','MarkerSize',5,'MarkerEdgeColor','k'); hold on;
% plotting median values
plot(edgeMidpoint,binMedian,'or','MarkerFaceColor','r','MarkerSize',10);
hold on;

% plotting the quadratic fit to the data
h = plot(myFit,'marker','none');
h(2).LineWidth = 3;           % set main line thick.
h(2).Color = 'b';             % set main line colour
set(h(3:4),'LineWidth',2.5,'Color','b'); % set thick. of confidence lines
title(['Age vs. ',name,'; No. of Data = ',num2str(length(x_val0)),...
    '; bin-size =',num2str(binsize),' Ma']);
grid off;
xlim([3200 3500]); ylim([25 65]); 
hold off; set(gca,'TickDir','out');
xlabel('Age (Ma)'); ylabel('LaYb thickness'); 
% set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);
legend off;
%%