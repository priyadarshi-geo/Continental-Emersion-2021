%% Analysing the compositional evolution of Singhbhum TTGs
clear; clc; clf; format shorte;
% Date: July, 2021
% Author: Priydarshi Chowdhury (Monash University)
% >> This code is a part of PNAS manuscript (Chowdhury et al.)
% >> Contact: pchowdhury59@gmail.com / priyadarshi.chowdhury@monash.edu

%% User defined details
name = 'SiO2 Data'; % name of the data file without extension
% Median plotting (1=yes, 0=no) and bin-size if yes
median_plot = 1; 
binsize = 10;
% bin size is 7 Ma for Sr/Y & La/Yb data; and 10 Ma for major element
% proxies (mentioned in Methods)

% up. and lr. age limit in Ma for binning
upX = 3605; lowX = 2995;

%% READING THE INPUT X-Y DATA
fname  = [name,'.txt']; fid = fopen(fname,'rt');
header = fscanf(fid,'%s',2); data = fscanf(fid,'%f',[2 inf]);
fclose(fid); data = data'; 
% assigning variables to the input data; these variables will be kept
% throughout the simulation
x_val0 = data(:,1); y_val0 = data(:,2);
fprintf('Total no. of input data = %d\n',length(y_val0));

% Remove all those x-y pairs either one of which is NaN or not defined
A = rmmissing([x_val0,y_val0]);
x_val0= A(:,1); y_val0= A(:,2);
fprintf('Data excluding missing values/NaNs = %d\n',length(x_val0));

%% CALCULATION PERFORMED ON THE WHOLE y-DATASET (NO GROUPPING OF DATA)

% STEP 1: REMOVE OUTLIERS FROM 'y_val0' DATA
[~,tf] = rmoutliers(y_val0,'median');
% Need to remove the corresponding 'x_val0' data
for i=1:length(y_val0)
    if tf(i)==0
        y_val(i,1) = y_val0(i); x_val(i,1) = x_val0(i);
    else
        y_val(i,1) = NaN; x_val(i,1) = NaN;
    end
end
clear A; A = rmmissing([x_val,y_val]);
x_val= A(:,1); y_val= A(:,2);
fprintf('Outliers = %d',length(y_val0)-length(y_val));

% STEP 2: KERNEL DENSITY DISTRIBUTION (KDD) CALCULATION
% Gridding of 'x_val-y_val' space
gridx = linspace(0.9*min(x_val),1.1*max(x_val),200);
gridy = linspace(0.0*min(y_val),1.5*max(y_val),200);
[x,y] = meshgrid(gridx, gridy);
% getting the x- & y-coord of all the nodes
XX = x(:); YY = y(:); node = [XX YY];
KDen_vec(:,i) = ksdensity([x_val,y_val],node);
KDen_grid     = reshape(KDen_vec(:,i),length(gridy),length(gridx));

% STEP 3: FITTING THE 'x_val-y_val' DATA & DISPLAY FITTING STATS
myFit = fitlm(x_val,y_val,'linear');
display(myFit)

%% CALCULATION PERFORMED AFTER BINNING THE y-DATASET

% STEP 1: CREATING BINS & BIN-EDGES
% Creating the bin edges according to 'binsize' for 'x_val'
binEdges = lowX:binsize:upX;
% Mid-point of 'x_val' bins for plotting
for i=1:length(binEdges)-1
    binEdges_Mid(i) = (binEdges(i)+ binEdges(i+1))/2;
end

% STEP 2: BINNING THE X-DATA
[N,edges,ind] = histcounts(x_val,binEdges);

% STEP 3: BINNING THE Y-DATA ACCORDING TO X-BINS
group = NaN(length(N),length(x_val0));
for i = 1:length(N)
    flag = (ind==i);
    binElements = y_val(flag); 
    % Storing bin-elements for each bin as a row of 'group' matrix
    if isempty(binElements)== 0
        for k=1:length(binElements)
            group(i,k)=binElements(k);
        end
    end
end

% performing the statistics for each row of 'group'
for k=1:length(N)
    if isempty(group(k,:))== 0
        % remove outliers from each bin using quartile method
        group_data     = group(k,:);
        group_data     = rmoutliers(group(k,:),'quartiles');
        binMedian (k)  = nanmedian(group_data);
    end
end
% STEP 4: FITTING THE MEDIAN DATA & DISPLAY THE FITTING STATS
myFit2 = fitlm(binEdges_Mid,binMedian,'linear');
display(myFit2)

%% PLOTTING FOR CALCULATION-1 (WITHOUT GROUPPING OF DATA)
figure(1)
set(gcf,'renderer','Painters')
% KDD plot as contoured colored field with the data overlayed
colormap(jet(20))
contourf(gridx,gridy,KDen_grid,19,'-',...
    'LineColor','w','LineStyle','-','ShowText','off');
hold on;
shading interp; grid on;

% Plotting all data
plot(x_val0,y_val0,'d','MarkerSize',10,'MarkerEdgeColor','r',...
    'MarkerFaceColor','w'); hold on;

% Plotting the curve fitting to the original dataset (myFit)
h = plot(myFit,'marker','none');
h(2).LineWidth = 3; h(2).Color = 'r'; % set main line thick & color
set(h(3:4),'LineWidth',2.5,'Color','r'); % set thick. of both confid. lines

xlim([3000 3600]); if strcmp(name, 'SiO2 Data')== 1, ylim([60 80]); end
title([name,' vs. Age ',';  No. of data = ',num2str(length(x_val0))]);
hold off; set(gca,'TickDir','out');
xlabel('X-DATA'); ylabel('Y-DATA'); % set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);
grid off; legend off


%% PLOTTING FOR CALCULATION-2 (GROUPPING OF DATA)
figure(2)
set(gcf,'renderer','Painters')

% KDD plot of the grouped data median
KDen_vec2(:,i) = ksdensity([binEdges_Mid',binMedian'],node);
KDen_grid2     = reshape(KDen_vec2(:,i),length(gridy),length(gridx));

% colormap for KDD
colormap(jet(20))
contourf(gridx,gridy,KDen_grid2,19,'-',...
    'LineColor','w','LineStyle','-','ShowText','off');
hold on;
shading interp; grid on;

% Plotting all data
%plot(x_val0,y_val0,'d','MarkerSize',8,'MarkerEdgeColor','w'); hold on;
% Plotting the median for each bin
plot(binEdges_Mid,binMedian,'o','MarkerSize',14,'MarkerEdgeColor','r',...
    'MarkerFaceColor','w');
hold on;

% Plotting the curve fitting to the median values (myFit2)
hh = plot(myFit2,'marker','none');
hh(2).LineWidth = 3; hh(2).Color = 'r'; % set main line thick & color
set(hh(3:4),'LineWidth',2.5,'Color','r'); % set thick. of both confid. lines
hold on;

title([name,' vs. Age',';  No. of data = ',num2str(length(x_val0)),...
    ';  Median binsize = ',num2str(binsize), ' Ma']);
xlim([3000 3600]); if strcmp(name, 'SiO2 Data')== 1, ylim([60 80]); end
hold off; set(gca,'TickDir','out');
xlabel('Age (Ma)'); ylabel('Y-DATA'); % set(gca,'xtick',[],'ytick',[]);
set(gca,'layer','top','linewidth',0.5);
grid off; legend off

% outfname= [name,'_res','.xlsx'];
% writematrix([binEdges_Mid',binMedian'],outfname,'Range','B2');
%%