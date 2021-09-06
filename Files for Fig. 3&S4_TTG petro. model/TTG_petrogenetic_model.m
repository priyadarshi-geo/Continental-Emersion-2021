%% Trace element partitioning between TTG-residua during melting
clear; clc; clf; format shorte;
% Date: July, 2021
% Author: Priydarshi Chowdhury (Monash University)
% >> This code is a part of PNAS manuscript (Chowdhury et al.)
% >> Contact: pchowdhury59@gmail.com / priyadarshi.chowdhury@monash.edu

%% Reading input files- Kd values and modebox data
fname = 'EAT750'; % four options: EAT750; EAT1000; DAT750; DAT1000
% these file names refer to EAT or DAT source compositions and geotherm of
% melting (750 or 1000 C/GPa); phase proportions taken from Kendrick and
% Yakymchuk, 2019, Geoscience Frontiers]

[data,row,column,phase_MB,Null] = read_datafile([fname,'.dat']);
phase_prop  = data(:,4:end); 
press = data(:,1); temp = data(:,2); meltF = data(:,3);

fname2 = 'Kd.dat'; % file containing Kd values
[Kd,rowKd,TEnum,elemKd,phaseKd] = read_datafile(fname2);

% source rock composition (See Methods for details and data source)
EAT = [10.00,12.80,0.78,13.00,30.00,17.00,450.00,4.00,8933.00,1.30,...
            3.80,0.70,4.20,30.00,2.30,2.20,0.38];
DAT = [4.0,5.0,0.29,3.6,9.2,6.6,100.0,2.0,5635.0,0.73,2.6,0.5,3.1,20.0,...
            2.0,1.9,0.31];
source = EAT; 
% these should be changed between EAT or DAT compositions and in accordance
% to the phase proportions used in Line 5. Like, source should be EAT when
% EAT750 or EAT1000 is used in Line 5. VERY IMPORTANT!!!
        
% File with data from the craton for comparison with modelled results
craton =  1; fname3 = 'Age-G2 TTG';
% the TTG data is given altogether in 'All TTG. dat' file, and also in
% 3 age-groups as shown in Fig. 3 and Fig. S4, namely 'Age-G1 TTG.dat', 
% 'Age-G2 TTG.dat' and 'Age-G3 TTG.dat' files
        
% chondrite norm. values
PM = [0.635,0.713,0.041,0.687,1.775,1.354,21.100,0.444,1300,0.168,...
            0.596,0.108,0.737,4.550,0.480,0.493,0.074];

%% Calculation begins...
for i=1:length(meltF)
    % renorm. the residual phases prop. to 100 for bulkD calc.
    sum_phase_prop = sum(phase_prop(i,:));
    Res(i,:) = (100.*phase_prop(i,:))./sum_phase_prop;
    % converting phase prop. to fraction for bulkD calc.
    Res(i,:) = Res(i,:)./100; F = meltF(i)/100;
    
    for j=1:TEnum
        bulkD(i,j) = Res(i,:)*Kd(:,j);
        meltTE(i,j) = source(j)/(F+bulkD(i,j)*(1-F));
    end
end

%% writing results to an excel file
outfname= [fname3,'-',fname,'.xlsx'];
elemKd(1,1) = {'melt(%)'};
res_heading= [{'P (kbar)'},elemKd];

% writing the melt-compositions
writecell(res_heading,outfname,'Sheet','Result','Range','B2');
writematrix([press,meltF,meltTE],outfname,...
    'Sheet','Result','Range','B3');

% writing the phase propotion used for each melt-fraction
writecell(phase_MB(1,3:end),outfname,'Sheet','Phase Prop.','Range','B2');
writematrix([meltF,Res],outfname,'Sheet','Phase Prop.','Range','B3');

% writing the bulk-D calcualted using Kd and phase proportions
writecell(elemKd,outfname,'Sheet','Bulk D','Range','B2');
writematrix([meltF,bulkD],outfname,'Sheet','Bulk D','Range','B3');

%% Plotting the modelled trace element data
clf; figure(1)
if craton == 1
    [measTE,data_num,TEnum2,elem,sample] = read_datafile([fname3,'.dat']);
    if isempty(setdiff(elem(1,2:end),elemKd(1,2:end)))==1 &&...
            TEnum == TEnum2
        for m=1:data_num
            measTE_N(m,:)= measTE(m,:)./PM;
        end
        patch([1:TEnum2 fliplr(1:TEnum2)],[min(measTE_N) ...
            fliplr(max(measTE_N))],[.95 .95 .95],"LineStyle","None",...
            'DisplayName','Range'); hold on;
        plot(nanmean(measTE_N),"LineWidth",3,"Color",'k',...
            'DisplayName','Avg.');hold on;
        boxplot(measTE_N,'OutlierSize',.5); hold on;
    end
end
set(gca,'YScale','log'); grid on;
for k=1:length(meltF)
    plot(meltTE(k,:)./PM,'LineWidth',1,'DisplayName',num2str(meltF(k)));
    hold on;
end
ylim([3e-1,3e2]); legend show; set(gca,'YScale','log'); 
grid on; hold off; set(gca,'TickDir', 'both');
set(gca, 'XTick', 1:TEnum, 'XTickLabel', elemKd(2:end));

%% Main code ends here

%% To read different form of input data
function [output,row,column,HEADINGS,LABELS] = read_datafile(f_name)
% Reads matrix form of data from a file & also gives its dimension
%
% FORMAT : read_datafile(f_name)
%
%    valid form of data: 
%       a. Matrix/Vector type with/without strings in first row and/or- column(HEADER 
%          STRINGS) with 'END OF LINE' character '\n' 
%       b. single array of numbers  without any string & 'END OF LINE' character '\n' 
%       c. lines with 'LINE SKIP' character '%'
%       d. matrix of data (with/without Header strings) with all column &/ rows with
%          different lengths  
format short

% opening & reading data from file
fid          = fopen(f_name,'rt');
if fid == -1
    fprintf('>> ERROR: opening of the file failed  \n');
else
    data_read    = textscan(fid,'%s','CommentStyle','%');
end
fclose(fid);

% determining the size of cell array
[~,c_data] = size(data_read);

% reading is valid iff 1-by-1 cell array  (it means entire data is read in the form of a column)
% Since a cell-array has always 1 row,therefore checking only the column number of cell-array

if c_data == 1
    
    % Transforming the read data into a matrix/vectror(array)form
    % Initializing counters for row and column indices
    m        = 1;                               % row index
    n        = 1;                               % column index
    
    % Starting the loop for checking the "END OF LINE" character and forming the matrix
    for k = 1:length(data_read{1})              % checking over the entire 1-by-1 cell-array
        
        % comparing each element of the cell, whether it is the 'END OF LINE' character or not..!!
        com = strcmp(data_read{1}(k),'\n');     
        
        if com == 0         % comparision is false, not similar to 'END OF LINE' character
            
            % storing the element in the matrix, & the position is defined by row-index = m and column-index = n
            data(m,n) = data_read{1}(k);
            % increasing 'n', because next element will be stored in next column (concatenating row-wise)
            n = n+1;
                
        else                % comparision is true,similar to 'END OF LINE' character
            
            % no need of storing this element of the cell, rather increaing
            % row index- 'm' by 1 and setting column index- 'n' = 1, as further storing will start from next
            % row but 1st column
            m = m+1;
            n = 1;
        end
    end
    
    % Storing the first row and first column of input data incase they are strings
    HEADINGS  = data(1,:);
    LABELS    = data(2:end,1);
end

% disp(data);

% Converting the data-matrix from string - cell array to double-type matrix
data     = str2double(data);

% Determining the size of the modified matrix
[~,c]= size(data);

% Checking whether first row & firat column of the matrix contains strings or numeric value
% since 'str2double' returns NaN for strings, therefore comparing element with NaN

% Checking for the first row
for l = 1:1
    % initializing a counter
    counter = 0;
    
    % starting the loop over entire row
    for k = 1:c
        com = isnan(data(l,k));
        % com=1 means the element is NaN
        if double(com) == 1
            % keeping a count of non-numeric elements
            counter = counter+1;
        end
    end
    
    % deleting the entire row iff all the elements of the row are non-numeric
    % therefore, checking the counter with no. of column
    if counter == c
        % deleting the entire row of the matrix
        data(l,:)=[];
    end
end

% Determining the size of the modified matrix
[r,~]= size(data);

% Checking for the first column
for k = 1:1
    % initializing a counter
    counter = 0;
    
    % starting the loop over entire column
    for l = 1:r
        com = isnan(data(l,k));
        % com=1 means the element is NaN
        if double(com) == 1
            % keeping a count of non-numeric elements
            counter = counter+1;
        end
    end
    
    % deleting the entire column iff all the elements of the column are non-numeric
    % therefore, checking the counter with no. of row
    if counter == r
        % deleting the entire column of the matrix
        data(:,k)=[];
    end
end

% determining the dimension of the numeric-data matrix
[row,column] = size(data);
output = data;
end
% End of function

