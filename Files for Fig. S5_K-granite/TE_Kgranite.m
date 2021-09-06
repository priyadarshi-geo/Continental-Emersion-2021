%% Trace element plot of Singhbhum K-granites
clear; clc; clf; format shorte;
% Date: July, 2021
% Author: Priydarshi Chowdhury (Monash University)
% >> This code is a part of PNAS manuscript (Chowdhury et al.)
% >> Contact: pchowdhury59@gmail.com / priyadarshi.chowdhury@monash.edu

%% Reading input files- Kd values and modebox data
fname3 = 'K-granites';
        
% Norm. values
PM = [0.635,0.713,0.041,0.687,1.775,1.354,21.100,0.444,1300,0.168,...
            0.596,0.108,0.737,4.550,0.480,0.493,0.074];

% Plotting the modelled trace element data
[measTE,data_num,TEnum2,elem,sample] = read_datafile([fname3,'.dat']);
for m=1:data_num
    measTE_N(m,:)= measTE(m,:)./PM;
end

% plotting the field
patch([1:TEnum2 fliplr(1:TEnum2)],[min(measTE_N) ...
    fliplr(max(measTE_N))],[.95 .95 .95],"LineStyle","None",...
    'DisplayName','Range'); hold on;
% plotting mean of the dataset per element
plot(nanmean(measTE_N),"LineWidth",3,"Color",'k',...
    'DisplayName','Avg.');hold on;
boxplot(measTE_N,'OutlierSize',.5); hold on;

set(gca,'YScale','log'); grid on;
ylim([3e-1,3e2]); legend show; set(gca,'YScale','log');
grid on; hold off; set(gca,'TickDir', 'both');
set(gca, 'XTick', 1:TEnum2, 'XTickLabel', elem(2:end));

% writing the normalized values into an excel file
% outfname= [fname3,'_res','.xlsx'];
% elem(1,1) = {'Sl. No.'};
% writecell(elem(2:end),outfname,'Sheet','Meas. TE-PM','Range','B2');
% writematrix(measTE_N,outfname,...
%     'Sheet','Meas. TE-PM','Range','B4');

%%
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

