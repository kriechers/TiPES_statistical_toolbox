%% DESCRIPTION
% MATLAB script to perform Recurrence Quantification Analysis (RQA) on the
% NGRIP d18O record (Rasmussen et al., 2014). This script calculates the
% recurrence matrix and the recurrence rate, and finds the most prominent
% RR peaks. The results are plotted and saved as a pdf file. This script
% creates Figure 8 of the paper:
% Bagniewski, W., Ghil, M., Rousseau, D.-D., (2021) Automatic detection of
% abrupt transitions in paleoclimate records. Chaos, 31(11):113129,
% https://doi.org/10.1063/5.0062543
%
% Author: Witold Bagniewski 
% Date: 18/01/2023
%
%% Load data

Fname = 'NGRIP';  ext = '.csv';
data = readcell([Fname ext], 'CommentStyle', '#');
data = cell2table(data(2:end,:), 'VariableNames', data(1,:));
Vname = 'd18O';  % Variable name used for saving the figure
t_column = 1;    % Column containing time
x_column = 2;    % Column containing variable
t = table2array(data(:,t_column));   % this assumes the first row contains variable names
x = table2array(data(:,x_column));

% % Load .txt file
% Fname = 'ngrip_d18o_20y';  ext = '.txt'; 
% data = load([Fname ext]);
% t = data(:,1);
% x = data(:,2);

%% Preprocessing the time series

% Delete NaN rows
idx = ~isnan(t+x);
t = t(idx);
x = x(idx);

% For dates that are equal find the mean value of x
[tt,~,idx] = unique(t,'stable');
xx = accumarray(idx,x,[],@mean);

%% Parameters for RQA and for plotting

epsl = 1.3; % epsilon (recurrence threshold)
w_min = 1;  % Size of the smallest window
w_max = 1;  % Size of the largest window
n_w = 1;    % Number of distinct window lengths

% updown = 1;    % plot the timeseries upside down
x_label = char(data.Properties.VariableNames(t_column)); % Set time label
y_label = char(data.Properties.VariableNames(x_column)); % Set variable label

% % Plot part of the time series
% t_lim = [0,250]; 
% tt = t(t_lim(1)<=t & t<=t_lim(2)); 
% xx = x(t_lim(1)<=t & t<=t_lim(2));
% t = tt; x = xx;

%% Use uniform sampling rate (in case it is nonuniform)

for i = 1:length(t)-1
    t_incr(i) = t(i+1)-t(i);     % calculate the time increments
end

if max(t_incr)-min(t_incr) > abs(mean(t_incr))*0.001   % test if uniform or not
    s_rate = 0.02;                 % set the new sampling rate
    i=1;
    for y = t(1):s_rate:t(end)
        r1=find(t > y-s_rate/2 & t <= y+s_rate/2);
        d2(i,1)=y;
        d2(i,2)=mean(x(r1));
        i=i+1;
    end
    d2(isnan(d2(:,1)),:) = [];	   % delete NaN rows
    d2(isnan(d2(:,2)),:) = [];	   % delete NaN rows
    xq = t(1):s_rate:t(end);
    vq1 = interp1(d2(:,1),d2(:,2),xq);
    clear t x;
    t=xq;
    x=vq1;
end

%% Calculate the recurrence matrix

N = length(x);
S = zeros(N,N);
for i = 1:N
    S(:,i) = abs(repmat(x(i),N,1) - x(:));
end

%% Calculate the Recurrence Rate (RR)

clear RR RRs p pks
if or(w_min == w_max , n_w==1)   % in case there is only one n_w window size
    wind = w_min;
    n_w = 1;
else
    for i=1:n_w
        wind(i)=w_min*(w_max/w_min)^((i-1)/(n_w-1));   % makes n_w window sizes
    end
end

for j=1:length(wind)
    for i=1:length(t)
        if or(t(i)-wind(j)/2 < min(t) , t(i)+wind(j)/2 > max(t))
            RRs(i,j) = NaN;
        else
            tt=find(t>=t(i)-wind(j)/2 & t<=t(i)+wind(j)/2);
            RRs(i,j)=length(find(S(tt,tt) < epsl))/(length(tt)*length(tt)) ;
        end
    end
end
if length(RRs(1,:)) > 1
    RR=nanmean(RRs')';
else
    RR=RRs;
end

% Find most prominent RR peaks (threshold = standard deviation of RR)
nn = isnan(RR);
RR(nn) = max(RR);   % this ensures that the first and last peaks can be found
[~,pks,~,p]=findpeaks(-RR,'MinPeakProminence',std(RR));
RR(nn) = NaN;

%% Plotting
% The plotting method is based on Matlab code used in the paper: 
% Trauth, M.H., Asrat, A., Duesing, W., Foerster, V., Kraemer, K.H.,
% Marwan, N., Maslin, M.A., Schaebitz, F. (2019) Classifying past climate
% change in the Chew Bahir basin, southern Ethiopia, using recurrence
% quantification analysis. Climate Dynamics, Springer Verlag GmbH Germany,
% https://doi.org/10.1007/s00382-019-04641-3

figure('Position',[400 600 400 600],'Color',[1 1 1])

% Time series.
h(1) = axes('Position',[0.055 0.819 0.9 0.16],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(x)-(max(x)-min(x))/50 max(x)+(max(x)-min(x))/50],...
     'LineWidth',0.75,...
     'XGrid','On',...
     'XMinorTick','on',...
     'Box','On',...
     'XAxisLocation','top'); hold on
line(t,x,'LineWidth',0.5,'Color',[0 0 0])
h(1).XRuler.TickLabelGapOffset = -1;
ylabel(y_label);
title(Fname,'FontWeight','Normal','Interpreter','none');
set(gca,'FontSize',10)
if exist('updown') && updown==1; axis ij; end

% Recurrence plot.
h(2) = axes('Position',[0.055 0.049 0.9 0.9],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(t) max(t)],...
     'LineWidth',0.75,...
     'Box','On',...
     'XTick',[],...
     'YTick',[]); hold on
axis square xy
imagesc(t, t, S < epsl)
axes('Position',[0.055 0.049 0.9 0.9],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(t) max(t)],...
     'LineWidth',0.75,...
     'Box','On',...
     'Color','none',...
     'XTick',[],...
     'TickDir','out',...
     'YMinorTick','on'), hold on
ylabel(x_label);
axis square xy
colormap([1 1 1; 0 0 0])
set(gca,'FontSize',10)

% RQA - Recurrence rate.
h(3) = axes('Position',[0.055 0.019 0.9 0.16],...
    'XLim',[min(t) max(t)],...
    'YLim',[0 1],...
    'LineWidth',0.75,...
    'XGrid','On',...
    'XMinorTick','on',...
    'Box','On'); hold on
line(h(3),t,RR,'LineWidth',1,'Color',[0 0.5 0.75])
h(3).XRuler.TickLabelGapOffset = -1;
hold on; scatter(t(pks),RR(pks),50,'x','m','LineWidth',1)
xlabel(x_label); ylabel('Recurrence rate');
set(gca,'FontSize',10)

% Print.
printname = [Fname '_' Vname '_RQA.pdf']; 
print(printname,'-dpdf', '-fillpage');

%% Save RQAjumps dates

clear par_n
if exist('w_min','var'); par_n{1} = ['w_min = ' num2str(w_min)]; end
if exist('w_max','var'); par_n{end+1} = [', w_max = ' num2str(w_max)]; end
if exist('n_w','var'); par_n{end+1} = [', n_w = ' num2str(n_w)]; end
if exist('epsl','var'); par_n{end+1} = [', epsilon = ' num2str(epsl)]; end

fileID = fopen([Fname '_' Vname '_RQAjumps.csv'],'w');
fprintf(fileID,'%s\n\n',['# RQA parameters: ' strcat(par_n{:})]);
fprintf(fileID,'time,RR prominence\n');
fprintf(fileID,'%g,%g\n',[t(pks),p]');
fclose(fileID);