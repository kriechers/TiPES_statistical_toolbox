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
% Author: Witold Bagniewski Date: 31/07/2022
%
%% Load data
Fname = 'ngrip_d18o_20y';  ext = '.txt'; 
data=load([Fname ext]);

% % remove equal dates, delete NaNs
[C,~,idx] = unique(data(:,1),'stable');
val = accumarray(idx,data(:,2),[],@mean);
data1 = [C val];  clearvars C ia idx val 
data1(isnan(data1(:,1)),:) = [];
data1(isnan(data1(:,2)),:) = [];
t = data1(:,1);
x = data1(:,2);

%% Parameters for RQA and for plotting

epsl = 1.3; % epsilon (recurrence threshold)
wind = 1;	% window size for RR
y_label = '\delta^{18}O [‰]';
x_label = 'Age (ky b2k)'; 
updown = 0;            % don't plot the timeseries upside down
% updown = 1;            % plot the timeseries upside down

%% Data options

t = t/1000;                % convert yr to kyr

% % use uniform sampling rate (in case it is nonuniform)
% s_rate=0.02; 
% i=1;
% for y = data1(1,1):s_rate:data1(end,1)
%     r1=find(data1(:,1) > y-s_rate/2 & data1(:,1) <= y+s_rate/2);
%     d2(i,1)=y;
%     d2(i,2)=mean(data1(r1,2));
%     i=i+1;
% end
% d2(isnan(d2(:,1)),:) = [];	   % delete NaN rows
% d2(isnan(d2(:,2)),:) = [];	   % delete NaN rows
% xq = data1(1,1):s_window:data1(end,1);
% vq1 = interp1(d2(:,1),d2(:,2),xq);
% clear data;
% t(:,1)=xq;
% x(:,2)=vq1;

%% Calculate the recurrence matrix

N = length(x);
S = zeros(N,N);
for i = 1:N
    S(:,i) = abs(repmat(x(i),N,1) - x(:));
end

%% Calculate the Recurrence Rate (RR)

clear RR
for i=1:length(t)
    if t(i)+wind <= t(end)
        tt=find(t(i)<=t & t<=t(i)+wind);
        RR(i,1)=mean(t(tt));
        RR(i,2)=length(find(S(tt,tt) < epsl))/(length(tt)*length(tt)) ;
    end
end

% Find most prominent RR peaks (threshold = standard deviation of RR)
[~,pks,~,p]=findpeaks(-RR(:,2),'MinPeakProminence',std(RR(:,2)));

%% Plotting
% The plotting method is based on Matlab code used in the paper: 
% Trauth, M.H., Asrat, A., Duesing, W., Foerster, V., Kraemer, K.H.,
% Marwan, N., Maslin, M.A., Schaebitz, F. (2019) Classifying past climate
% change in the Chew Bahir basin, southern Ethiopia, using recurrence
% quantification analysis. Climate Dynamics, Springer Verlag GmbH Germany,
% https://doi.org/10.1007/s00382-019-04641-3

figure('Position',[400 600 400 600],'Color',[1 1 1])

% Time series.
h(1) = axes('Position',[0.05 0.83 0.9 0.16],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(x)-0.1 max(x)+0.1],...
     'LineWidth',0.75,...
     'XGrid','On',...
     'XMinorTick','on',...
     'Box','On',...
     'XAxisLocation','top'); hold on
line(t,x,'LineWidth',0.5,'Color',[0 0 0])
ylabel(y_label);
if updown==1;axis ij;end

% Recurrence plot.
h(2) = axes('Position',[0.05 0.06 0.9 0.9],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(t) max(t)],...
     'LineWidth',0.75,...
     'Box','On',...
     'XTick',[],...
     'YTick',[]); hold on
axis square xy
imagesc(t, t, S < epsl)
axes('Position',[0.05 0.06 0.9 0.9],...
     'XLim',[min(t) max(t)],...
     'YLim',[min(t) max(t)],...
     'LineWidth',0.75,...
     'Box','On',...
     'Color','none',...
     'XTick',[],...
     'YTick',[]), hold on
axis square xy
colormap([1 1 1; 0 0 0])

% RQA - Recurrence rate.
h(3) = axes('Position',[0.05 0.03 0.9 0.16],...
    'XLim',[min(t) max(t)],...
    'YLim',[0 1],...
    'LineWidth',0.75,...
    'XGrid','On',...
    'XMinorTick','on',...
    'Box','On'); hold on
line(h(3),RR(:,1),RR(:,2),'LineWidth',1,'Color',[0 0.5 0.75])
hold on; scatter(RR(pks,1),RR(pks,2),50,'x','m','LineWidth',1)
xlabel(x_label); ylabel('Recurrence rate');

% Print.
printname = [Fname '_RR.pdf']; 
print(printname,'-dpdf', '-fillpage');