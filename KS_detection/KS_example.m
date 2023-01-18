%% DESCRIPTION
% MATLAB script to detect abrupt transitions in the NGRIP d18O record
% (Rasmussen et al., 2014). This script identifies the transitions as well
% as the glacial/interglacial boundaries. The results are plotted and saved
% as a pdf file. Dates of the transitions are saved to a csv file.
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

%% Preprocessing the time series (optional)

% % Delete NaN rows
idx = ~isnan(t+x);
t = t(idx);
x = x(idx);

% % For dates that are equal find the mean value of x
[t,~,idx] = unique(t,'stable');
x = accumarray(idx,x,[],@mean);

% % Data options
t1=t; x1=x;                 % these will be used to plot the unaltered data
% x = log10(x); logdata = 1;  % logarithmic scale
% x = -x; updown = 1;         % flip upside down

%% Plotting options

plot_gigs = 1;                    % show stadial/interstadial boundaries
xtck = 0:2:5000;                  % x ticks
xtck_m = 0:0.2:5000;              % small x ticks
t_lim = [min(t1),max(t1)];        % Plot the entire time series
%t_lim = [7.8,122.3];              % Plot only part of the time series
x_label = char(data.Properties.VariableNames(t_column)); % Set time label
y_label = char(data.Properties.VariableNames(x_column)); % Set variable label
% x_label = 'Age [kyr b2k]';       % Use a custom time label
% y_label = '\delta^{18}O [‰]';   % Use a custom variable label

%% Run the augmented KS test

% % KS_detection parameters
w_min = 0.12;
w_max = 2.5;
n_w   = 15;
d_c   = 0.77;
n_c   = 3;
s_c   = 2;
x_c   = 0.8;

% % Run KS_detection
jumps = KS_detection(t,x,w_min,w_max,n_w,d_c,n_c,s_c,x_c);

%% Find "stadial - interstadial" boundaries (optional, may take a long time)

if exist('plot_gigs','var') && plot_gigs == 1 && any(jumps(:,2)==1) && any(jumps(:,2)==-1)

clear mmmGIGS gi_gs G_I G_S w_gigs wgigs_min wgigs_max n_wgigs jumps_gigs;
wgigs_min = 5;   % min window size
wgigs_max = 60;  % max window size
n_wgigs = 10;    % number of windows

for i=1:n_wgigs
    w_gigs(i)=wgigs_min*(wgigs_max/wgigs_min)^((i-1)/(n_wgigs-1));   % makes n_w window sizes
end
kk=1;
for i=w_gigs
    [G_I(:,kk),G_S(:,kk)] = findgigs(t,x,i);
    mmmGIGS(:,kk)=(G_I(:,kk)+G_S(:,kk))/2;
    kk=kk+1;
end
mmmGIGS=mean(mmmGIGS');

gi_gs=[];
jumps_gigs=jumps;

if mean(x(t>jumps_gigs(end,1)))>mmmGIGS(end)	% Check if it starts with glacial or interglacial
    mark=1; gi_gs=t(end); k=2; jumps_gigs(end+1,:)=[t(end),1];
else
    mark=0; k=1; jumps_gigs(end+1,:)=[t(end),0];
end

for i = length(jumps_gigs):-1:1
    if i>1
        r2=find(t >= jumps_gigs(i-1,1) & t < jumps_gigs(i,1)); 
    else
        r2=find(t >= t(1) & t < jumps_gigs(i,1)); 
    end
    m2 = mean(x(r2));
    if m2 > mmmGIGS(max(r2)+1) && jumps_gigs(i,2) == 1 && mark == 0 
        gi_gs(k) = jumps_gigs(i,1); mark = 1; k=k+1;
    elseif m2 < mmmGIGS(max(r2)+1) && jumps_gigs(i,2) == -1 && mark == 1
        gi_gs(k) = jumps_gigs(i,1); mark = 0; k=k+1;
    end
end
if jumps_gigs(jumps_gigs(:,1)==gi_gs(end),2) == 1
    gi_gs(end+1)=t(1);
end
gi_gs=flip(gi_gs);

else
    gi_gs = [];
end

%% Plotting

jump_u = jumps(jumps(:,2)==1);
jump_d = jumps(jumps(:,2)==-1);

time3 = (t_lim(2)-t_lim(1))/3;
xl1=[t_lim(1),t_lim(1)+time3*1.026];
xl2=[t_lim(1)+time3*0.987,t_lim(1)+time3*2.013];
xl3=[t_lim(1)+time3*1.974,t_lim(1)+time3*3];
xtlim = x1(t_lim(1)<=t1 & t1<=t_lim(2)); 
yl = [min(xtlim), max(xtlim)]; 

figure; subaxis(3,1,1,1,'Spacing',0,'SpacingHoriz',0.3,'Margin',0.002,'MarginBottom',0.004,'MarginTop',0.003,'MarginLeft',0.075,'Paddingbottom',0.037);
if exist('logdata','var') && logdata == 1
    semilogy(0,0)
end
if exist('plot_gigs','var') && exist('gi_gs','var') && plot_gigs == 1
    for i=1:2:length(gi_gs)-1
        rectangle('Position',[gi_gs(i),min(x1),gi_gs(i+1)-gi_gs(i),max(x1)-min(x1)],'FaceColor', [0.73,0.73,0.73],'LineStyle','none'); hold on
    end
end
if exist('logdata','var') && logdata == 1
    semilogy(t1,x1,'k')
else
    plot(t1,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl1); ylim(yl);
if exist('updown','var') && updown==1; axis ij; end
xticks(xtck); ylabel(y_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
if exist('logdata','var') && logdata == 1 && length(yticks)==1; hA.YAxis.TickValues = [hA.YAxis.TickValues/5 hA.YAxis.TickValues hA.YAxis.TickValues*5]; end
title(Fname,'FontWeight','Normal','Interpreter','none');
set(gca,'FontSize',11)

subaxis(3,1,1,2)
if exist('logdata','var') && logdata == 1
    semilogy(0,0)
end
if exist('plot_gigs','var') && exist('gi_gs','var') && plot_gigs == 1
    for i=1:2:length(gi_gs)-1
        rectangle('Position',[gi_gs(i),min(x1),gi_gs(i+1)-gi_gs(i),max(x1)-min(x1)],'FaceColor', [0.73,0.73,0.73],'LineStyle','none'); hold on
    end
end
if exist('logdata','var') && logdata == 1
    semilogy(t1,x1,'k')
else
    plot(t1,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl2); ylim(yl); 
if exist('updown','var') && updown==1; axis ij; end
xticks(xtck); ylabel(y_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
if exist('logdata','var') && logdata == 1 && length(yticks)==1; hA.YAxis.TickValues = [hA.YAxis.TickValues/5 hA.YAxis.TickValues hA.YAxis.TickValues*5]; end
set(gca,'FontSize',11)

subaxis(3,1,1,3)
if exist('logdata','var') && logdata == 1
    semilogy(0,0)
end
if exist('plot_gigs','var') && exist('gi_gs','var') && plot_gigs == 1
    for i=1:2:length(gi_gs)-1
        rectangle('Position',[gi_gs(i),min(x1),gi_gs(i+1)-gi_gs(i),max(x1)-min(x1)],'FaceColor', [0.73,0.73,0.73],'LineStyle','none'); hold on
    end
end
if exist('logdata','var') && logdata == 1
    semilogy(t1,x1,'k')
else
    plot(t1,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl3); ylim(yl);
if exist('updown','var') && updown==1; axis ij; end
xticks(xtck); ylabel(y_label); xlabel(x_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
if exist('logdata','var') && logdata == 1 && length(yticks)==1; hA.YAxis.TickValues = [hA.YAxis.TickValues/5 hA.YAxis.TickValues hA.YAxis.TickValues*5]; end
set(gca,'FontSize',11)

set(gcf,'PaperUnits','inches','PaperPosition',[-0.05 0 8.2 11.5])
print('-dpdf','-painters',[Fname '_' Vname '_KS'])

%% Save KSjumps dates

clear par_n
if exist('w_min','var'); par_n{1} = ['w_min = ' num2str(w_min)]; end
if exist('w_max','var'); par_n{end+1} = [', w_max = ' num2str(w_max)]; end
if exist('n_w','var'); par_n{end+1} = [', n_w = ' num2str(n_w)]; end
if exist('d_c','var'); par_n{end+1} = [', d_c = ' num2str(d_c)]; end
if exist('n_c','var'); par_n{end+1} = [', n_c = ' num2str(n_c)]; end
if exist('s_c','var'); par_n{end+1} = [', s_c = ' num2str(s_c)]; end
if exist('x_c','var'); par_n{end+1} = [', x_c = ' num2str(x_c)]; end

fileID = fopen([Fname '_' Vname '_KSjumps.csv'],'w');
fprintf(fileID,'%s\n',['# KS test parameters: ' strcat(par_n{:})]);
fprintf(fileID,'# The sign of a jump indicates direction: 1 = up, -1 = down, e.g. abrupt warming vs. abrupt cooling\n\n');
fprintf(fileID,'time,jump\n');
fprintf(fileID,'%g,%g\n',jumps');
fclose(fileID);