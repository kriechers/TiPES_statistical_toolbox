%% DESCRIPTION
% MATLAB script to detect abrupt transitions in the NGRIP d18O record
% (Rasmussen et al., 2014). This script identifies the transitions as well
% as the glacial/interglacial boundaries. The results are plotted and saved
% as a pdf file.
%
% Author: Witold Bagniewski Date: 31/07/2022
%
%% Load data
Fname = 'ngrip_d18o_20y';  ext = '.txt'; 
data=load([Fname ext]);
t=data(:,1);
x=data(:,2);
x1=x;

% % Data options
t = t/1000;                % convert yr to kyr
logdata = 0;               % default (no change to data)
% x = log10(x); logdata=1;   % logarithmic scale
updown = 0;                % default (no change to data)
% x = -x; updown=1;          % flip upside down

% % Plotting options
y_label = '\delta^{18}O [‰]';
x_label = 'Age (ky b2k)';
xtck = 0:2:5000;       % x ticks
xtck_m = 0:0.2:5000;   % small x ticks
LGC = 1;               % Plot only the last glacial cycle (7.8 to 122.3 ky b2k).
% LGC = 0;               % Plot the entire time series
plot_gigs = 1;         % show stadial/interstadial boundaries
% plot_gigs = 0;         % don't show stadial/interstadial boundaries

% % KS_detection parameters
min_w = 0.12;
max_w = 2.5;
n_w   = 12;
d_c   = 0.77;
n_c   = 3;
s_c   = 2.2;
x_c   = 0.8;

%% Run the augmented KS test

[jump_u,jump_d]=KS_detection(t,x,min_w,max_w,n_w,d_c,n_c,s_c,x_c);

%% Find "stadial - interstadial" boundaries (optional, takes a long time)

if plot_gigs == 1

clear mmmGIGS G_I G_S w_gigs min_wgigs max_wgigs n_wgigs;
min_wgigs = 5;   % min window size
max_wgigs = 60;  % max window size
n_wgigs = 10;    % number of windows

for i=1:n_wgigs
    w_gigs(i)=min_wgigs*(max_wgigs/min_wgigs)^((i-1)/(n_wgigs-1));   % makes n_w window sizes
end
kk=1;
for i=w_gigs
    [G_I(:,kk),G_S(:,kk)] = findgigs(t,x,i);
    mmmGIGS(:,kk)=(G_I(:,kk)+G_S(:,kk))/2;
    kk=kk+1;
end
mmmGIGS=mean(mmmGIGS');

gi_gs=[]; 
jump_d2=jump_d;
jump_d2(:,2)=-1;
jump_u2=jump_u;
jump_u2(:,2)=1;
jumps=sortrows([jump_u2;jump_d2]);

if mean(x(t>jumps(end,1)))>mmmGIGS(end)	% Check if it starts with glacial or interglacial
    mark=1; gi_gs=t(end); k=2; jumps(end+1,:)=[t(end),1];
else
    mark=0; k=1; jumps(end+1,:)=[t(end),0];
end

kk=1;   gi_gs_sm_dist = 0.8; 
for i = length(jumps):-1:1
    if i>1
        r2=find(t >= jumps(i-1,1) & t < jumps(i,1)); 
    else
        r2=find(t >= t(1) & t < jumps(i,1)); 
    end
    m2 = mean(x(r2));
    if m2 > mmmGIGS(max(r2)+1) && jumps(i,2) == 1 && mark == 0 
        gi_gs(k) = jumps(i,1); mark = 1; k=k+1;
    elseif m2 < mmmGIGS(max(r2)+1) && jumps(i,2) == -1 && mark == 1
        gi_gs(k) = jumps(i,1); mark = 0; k=k+1;
    end
end
if jumps(find(jumps(:,1)==gi_gs(end)),2) == 1
    gi_gs(end+1)=t(1);
end
gi_gs=flip(gi_gs);

else
    gi_gs = [];
end

%% Plotting

if LGC == 1
    xl1=[7.8,47.1]; xl2=[45.8,85.1]; xl3=[83.8,122.3];
else
    time3 = (max(t)-min(t))/3;
    xl1=[min(t),min(t)+time3*1.026];
    xl2=[min(t)+time3*0.987,min(t)+time3*2.013];
    xl3=[min(t)+time3*1.974,min(t)+time3*3];
end

yl=[min(x(x<=xl3(2))),max(x(x<=xl3(2)))];
figure; subaxis(3,1,1,1,'Spacing',0,'SpacingHoriz',0.3,'Margin',0.003,'MarginBottom',0.002,'MarginTop',0.001,'MarginLeft',0.075,'Paddingbottom',0.04);
for i=1:2:length(gi_gs)-1
    rectangle('Position',[gi_gs(i),-1000,gi_gs(i+1)-gi_gs(i),2000],'FaceColor',[.73 .73 .73],'LineStyle','none'); hold on
end
if logdata == 1
    semilogy(t,x1,'k')
else
    plot(t,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl1); ylim(yl);
if updown==1;axis ij;end
xticks(xtck); ylabel(y_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
set(gca,'FontSize',11)

subaxis(3,1,1,2)
for i=1:2:length(gi_gs)-1
    rectangle('Position',[gi_gs(i),-1000,gi_gs(i+1)-gi_gs(i),2000],'FaceColor',[.73 .73 .73],'LineStyle','none'); hold on
end
if logdata == 1
    semilogy(t,x1,'k')
else
    plot(t,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl2); ylim(yl); 
if updown==1;axis ij;end
xticks(xtck); ylabel(y_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
set(gca,'FontSize',11)

subaxis(3,1,1,3)
for i=1:2:length(gi_gs)-1
    rectangle('Position',[gi_gs(i),-1000,gi_gs(i+1)-gi_gs(i),2000],'FaceColor',[.73 .73 .73],'LineStyle','none'); hold on
end
if logdata == 1
    semilogy(t,x1,'k')
else
    plot(t,x1,'k')
end
vline(jump_u,'r');
vline(jump_d,'b');
xlim(xl3); ylim(yl);
if updown==1;axis ij;end
xticks(xtck); ylabel(y_label); xlabel(x_label); set(gca,'layer','top');
hA=gca; hA.XAxis.MinorTick = 'on'; hA.XAxis.MinorTickValues = xtck_m;
set(gca,'FontSize',11)

set(gcf,'PaperUnits','inches','PaperPosition',[-0.05 0 8.26 11.5])
print('-dpdf','-painters',[Fname '_KS'])