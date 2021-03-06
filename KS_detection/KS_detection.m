function [jump_u,jump_d] = KS_detection(t,x,min_w,max_w,varargin) 
%% DESCRIPTION
% The KS_detection function detects abrupt transitions (jumps) in climate
% proxy records using an augmented Kolmogorov–Smirnov (KS) test. This
% function is based on methodology described in the paper:
% Bagniewski, W., Ghil, M., Rousseau, D.-D., (2021) Automatic detection of
% abrupt transitions in paleoclimate records. Chaos, 31(11):113129,
% https://doi.org/10.1063/5.0062543
%
% Author: Witold Bagniewski
% Date: 31/07/2022 
%
%% USAGE 
% Syntax
% [jump_u,jump_d] = KS_detection(t,x,min_w,max_w,n_w,d_c,n_c,s_c,x_c)
%
% Required Inputs:
% t     : Time. Ordered from youngest to oldest.
% x     : Proxy value. Note: t and x must have the same length.
% min_w : Size of the smallest window
% max_w : Size of the largest window
%
% Optional Inputs:
% n_w   : Number of distinct window lengths
% d_c   : Cut-off threshold D_c for the KS statistic
% n_c   : Minimum sample size threshold n_c
% s_c   : Cut-off threshold s_c for the standard deviation
% x_c   : Cut-off threshold x_c for the change in proxy value
%
% Outputs:
% jump_u : Dates of jumps that go "up", e.g. abrupt warming events
% jump_d : Dates of jumps that go "down", e.g. abrupt cooling events
%
%% EXAMPLE
% % Generate some data:
% tt=1:5000;
% xx=0;
% for i=2:5000;
%     xx(i)=xx(i-1)+rand*1.3-0.65;
%     xx(i)=abs(xx(i))^0.5*sign(xx(i));
% end
% 
% % Find jumps:
% [jump_ux,jump_dx] = KS_detection(tt,xx,20,200,5);
% 
% % Plot the data and the jumps:
% figure; plot(tt,xx,'k'); hold on; vline(jump_ux,'r'); vline(jump_dx,'b');

%% Default values

n_w = 10;
d_c = 0.75;
n_c = 3;
s_c = 1.5;
x_c = std(x)*0.15 + (max(x)-min(x))*0.03;

try
    n_w = varargin{1};
    d_c = varargin{2};
    n_c = varargin{3};
    s_c = varargin{4};
    x_c = varargin{5};
end

%% Preprocessing the time series

data(:,1) = t;
data(:,2) = x;

% For dates that are equal find the mean value of x
[C,~,idx] = unique(data(:,1),'stable');
val = accumarray(idx,data(:,2),[],@mean);
data1 = [C val];  clearvars C ia idx val 

% Delete NaN rows
data1(isnan(data1(:,1)),:) = [];
data1(isnan(data1(:,2)),:) = [];

tt = data1(:,1);
xx = data1(:,2);

% Make sure all values have the same sign
if sign(max(xx))>sign(min(xx))
    xx=xx-min(xx)+0.1;
end

%% Varying window size

if or(min_w == max_w , n_w==1)   % in case there is only one n_w window size
    kswindow = min_w;
    n_w = 1;
else
    for i=1:n_w
        kswindow(i)=min_w*(max_w/min_w)^((i-1)/(n_w-1));   % makes n_w window sizes
    end
end

%% Run the KS test

for j = 1:n_w
    for i = 1:length(tt)
        
        % locate dates for two sliding windows
        y=tt(i);
        r1=find(tt > y-kswindow(j) & tt <= y);
        r2=find(tt > y & tt <= y+kswindow(j));
        
        % run KS test for the two windows
        if ~isempty(r1) && ~isempty(r2) && min(tt)<=y-kswindow(j) && max(tt)>=y+kswindow(j)
            [~,~,c]=kstest2(xx(r1),xx(r2));
            cc(i)=c;            
        else
            cc(i)=0;
        end
        
        % calculate additional variables
        changes(i)=mean(xx(r1))-mean(xx(r2));
        len1(i)=length(xx(r1)); 
        len2(i)=length(xx(r2));
        std1(i)=std(xx(r1)); 
        std2(i)=std(xx(r2));
    end
    
    ksstat(:,j)=cc.*sign(changes);     % KS statistic (D_KS)+ sign
    change_ks(:,j)=changes;         % difference between the two samples
    kslen1(:,j)=len1;               % size of sample1
    kslen2(:,j)=len2;               % size of sample2
    ksstd1(:,j)=std1;               % st. dev. of sample1
    ksstd2(:,j)=std2;               % st. dev. of sample2
end

%% Find the jumps

% Apply d_c, n_c, s_c, x_c threshold parameters
clearvars -except kswindow data1 tt xx n_w ksstat change_ks kslen1 kslen2 ksstd1 ksstd2 d_c n_c s_c x_c
th_si = kslen1>=n_c & kslen2>=n_c;	% sample size threshold
th_st = s_c*ksstd1<=abs(change_ks) & s_c*ksstd2<=abs(change_ks);  % sample std threshold
th_change = x_c<=abs(change_ks); 
ksstat2 = 1-(1-abs(ksstat))./(1-((kslen1+kslen2)./(kslen1.*kslen2)).^0.5); 
ksstat2(ksstat2 < 0)=0; 
ksstat2 = ksstat2.*sign(change_ks); 
ksstat2([1,end-1],:)=0;
th_ks = abs(ksstat2)>=d_c;
ks_all = ksstat2.*th_st.*th_si.*th_change.*th_ks;

% Scale by mean D_KS of smaller windows
if length(ks_all(1,:))>2
    for i=length(ks_all(1,:)):-1:3
        ks_all(:,i)=ks_all(:,i).*0.67+mean(ksstat2(:,1:i-1)')'.*0.33;
    end
    ks_all(:,2)=ks_all(:,2).*0.67+ksstat2(:,1).*0.33;
end

% Append mean of all windows
ksstat_m=ksstat; 
if length(ks_all(1,:))>1
    ksstat_m(:,end+1) = mean(ksstat_m')';
    ksstat2(:,end+1) = mean(ksstat2')';
    ks_all(:,end+1) = mean(ks_all')';
    ks_all(abs(ks_all)<d_c) = 0;
end

% Find D_KS peaks
for j=1:length(ksstat2(1,:))
    k=1; kk=1;                                           % k = index of sign changes, kk = index of peaks
    for i=2:length(ksstat2(:,1))
        if sign(ksstat_m(i-1,j)) ~= sign(ksstat_m(i,j))          % check if sign of ks_all changes
            if any(ks_all(k:i-1,j) ~= 0) && any(~isnan(ks_all(k:i-1,j)))
                mmm=ks_all(k:i-1,j); 
                [~,I] = max(abs(mmm));                   % maximum of ks_all2
                iii=find(mmm==mmm(I))+k-1;  
                [~,II] = max(abs(sum(ksstat2(iii,:)')));  % maximum of ks_all, in case ks_all2 are equal
                ks_peak{j}(kk,2)=ks_all(iii(II(end)),j);
                ks_peak{j}(kk,1)=iii(II(end));
                kk=kk+1;
            end
            k=i;
        end
    end
end
peak_up=ks_peak;
peak_do=ks_peak;
for i=1:length(ks_peak)
    peak_up{i}(peak_up{i}(:,2)<0,:)=[];  % Split ks_peak into positive and negative
    peak_do{i}(peak_do{i}(:,2)>0,:)=[];  % Split ks_peak into positive and negative
end

% delete jumps found with smaller windows that correspond to those found with larger windows
jump_u=(tt(peak_up{length(ks_all(1,:))}(:,1))+tt(peak_up{length(ks_all(1,:))}(:,1)+1))/2;  % use middle date between sample1 and sample2 as the transition date 
jump_d=(tt(peak_do{length(ks_all(1,:))}(:,1))+tt(peak_do{length(ks_all(1,:))}(:,1)+1))/2;
if length(ks_all(1,:))>1
    for j=length(ks_all(1,:))-1:-1:1
        if ~isempty(peak_up{j})
            xxx_u=(tt(peak_up{j}(:,1))+tt(peak_up{j}(:,1)+1))/2;
        else
            xxx_u=[];
        end
        if ~isempty(peak_do{j})
            xxx_d=(tt(peak_do{j}(:,1))+tt(peak_do{j}(:,1)+1))/2;
        else
            xxx_d=[];
        end
        for i=1:length(jump_u)
            xxx_u(xxx_u<jump_u(i)+kswindow(j) & xxx_u>jump_u(i)-kswindow(j))=[];   % delete peaks that are close to those found with larger windows
        end
        for i=1:length(jump_d)
            xxx_d(xxx_d<jump_d(i)+kswindow(j) & xxx_d>jump_d(i)-kswindow(j))=[];
        end
        jump_u=[jump_u; xxx_u];
        jump_d=[jump_d; xxx_d];
    end
end
jump_u=sort(jump_u);
jump_d=sort(jump_d);

end