function [G_I,G_S] = findgigs(year,proxy,window)
%% DESCRIPTION
% A function to extract "stadial" and "interstadial" time series from a
% paleoproxy time series. For each sampling window, the samples are divided
% into two sets, below and above the mean value of the series. The
% "stadial" values are then calculated from the 25th percentile of the set
% below the mean, while "interstadial" values are calculated from the 75th
% percentile of the set above the mean.
%
% Author: Witold Bagniewski
% Date: 31/07/2022 
%%
intxx=min(year):(max(year)-min(year))/(length(year)-1):max(year);
intyy=interp1(year,proxy,intxx);
for i = 1:length(intxx)
    r1=find(intxx >= intxx(i)-window/2 & intxx < intxx(i)+window/2);	% window around
    xxx=intyy(r1);
    Gii(i) = quantile(xxx(xxx >= mean(xxx)),0.75); 
    Gss(i) = quantile(xxx(xxx <= mean(xxx)),0.25); 
end
G_I=interp1(intxx,Gii,year);
G_S=interp1(intxx,Gss,year);