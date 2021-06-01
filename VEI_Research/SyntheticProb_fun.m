%% Function for creating the synthetic data for the baseline 

function [VEI_synth]=SyntheticProb_fun(VEIdata,NumVolc);
% NumVolc = number of synthetic VEI to produce
% VEIdata= a list of VEIs from actual data to define distribution

% the number of bins change depending on the threshold, is this okay that
% it is always 0:7? 
x_num=[0:7];
[HVEI,EVEI]=hist(VEIdata,x_num); %histogram

pdfVEI=HVEI./(sum(HVEI)); % probability density function

cdfVEI=cumsum(pdfVEI); % cumulative distribution function


%Creating a Four loop to produce values that can be used for a histograph
%of synthetic data point based on the cdf of the Last VEI
p=rand(NumVolc,1);
VEI_synth1=[];
for i=1:length(p)
    I=find(p(i)<=cdfVEI,1,'first');
    VEI_synth(i)=EVEI(I);
end
end

