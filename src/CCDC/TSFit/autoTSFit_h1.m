function [fit_cft,rmse,yhat]=autoTSFit_h1(x,y)
% Using lasso for timeseries modeling (01/27/2013)
% Auto Trends and Seasonal Fit between breaks
% INPUTS:
% x - Julian day [1; 2; 3];
% y - predicted reflectances [0.1; 0.2; 0.3];
%
% OUTPUTS:
% fit_cft - fitted coefficients;
% General model TSModel:
% f(x) =  a0 + b0*x + a1*cos(x*w) + b1*sin(x*w)

n=length(x); % number of clear pixels
w=2*pi/365.25; % anual cycle
df=4; %degree of freedom

% build X
X=zeros(n,df-1);
X(:,1)=x;

X(:,2)=cos(w*x);
X(:,3)=sin(w*x);

fit=glmnet(X,y);
fit_cft=[fit.a0;fit.beta];

yhat=autoTSPred_h1(x,fit_cft);
rmse=norm(y-yhat)/sqrt(n);
end
