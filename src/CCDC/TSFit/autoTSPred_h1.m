function outfity=autoTSPred_h1(outfitx,fit_cft)
% Auto Trends and Seasonal Predict
% INPUTS:
% outfitx - Julian day [1; 2; 3];
% fit_cft - fitted coefficients;
% OUTPUTS:
% outfity - predicted reflectances [0.1; 0.2; 0.3];
% General model TSModel:
% f(x) =  a0 + b0*x + a1*cos(x*w) + b1*sin(x*w) 

% anual cycle
w=2*pi/365.25; 

outfity=[ones(size(outfitx)),outfitx,...% overall ref + trending
        cos(w*outfitx),sin(w*outfitx)] * fit_cft; % add seasonality
end
