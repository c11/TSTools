function rec_cg=TrendSeasonalFit_v9_QGIS_max(sdate, line_t, n_times, conse, T_cg, num_c, B_detect)
% function rec_cg=TrendSeasonalFit_v9_QGIS_max(sdate, line_t, n_times, conse, T_cg, num_c, B_detect)
% CCDC 9.3+ version - Zhe Zhu, Boston University
% It is based on 7-bands fitting for Iterative Seasonal, Linear, and Break Models
% This function works for analyzing one line of time series pixel
%
% Revisions: $ Date: 03/31/2013 $ Copyright: Zhe Zhu
%
% MODIFICATIONS FROM CCDCv9.3:
%   Use FMask 3.2 snow mask
%   Only fit more than 50% of Landat images overlap area (08/28/2013)
%   Fixed bugs in calculating temporally adjusted rmse (08/01/2013)
%
% Fit curve again after one year (03/28/2013)
% Use mini rmse = T_const/T_cg for small rmse cases (03/26/2013)
% Remove out of range pixels before time series analysis (02/09/2013)
% Using 8 coefficients instead of 4 coefficients (02/01/2013)
% Using lasso for timeseries modeling: "V8->V9" (01/27/2013)
% Use max v_slope instead of average v_slope  (01/16/2013)
% Start initialization when time_span > 1 year (01/16/2013)
% Fixed bugs in not fitting models at the begining (01/16/2013)
% Do not show clouds/cloud shadows/snows in the plot (01/15/2013)
% Fixed bugs in counting "i" and "i_span"(01/13/2013)
% Temporally changing RMSE: "V7->V8" (01/09/2013) 
% Fix bugs in model intialization (07/11/2012)
% Check abnormal pixels at the end of model intialization (07/01/2012)
% Including thermal band for change detectiona & classification (03/25/2012)
% Basic requrirements only has i_span (03/13/2012)
% Check abnormal pixels at the begining of the model initialization (03/09/2012)
% Use only 2 parameters for capturing phenology (03/08/2012)
% Use Z-score statistical thresholds (03/01/2012)
%
% Inputs:
% stk_n='calmatch'; stack image name
% ncols = 8021; % number of pixels processed per line
% nrows=1; % the nrowsth lines
% for example    1 2 3 4 5
%                6 7 8 9 10
%
% Outputs:
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
%
fprintf('\n\n');
fprintf('Running CCDC v9.3+ with options:\n');
fprintf('   n_times:    %d\n', n_times);
fprintf('   conse:      %i\n', conse);
fprintf('   T_cg:       %d\n', T_cg);
fprintf('   num_coef:   %i\n', num_c);
disp(strrep(['   B_detect:  [' sprintf(' %i,', B_detect) ']'], ',]', ' ]'))
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
% Vary number of coefficients
if num_c == 2
    autoTSFit = @autoTSFit_h0;
    autoTSPred = @autoTSPred_h0;
elseif num_c == 4
    autoTSFit = @autoTSFit_h1;
    autoTSPred = @autoTSPred_h1;
elseif num_c == 6
    autoTSFit = @autoTSFit_h2;
    autoTSPred = @autoTSPred_h2;
elseif num_c == 8
    autoTSFit = @autoTSFit_h3;
    autoTSPred = @autoTSPred_h3;
end

%% Constants
% NUM of Functional Curves (num_fc)
num_fc=0;
% number of days per year
num_yrs=365.25;
% number of bytes: int16
num_byte=2;
% number of bands
nbands = 8;

%%% total bands (1-5,7,6,Fmask)
% Band for multitemporal cloud/snow detection (Green)
num_B1=2;
% Band for multitemporal shadow/snow shadow detection (SWIR)
num_B2=5;
% Data bands used for change test
num_detect = length(B_detect);

% Threshold for cloud, shadow, and snow detection.
T_const=400;
% mininum rmse
% mini = T_const / T_cg;
% CEHOLDEN change:  10.9 mini RMSE
mini = 150;
% mininum years
mini_yrs=1;

% initialize the struct data of RECording of ChanGe (rec_cg)
rec_cg=struct('t_start',[], ...
    't_end', [], ...
    't_break', [], ...
    'coefs', [], ...
    'rmse', [], ...
    'pos', [], ...
    'band', [], ...
    'v_dif', []);

    % mask data
    line_m=line_t(:, nbands);

    % CEHOLDEN change: code upgrade to skip pixels with <50% of data
    % Only run CCDC for places where more than 50% of images has data
    idexist = line_m < 255;
    overlap_pct = 100 * sum(idexist) / size(sdate, 1);
    if overlap_pct < 50
        fprintf('At least 50% of pixel is NODATA\n');
        return;
    end

    % CEHOLDEN change: FMask 3.2 snow mask
    % improved snow mask (will be removed if using Fmask 3.2 version 
    line_m((line_t(:, 2) - ... % B2 - B5
        line_t(:, 5)) ./ ... % over
        (line_t(:, 2) + ... % B2 + B5
        line_t(:, 5)) > 0.15 & ... % NDSI > 0.15
        line_t(:, 4) > 1100 & ... % B4 > 1100
        line_t(:, 2) > 1000 & ... % B2 > 1000
        line_m <= 1) = 3; % and is already clear --> set to snow

    % clear land or water
    idclear=line_m==0|line_m==1;
    % clear pixel should have reflectance between 0 and 1 and temperature 
    % between -89 to 58 celsius degree
    idclear=idclear & ...
        line_t(:, 1) > 0 & line_t(:, 1) < 10000 & ...
        line_t(:, 2) > 0 & line_t(:, 2) < 10000 & ...
        line_t(:, 3) > 0 & line_t(:, 3) < 10000 & ...
        line_t(:, 4) > 0 & line_t(:, 4) < 10000 & ...
        line_t(:, 5) > 0 & line_t(:, 5) < 10000 & ...
        line_t(:, 6) > 0 & line_t(:, 6) < 10000 & ...
        line_t(:, 7) > -8900 & line_t(:, 7) < 5800;
    
    % Xs & Ys for computation
    clrx=sdate(idclear);
    clry=line_t(idclear, 1:nbands - 1);
    
    % defining computed variables
    fit_cft=zeros(num_c,nbands-1);
    
    % start with the miminum requirement of clear obs
    i=n_times*num_c;
    
    % enough clear observation for fitting or not
    if length(clrx)<i+conse % not enough clear obs
        fprintf('Not enough clear observations\n');
        return;
    else
        % initializing variables
        % the first observation for TSFit
        i_start=1;
        % record the start of the model(0=>initial;1=>done)
        BL_train=0;
        % identified and move on for the next curve
        num_fc=num_fc+1; % NUM of Fitted Curves (num_fc)
    end
    
    % while loop - process till the last clear observation - conse
    while i<=length(clrx)-conse
        % span of "i"
        i_span=i-i_start+1;
        % span of time (num of years)
        time_span=(clrx(i)-clrx(i_start))/num_yrs;
        
        % basic requrirements: 1) enough observations; 2) enough time
        if i_span>=n_times*num_c&&time_span>=mini_yrs
            % initializing model
            if BL_train==0
                % Step 1: noise removal
                % multi-temporal cloud, cloud shadow and snow masking
                
                % blIDs = [0 1]; 0 => clear & 1 => cloud/snow/shadow
                blIDs=autoMask(clrx(i_start:i+conse), ...
                        clry(i_start:i+conse,[num_B1,num_B2]), ...
                        time_span,T_const);
                
                % update i_span after noise removal
                i_span=sum(~blIDs(1:end-conse));
                
                % check if there is enough observation
                if i_span<n_times*num_c
                    % move forward to the i+1th clear observation
                    i=i+1;
                    % not enough clear observations
                    continue;
                else
                    IDs=i_start:i;% all ids
                    rmIDs=IDs(blIDs(1:end-conse)==true);% IDs to be removed
                    
                    % copy x & y
                    cpx=clrx;
                    cpy=clry;
                    
                    % remove noise pixels between i_start & i
                    cpx(rmIDs)=[];
                    cpy(rmIDs,:)=[];
                    
                    % record i before noise removal
                    % This is very important as if model is not initialized
                    % the multitemporal masking shall be done again instead
                    % of removing outliers in every masking
                    i_rec=i;
                    
                    % update i afer noise removal (i_start stays the same)
                    i=i_start+i_span-1;
                    % update span of time (num of years)
                    time_span=(cpx(i)-cpx(i_start))/num_yrs;
                    
                    % check if there is enough time
                    if time_span < mini_yrs
                        % keep the original i
                        i=i_rec;
                        % move forward to the i+1th clear observation
                        i=i+1;
                        % not enough time
                        continue;
                        % Step 2: model fitting
                    else
                        % initialize model testing variables
                        % CEHOLDEN change:  store test variables for all bands
                        v_dif = zeros(size(B_detect, 2)); % summation of all tests divided by rmse
                        
%                        for i_B=1:nbands-1
                        % CEHOLDEN change:  move to B_detect subset of all bands
                        for i_B = 1:size(B_detect, 2)
                            % initial model fit
                            [fit_cft(:, B_detect(i_B)), rmse] = ...
                                autoTSFit(cpx(i_start:i), ...
                                cpy(i_start:i, B_detect(i_B)));
                            % prevent ideal fit
                            mini_rmse = max(rmse, mini);

                            % CEHOLDEN change:  update how slope is dealt with
                            v_slope = abs(fit_cft(2, B_detect(i_B)) * clrx(i) - clrx(i_start));
                            % compare the first clear obs
                            v_start = abs(clry(i_start, B_detect(i_B)) - ...
                                autoTSPred(clrx(i_start), fit_cft(:, B_detect(i_B))));
                            % compare the last clear observation
                            v_end = abs(clry(i, B_detect(i_B)) - ...
                                autoTSPred(clrx(i), fit_cft(:, B_detect(i_B))));

                            % Total differences
                            v_dif(i_B) = (v_slope + v_start + v_end) / mini_rmse;
                        end
                        
                        % CEHOLDEN change:  use maximum of all bands
                        % find stable start for each curve
                        if max(v_dif) > T_cg 
                            % start from next clear obs
                            i_start=i_start+1;
                            % keep the original i
                            i=i_rec;
                            % move forward to the i+1th clear observation
                            i=i+1;
                            % keep all data and move to the next obs
                            continue;
                        else
                            % model ready!
                            BL_train=1;
                            % count difference of i for each iteration
                            i_count=0;
                            % remove noises
                            clrx=cpx;
                            clry=cpy;
                        end
                    end
                end
            end % end of initializing model
            
            % continuous monitoring started!!!
            if BL_train==1
                % all IDs
                IDs=i_start:i;
                
                % dynamic model estimation
                % if (i-i_start+1)-i_count >= (i-i_start+1)/3
                if clrx(i)-clrx(i_start) >= i_count + num_yrs   
                    rmse = zeros(nbands - 1, 1);
                    for i_B = 1:(nbands - 1)
                        [fit_cft(:, i_B), rmse(i_B)] = ...
                            autoTSFit(clrx(IDs), clry(IDs, i_B));
                    end
                    % update i_count at each interation
                    % i_count=i-i_start+1;
                    i_count = clrx(i) - clrx(i_start);
                    
                    % updating information at each iteration
                    % record time of curve start
                    rec_cg(num_fc).t_start = clrx(i_start);
                    % record time of curve end
                    rec_cg(num_fc).t_end = clrx(i);
                    % record break time
                    rec_cg(num_fc).t_break = 0; % no break at the moment
                    % record postion of the pixel
                    rec_cg(num_fc).pos = 0;
                    % record fitted coefficients
                    rec_cg(num_fc).coefs = fit_cft;
                    % record rmse of the pixel
                    rec_cg(num_fc).rmse = rmse;
                    % record bands exceeding T_cg (none for now)
                    rec_cg(num_fc).bands = 0;
                    % record critical value matrix (all 0 for now)
                    rec_cg(num_fc).v_dif = zeros(conse, size(B_detect, 2));
                else
                    % update information without iteration
                    % record time of curve end
                    rec_cg(num_fc).t_end=clrx(i);
                end
                
                % change format to year/month/day
                vec_x=datevecmx(clrx(IDs));
                % value of difference for conse observations
                % CEHOLDEN change:  set v_dif as matrix to store conse
                %                   values for all bands
                v_dif=zeros(conse, size(B_detect, 2));
                % num of observations for calculating RMSE
                n_rmse=n_times*num_c;
                
                for i_conse=1:conse
                    % sorted index of the kth minimum distant dates
                    sorted_indx=zeros(n_rmse,1);
                    % vector format of i+conse julian date
                    vec_j=datevecmx(clrx(i+conse));
                    % convert julian date to date of day and then the 
                    % difference of days
                    % CEHOLDEN bugfix:  temporally changing RMSE fit
                    %                   (i+conse) to (i+i_conse)
                    d_yr=abs((clrx(IDs) - datenummx(vec_x(:,1), 1, 0)) ...
                            - (clrx(i + i_conse) - datenummx(vec_j(1), 1, 0)));
                    
                    % get the list of n_rmse smallest elements
                    for j=1:n_rmse
                        [~,sorted_indx(j)]=min(d_yr);
                        d_yr(sorted_indx(j))=num_yrs; % give max value to exclude the smallest elements
                    end
                    
                    % CEHOLDEN change:  move to B_detect subset of all bands
%                    for i_B=1:nbands-1
                    for i_B = 1:size(B_detect, 2)
                        % temporally changing RMSE
                        rmse=norm(clry(IDs(sorted_indx), B_detect(i_B)) - ...
                                autoTSPred(clrx(IDs(sorted_indx)), ...
                                    fit_cft(:, B_detect(i_B)))) ...
                                / sqrt(n_rmse);
                        
                        % prevent ideal fit
                        rmse = max(rmse, mini);
                        
                        % difference of all spectral bands for conse observations
                        v_dif(i_conse, i_B) = abs(clry(i+i_conse, B_detect(i_B)) - ... 
                            autoTSPred(clrx(i + i_conse), fit_cft(:, B_detect(i_B)))) ...
                            / (T_cg * rmse);
                    end
                end
                
                % average spectral difference of three more clear obs
%                v_dif = v_dif / num_detect;
%                 v_dif=v_dif / (nbands - 1);
                
                % CEHOLDEN change: get max for each consecutive for any band
                max_v_dif = max(v_dif, [], 2);

                % change detected
                if min(max_v_dif) > 1
                    % t_break > 0 record break time
                    rec_cg(num_fc).t_break = clrx(i + 1);

                    % CEHOLDEN change: record which bands exceeded critical value
                    [~, band_crit] = find(v_dif > 1);
                    rec_cg(num_fc).band = unique(B_detect(band_crit));

                    % CEHOLDEN change: store the v_dif matrix
                    rec_cg(num_fc).v_dif = v_dif;

                    % identified and move on for the next functional curve
                    num_fc = num_fc + 1;
                    % start from i+1 for the next functional curve
                    i_start = i + 1;
                    % start training again
                    BL_train = 0;
                    % false change
                elseif max_v_dif(1) > 1
                    % remove noise
                    clrx(i + 1,:)=[];
                    clry(i + 1,:)=[];
                    i = i - 1; % stay & check again after noise removal
                end
            end % end of continuous monitoring
        end % end of checking basic requrirements
        
        % move forward to the i+1th clear observation
        i=i+1;
    end % end of while iterative

end % end of function
