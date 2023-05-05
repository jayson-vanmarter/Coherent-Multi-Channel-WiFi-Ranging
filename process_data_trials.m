%% Copyright(C) 2023 The University of Texas at Dallas
% Developed by: Jayson P. Van Marter
% Advisor: Prof. Murat Torlak
% Department of Electrical and Computer Engineering
%
%  This work was supported by the Semiconductor Research Corporation (SRC)
%  task 2810.076 through The University of Texas at Dallas' Texas Analog
%  Center of Excellence (TxACE).
%
%  Redistributions and use of source must retain the above copyright notice

%%
% This script evaluates ranging performance for a given set of trial
% locations, a given packet type, and a given bandwidth. If the bandwidth
% is 20 MHz, multi-channel stitches are also evaluated.

clear
close all

%% Set the file paths

% File path for the data sets
mainpath = '.\savedata\';

% Path for the calibration measurements
filepath_cal = '.\savedata\Calibration\';

% Get the trial names and true distances for a given data set or across
% multiple data sets
[trialnames1,dist_true1] = get_trialdir_info("lab_room");
[trialnames2,dist_true2] = get_trialdir_info("roomtoroom_LOS");
[trialnames3,dist_true3] = get_trialdir_info("classroom");
[trialnames4,dist_true4] = get_trialdir_info("indoor_open_space");
trialnames = [trialnames1;trialnames2;trialnames3;trialnames4]; % Concatenate trial names
dist_true = [dist_true1;dist_true2;dist_true3;dist_true4]; % Concatenate distances
% trialnames = trialnames1(1:3); % Test a subset of trials
% dist_true = dist_true1(1:3); % Get the corresponding subset of distances

num_subimp = 2; % One LTF for HE/VHT and one LTF for Legacy for each packet
td = 2; % Number of measurements across time per LTF type per packet
num_locs = length(trialnames); % Number of unique locations evaluated

%% Set packet type and bandwidth to evaluate.
% Packet Types: HESU (High Efficiency Single User) or VHT (Very High
% Throughput)
% Bandwidths: 20, 40, 80, 160 MHz

% filename = 'HESU_20'; % HESU, 20 MHz BW
% filename = 'HESU_40'; % HESU, 40 MHz BW
% filename = 'HESU_80'; % HESU, 80 MHz BW
% filename = 'HESU_160'; % HESU, 160 MHz BW
% filename = 'VHT_20'; % VHT, 20 MHz BW
% filename = 'VHT_40'; % VHT, 40 MHz BW
% filename = 'VHT_80'; % VHT, 80 MHz BW
filename = 'VHT_160'; % VHT, 160 MHz BW

% Get the corresponding bandwidth in MHz
BW = str2double(erase(filename,{'HESU_','VHT_'})); 

%% Set Parameters

% Define the propagation speed
c = 3e8; % For free space measurements
% c = 3e8 / sqrt(1.4); % For wired measurements (dielectric_const = 1.4)

% Define the signal subspace sizes for MUSIC and the downsampling rates
% depending on the bandwidth.
[M,M_320,M_745,...
    Dr_VHT,Dr_HE,Dr_Legacy,...
    Dr_VHT_320,Dr_HE_320,Dr_Legacy_320,...
    Dr_VHT_745,Dr_HE_745,Dr_Legacy_745] = ...
    get_MandDr(BW);

%% Evaluate trial location for the given packet type and bandwidth

% For each location
for p = 1:num_locs

    % Load in the data for the given implementation and bandwidth
    sve = load([mainpath,trialnames{p},'\',filename]);
    sve_cal = load([filepath_cal,filename]);
    sve.dist_true = dist_true(p); % Set the true distance

    % Set parameters
    [Delta_f,~,~,~,~,...
        Delta_f_LLTF,~,~,~,~] = get_impparams(sve);
    num_ch = length(sve.fc); % Number of channels (trials) scanned for the bandwidth

    % Set downsampling rates
    switch Delta_f
        case 78.125e3 % HE-LTF case
            Dr = Dr_HE;
            Dr_320 = Dr_HE_320;
            Dr_745 = Dr_HE_745;
        case 312.5e3 % VHT-LTF case
            Dr = Dr_VHT;
            Dr_320 = Dr_VHT_320;
            Dr_745 = Dr_VHT_745;
    end

    % Format, calibration, and interpolate for null subcarriers within each
    % channel
    [chanEst,chanEst_LLTF,fInterp,fInterp_LLTF,SNR_est,SNR_est_LLTF] = ...
        format_calibrate(sve,sve_cal,c);
    
    % Loop once for HE/VHT (chanEst), then once for Legacy (chanEst_LLTF)
    for u = 1:num_subimp

        % Set variables depending on the LTF type
        if u == 1 % HE or VHT case
            chanEst_in = chanEst;
            fInterp_in = fInterp;
            SNR_est_in = SNR_est;
            Delta_f_in = Delta_f;
            Dr_in = Dr;
            Dr_320_in = Dr_320;
            Dr_745_in = Dr_745;
        elseif u == 2 % Legacy case
            chanEst_in = chanEst_LLTF;
            fInterp_in = fInterp_LLTF;
            SNR_est_in = SNR_est_LLTF;
            Delta_f_in = Delta_f_LLTF;
            Dr_in = Dr_Legacy;
            Dr_320_in = Dr_Legacy_320;
            Dr_745_in = Dr_Legacy_745;
        end

        % Loop for each individual channel
        for k = 1:num_ch

            % Estimate distance per time measurement per channel
            dist_est = nan(td,1);
            for j = 1:td
                dist_est_all(p,k,u,j) = singlech_distest(...
                    fInterp_in(:,k),squeeze(chanEst_in(:,k,:,j)),Delta_f_in,M,c,Dr_in);
            end
    
            % Average SNR estimates across antennas and time for the
            % channel
            SNR_est_ch(p,k,u) = mean(SNR_est_in(k,:,:),'all'); 

            % Console output
            disp(['Location ',num2str(p),'/',num2str(num_locs),...
                ', Channel ',num2str(k),'/',num2str(length(sve.fc)),...
                ', Dist Errs: ',num2str(dist_est_all(p,k,u,1)-dist_true(p)),', ',num2str(dist_est_all(p,k,u,2)-dist_true(p))]);
        
        end
    
        % Estimate distance per set of CFRs measured in time
        for j = 1:td
            if sve.fs == 20e6 % Only do channel stitching for 20 MHz channel measurements
                dist_est_chstitch320_all(p,u,j) = chstitch320M_distest(...
                    fInterp_in,squeeze(chanEst_in(:,:,:,j)),Delta_f_in,M_320,c,Dr_320_in); % 320 MHz multi-channel stitch
                dist_est_chstitch745_all(p,u,j) = chstitch745M_distest(...
                    fInterp_in,squeeze(chanEst_in(:,:,:,j)),Delta_f_in,M_745,c,Dr_745_in); % 745 MHz multi-channel stitch
            else
                dist_est_chstitch320_all(p,u,j) = nan;
                dist_est_chstitch745_all(p,u,j) = nan;
            end
        end

        % Console output
        if sve.fs == 20e6
            disp('---')
            disp(['Location ',num2str(p),'/',num2str(num_locs),...
                ', 320 MHz Stitch Dist Errs: ',num2str(dist_est_chstitch320_all(p,u,1)-dist_true(p)),...
                ', ',num2str(dist_est_chstitch320_all(p,u,2)-dist_true(p))]);
            disp(['Location ',num2str(p),'/',num2str(num_locs),...
                ', 745 MHz Stitch Dist Errs: ',num2str(dist_est_chstitch745_all(p,u,1)-dist_true(p)),...
                    ', ',num2str(dist_est_chstitch745_all(p,u,2)-dist_true(p))]);
                
            disp('---')
        end

    end

end

%% Process the results and compute statistics

% Separate out the HE or VHT and Legacy measurements
dist_est = squeeze(dist_est_all(:,:,1,:)); % HE or VHT
dist_est_legacy = squeeze(dist_est_all(:,:,2,:)); % Legacy
dist_est_chstitch320 = squeeze(dist_est_chstitch320_all(:,1,:)); % HE or VHT
dist_est_chstitch320_legacy = squeeze(dist_est_chstitch320_all(:,2,:)); % Legacy
dist_est_chstitch745 = squeeze(dist_est_chstitch745_all(:,1,:)); % HE or VHT
dist_est_chstitch745_legacy = squeeze(dist_est_chstitch745_all(:,2,:)); % Legacy

% Calculate error
dist_err = dist_est - dist_true;
dist_err_legacy = dist_est_legacy - dist_true;
dist_err_chstitch320 = dist_est_chstitch320 - dist_true;
dist_err_chstitch320_legacy = dist_est_chstitch320_legacy - dist_true;
dist_err_chstitch745 = dist_est_chstitch745 - dist_true;
dist_err_chstitch745_legacy = dist_est_chstitch745_legacy - dist_true;

% Calculate RMSE
dist_rmse = sqrt(mean(dist_err(:).^2));
dist_rmse_legacy = sqrt(mean(dist_err_legacy(:).^2));
dist_rmse_chstitch320 = sqrt(mean(dist_err_chstitch320(:).^2));
dist_rmse_chstitch320_legacy = sqrt(mean(dist_err_chstitch320_legacy(:).^2));
dist_rmse_chstitch745 = sqrt(mean(dist_err_chstitch745(:).^2));
dist_rmse_chstitch745_legacy = sqrt(mean(dist_err_chstitch745_legacy(:).^2));

%% Load saved data

% LOS trials
% load('.\savedata\HESU_20_AllLOS')
% load('.\savedata\HESU_40_AllLOS')
% load('.\savedata\HESU_80_AllLOS')
% load('.\savedata\HESU_160_AllLOS')
% load('.\savedata\VHT_20_AllLOS')
% load('.\savedata\VHT_40_AllLOS')
% load('.\savedata\VHT_80_AllLOS')
% load('.\savedata\VHT_160_AllLOS')

%% Plot CDFs
close all

p_lim = 2; % Plot x-axis limit

% HE or VHT CDFs per bandwidth -----------------------------------------%
figure(1);
subplot(2,1,1);
hold on
cdfp(1) = cdfplot(dist_err(:));
if sve.fs == 20e6
    cdfp(2) = cdfplot(dist_err_chstitch320(:));
    cdfp(3) = cdfplot(dist_err_chstitch745(:));
end
set(cdfp,'linewidth',2)
hold off
xlim([-p_lim p_lim])
grid on
title([erase(filename,"_"),' CDFs'],'interpreter','tex')
xlabel('Error (m)')
ylabel('CDF')

subplot(2,1,2);
hold on
abscdfp(1) = cdfplot(abs(dist_err(:)));
if sve.fs == 20e6
    abscdfp(2) = cdfplot(abs(dist_err_chstitch320(:)));
    abscdfp(3) = cdfplot(abs(dist_err_chstitch745(:)));
end
set(abscdfp,'linewidth',2)
xlim([0 p_lim])
hold off
grid on
title('')
xlabel('Absolute Error (m)')
ylabel('CDF')
legend(['Per Channel, RMSE: ',num2str(dist_rmse),' m'],...
    ['320 MHz Stitch, RMSE: ',num2str(dist_rmse_chstitch320),' m'],...
    ['745 MHz Stitch, RMSE: ',num2str(dist_rmse_chstitch745),' m'],...
    'location','best','fontsize',12)

p_width = 680;
p_height = 600;
set(1,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])

% Legacy CDFs per bandwidth ----------------------------------------------%
figure(2);
subplot(2,1,1);
hold on
cdfp(1) = cdfplot(dist_err_legacy(:));
if sve.fs == 20e6
    cdfp(2) = cdfplot(dist_err_chstitch320_legacy(:));
    cdfp(3) = cdfplot(dist_err_chstitch745_legacy(:));
end
set(cdfp,'linewidth',2)
hold off
xlim([-p_lim p_lim])
grid on
title([erase(filename,"_"),' Legacy CDFs'])
xlabel('Error (m)')
ylabel('CDF')

subplot(2,1,2);
hold on
abscdfp(1) = cdfplot(abs(dist_err_legacy(:)));
if sve.fs == 20e6
    abscdfp(2) = cdfplot(abs(dist_err_chstitch320_legacy(:)));
    abscdfp(3) = cdfplot(abs(dist_err_chstitch745_legacy(:)));
end
set(abscdfp,'linewidth',2)
xlim([0 p_lim])
hold off
grid on
title('')
xlabel('Absolute Error (m)')
ylabel('CDF')
legend(['Per Channel, RMSE: ',num2str(dist_rmse_legacy),' m'],...
    ['320 MHz Stitch, RMSE: ',num2str(dist_rmse_chstitch320_legacy),' m'],...
    ['745 MHz Stitch, RMSE: ',num2str(dist_rmse_chstitch745_legacy),' m'],...
    'location','best','fontsize',12)

p_width = 680;
p_height = 600;
set(2,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])
