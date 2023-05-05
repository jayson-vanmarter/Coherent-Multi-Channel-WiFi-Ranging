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
% This script uses saved data from process_data_trials.m to plot results
% across different bandwidths. The saved data provided is from evaluating
% all measurements taken in LOS conditions.

clear
close all

%% Load in the data

% Single-channel results
dist_err_HESU20 = load('.\savedata\HESU_20_AllLOS').dist_err;
dist_err_HESU40 = load('.\savedata\HESU_40_AllLOS').dist_err;
dist_err_HESU80 = load('.\savedata\HESU_80_AllLOS').dist_err;
dist_err_HESU160 = load('.\savedata\HESU_160_AllLOS').dist_err;
dist_err_VHT20 = load('.\savedata\VHT_20_AllLOS').dist_err;
dist_err_VHT40 = load('.\savedata\VHT_40_AllLOS').dist_err;
dist_err_VHT80 = load('.\savedata\VHT_80_AllLOS').dist_err;
dist_err_VHT160 = load('.\savedata\VHT_160_AllLOS').dist_err;
dist_err_Legacy20 = cat(3,load('.\savedata\HESU_20_AllLOS').dist_err_legacy,load('.\savedata\VHT_20_AllLOS').dist_err_legacy);
dist_err_Legacy40 = cat(3,load('.\savedata\HESU_40_AllLOS').dist_err_legacy,load('.\savedata\VHT_40_AllLOS').dist_err_legacy);
dist_err_Legacy80 = cat(3,load('.\savedata\HESU_80_AllLOS').dist_err_legacy,load('.\savedata\VHT_80_AllLOS').dist_err_legacy);
dist_err_Legacy160 = cat(3,load('.\savedata\HESU_160_AllLOS').dist_err_legacy,load('.\savedata\VHT_160_AllLOS').dist_err_legacy);

% 320 MHz Multi-channel stitch
dist_err_chstitch320_HESU = load('.\savedata\HESU_20_AllLOS').dist_err_chstitch320;
dist_err_chstitch320_VHT = load('.\savedata\VHT_20_AllLOS').dist_err_chstitch320;
dist_err_chstitch320_Legacy = cat(3,load('.\savedata\HESU_20_AllLOS').dist_err_chstitch320_legacy,load('.\savedata\VHT_20_AllLOS').dist_err_chstitch320_legacy);

% 745 MHz Multi-channel stitch
dist_err_chstitch745_HESU = load('.\savedata\HESU_20_AllLOS').dist_err_chstitch745;
dist_err_chstitch745_VHT = load('.\savedata\VHT_20_AllLOS').dist_err_chstitch745;
dist_err_chstitch745_Legacy = cat(3,load('.\savedata\HESU_20_AllLOS').dist_err_chstitch745_legacy,load('.\savedata\VHT_20_AllLOS').dist_err_chstitch745_legacy);

%% Concatenate results across LTF types

dist_err_20 = cat(3,dist_err_HESU20,dist_err_VHT20,dist_err_Legacy20);
dist_err_40 = cat(3,dist_err_HESU40,dist_err_VHT40,dist_err_Legacy40);
dist_err_80 = cat(3,dist_err_HESU80,dist_err_VHT80,dist_err_Legacy80);
dist_err_160 = cat(3,dist_err_HESU160,dist_err_VHT160,dist_err_Legacy160);
dist_err_chstitch320 = cat(3,dist_err_chstitch320_HESU,dist_err_chstitch320_VHT,dist_err_chstitch320_Legacy);
dist_err_chstitch745 = cat(3,dist_err_chstitch745_HESU,dist_err_chstitch745_VHT,dist_err_chstitch745_Legacy);

%% Compute RMSEs

dist_rmse_20 = sqrt(mean(dist_err_20(:).^2));
dist_rmse_40 = sqrt(mean(dist_err_40(:).^2));
dist_rmse_80 = sqrt(mean(dist_err_80(:).^2));
dist_rmse_160 = sqrt(mean(dist_err_160(:).^2));
dist_rmse_chstitch320 = sqrt(mean(dist_err_chstitch320(:).^2));
dist_rmse_chstitch745 = sqrt(mean(dist_err_chstitch745(:).^2));

%% Plot

p_lim = 2;

try close 1
catch end
figure(1);
hold on
cdfp(1) = cdfplot(abs(dist_err_20(:)));
cdfp(2) = cdfplot(abs(dist_err_40(:)));
cdfp(3) = cdfplot(abs(dist_err_80(:)));
cdfp(4) = cdfplot(abs(dist_err_160(:)));
cdfp(5) = cdfplot(abs(dist_err_chstitch320(:)));
cdfp(6) = cdfplot(abs(dist_err_chstitch745(:)));
set(cdfp,'linewidth',2.5)
xlim([0 p_lim])
grid on
grid minor
title('')
xlabel('Absolute Error (m)')
ylabel('CDF')
legend(['20 MHz, RMSE: ',num2str(round(dist_rmse_20,3)),' m'],...
    ['40 MHz, RMSE: ',num2str(round(dist_rmse_40,3)),' m'],...
    ['80 MHz, RMSE: ',num2str(round(dist_rmse_80,3)),' m'],...
    ['160 MHz, RMSE: ',num2str(round(dist_rmse_160,3)),' m'],...
    ['325 MHz Stitch, RMSE: ',num2str(round(dist_rmse_chstitch320,3)),' m'],...
    ['745 MHz Stitch, RMSE: ',num2str(round(dist_rmse_chstitch745,3)),' m'],...
    'location','best','fontsize',12)

p_width = 500;
p_height = 400;
set(1,'Position',[1920/2-p_width/2,1080/2-p_height/2,p_width,p_height])