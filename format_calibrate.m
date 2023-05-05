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
% This function formats (reshaping) the data saved from the testbed,
% calibrates, obtains the two-way measurements, and interpolates for null
% subcarrier gaps within each channel.
%
% Inputs:
%   sve : Struct of data for a given packet type and bandwidth
% 
%   sve_cal : Struct of calibraiton data for the same packet type and
%   bandwidth
%
%   c : Propagation speed (3e8 for wireless and 3e8/sqrt(1.4) for wired)
%
% Outputs:
%   chanEst_2W : Formatted and calibrated two-way CFR measurements for the
%   VHT-LTF or HE-LTF.
%   Dimensions:
%       num subcarriers (interpolated) x num channels x num antenna pairs x num time measurements
%       Example: HE packet type 20 MHz: 245 (3x subcarriers interpolated across DC) x 42 x 8 x 2
%       Example: VHT packet type 20 MHz: 57 (1x subcarriers interpolated across DC) x 42 x 8 x 2
% 
%   chanEst_LLTF_2W : Formatted and calibrated two-way CFR measurements for
%   the L-LTF.
%   Dimensions:
%       num subcarriers (interpolated) x num channels x num antenna pairs x num time measurements
%       Example: HE packet type 20 MHz: 53 (1x subcarriers interpolated across DC) x 42 x 8 x 2
%       Example: VHT packet type 20 MHz: 53 (1x subcarriers interpolated across DC) x 42 x 8 x 2
%
%   fInterp : Interpolated channel frequency measurements for the VHT-LTF
%   or HE-LTF.
%   Dimensions: num subcarriers (interpolated) x num channels
%
%   fInterp_LLTF : Interpolated channel frequency measurements for the
%   L-LTF.
%   Dimensions: num subcarriers (interpolated) x num channels
% 
%   SNREst : SNR estimate in dB using the VHT-LTF or HE-LTF.
%   Dimensions: num channels x 2 (reflector/initiator) x num antenna pairs
%
%   SNREst_LLTF : SNR estimate in dB using the L-LTF
%   Dimensions: num channels x 2 (reflector/initiator) x num antenna pairs


function [chanEst_2W,chanEst_LLTF_2W,fInterp,fInterp_LLTF,SNREst,SNREst_LLTF] = ...
format_calibrate(sve,sve_cal,c)

%% Reshape the inputs of interest to separate out antenna pair measurements and leakage measurements
num_ant_pairs = size(sve.ant_pairs,1)-1; % Number of antenna pair measurements
num_perCh = 2*(num_ant_pairs + 1); % Number of measurements per channel. Includes leakage measurements and different measurements on reflector/initiator.

% Separate out the time captures of Legacy and CFRs to the 5th dim for the
% real-world measurements and to the 4th dim for the calibration
% measurements which do not have multiple antennas
td = 2; % Number of time captures per LTF per packet
sve.chanEst = cat(5,sve.chanEst(:,1:end/2),sve.chanEst(:,end/2+1:end));
sve.chanEst_LLTF = cat(5,sve.chanEst_LLTF(:,1:end/2),sve.chanEst_LLTF(:,end/2+1:end));
sve_cal.chanEst = cat(4,sve_cal.chanEst(:,1:end/2),sve_cal.chanEst(:,end/2+1:end));
sve_cal.chanEst_LLTF = cat(4,sve_cal.chanEst_LLTF(:,1:end/2),sve_cal.chanEst_LLTF(:,end/2+1:end));

% Separate the leakage measurements into sve_L.
% SNR Dimensions: channel measurements x reflector/initiator
% CFR Dimensions: subcarriers x channel measurements x reflector/initiator x 1 (single link) x time measurements
sve_L = sve; % Initialize
sve_L.SNREst = [sve.SNREst(2*num_ant_pairs+1:num_perCh:end),sve.SNREst(2*num_ant_pairs+2:num_perCh:end)];
sve_L.SNREst_LLTF = [sve.SNREst_LLTF(2*num_ant_pairs+1:num_perCh:end),sve.SNREst_LLTF(2*num_ant_pairs+2:num_perCh:end)];
sve_L.chanEst = cat(3,sve.chanEst(:,2*num_ant_pairs+1:num_perCh:end,:,:,:),sve.chanEst(:,2*num_ant_pairs+2:num_perCh:end,:,:,:));
sve_L.chanEst_LLTF = cat(3,sve.chanEst_LLTF(:,2*num_ant_pairs+1:num_perCh:end,:,:,:),sve.chanEst_LLTF(:,2*num_ant_pairs+2:num_perCh:end,:,:,:));

% Remove the leakage measurements from the sve and sve_cal structs
sve.SNREst(2*num_ant_pairs+1:num_perCh:end) = [];
sve.SNREst(2*num_ant_pairs+1:num_perCh-1:end) = [];
sve.SNREst_LLTF(2*num_ant_pairs+1:num_perCh:end) = [];
sve.SNREst_LLTF(2*num_ant_pairs+1:num_perCh-1:end) = [];
sve.chanEst(:,2*num_ant_pairs+1:num_perCh:end,:,:,:) = [];
sve.chanEst(:,2*num_ant_pairs+1:num_perCh-1:end,:,:,:) = [];
sve.chanEst_LLTF(:,2*num_ant_pairs+1:num_perCh:end,:,:,:) = [];
sve.chanEst_LLTF(:,2*num_ant_pairs+1:num_perCh-1:end,:,:,:) = [];
sve_cal.SNREst(2+1:4:end) = [];
sve_cal.SNREst(2+1:4-1:end) = [];
sve_cal.SNREst_LLTF(2+1:4:end) = [];
sve_cal.SNREst_LLTF(2+1:4-1:end) = [];
sve_cal.chanEst(:,2+1:4:end,:,:,:) = [];
sve_cal.chanEst(:,2+1:4-1:end,:,:,:) = [];
sve_cal.chanEst_LLTF(:,2+1:4:end,:,:,:) = [];
sve_cal.chanEst_LLTF(:,2+1:4-1:end,:,:,:) = [];

% Seperate reflector and initiator measurements for sve and sve_cal.
% SNR Dimensions: channel measurements x reflector/initiator x antenna pairs
% CFR Dimensions: subcarriers x channel measurements x reflector/initiator x antenna pairs x time measurements
sve.SNREst = [sve.SNREst(1:2:end),sve.SNREst(2:2:end)];
sve.SNREst_LLTF = [sve.SNREst_LLTF(1:2:end),sve.SNREst_LLTF(2:2:end)];
sve.chanEst = cat(3,sve.chanEst(:,1:2:end,:,:,:),sve.chanEst(:,2:2:end,:,:,:));
sve.chanEst_LLTF = cat(3,sve.chanEst_LLTF(:,1:2:end,:,:,:),sve.chanEst_LLTF(:,2:2:end,:,:,:));
sve_cal.SNREst = [sve_cal.SNREst(1:2:end),sve_cal.SNREst(2:2:end)];
sve_cal.SNREst_LLTF = [sve_cal.SNREst_LLTF(1:2:end),sve_cal.SNREst_LLTF(2:2:end)];
sve_cal.chanEst = cat(3,sve_cal.chanEst(:,1:2:end,:,:,:),sve_cal.chanEst(:,2:2:end,:,:,:));
sve_cal.chanEst_LLTF = cat(3,sve_cal.chanEst_LLTF(:,1:2:end,:,:,:),sve_cal.chanEst_LLTF(:,2:2:end,:,:,:));

% Separate out individual antenna pair measurements
SNREst_temp = sve.SNREst;
SNREst_LLTF_temp = sve.SNREst_LLTF;
chanEst_temp = sve.chanEst;
chanEst_LLTF_temp = sve.chanEst_LLTF;
sve.SNREst = [];
sve.SNREst_LLTF = [];
sve.chanEst = [];
sve.chanEst_LLTF = [];
for i = 1:num_ant_pairs
    sve.SNREst = cat(3,sve.SNREst,SNREst_temp(i:num_ant_pairs:end,:));
    sve.SNREst_LLTF = cat(3,sve.SNREst_LLTF,SNREst_LLTF_temp(i:num_ant_pairs:end,:));
    sve.chanEst = cat(4,sve.chanEst,chanEst_temp(:,i:num_ant_pairs:end,:,:,:));
    sve.chanEst_LLTF = cat(4,sve.chanEst_LLTF,chanEst_LLTF_temp(:,i:num_ant_pairs:end,:,:,:));
end

%% Get the SNR Outputs
SNREst = sve.SNREst;
SNREst_LLTF = sve.SNREst_LLTF;

%% Obtain TX and RX gains per channel
channel_txgain_user1 = sve.channel_txgains_user1.';
channel_txgain_user2 = sve.channel_txgains_user2.';
channel_rxgain_user1 = sve.channel_rxgains_user1.';
channel_rxgain_user2 = sve.channel_rxgains_user2.';

channel_cal_txgain_user1 = sve_cal.channel_txgains_user1.';
channel_cal_txgain_user2 = sve_cal.channel_txgains_user2.';
channel_cal_rxgain_user1 = sve_cal.channel_rxgains_user1.';
channel_cal_rxgain_user2 = sve_cal.channel_rxgains_user2.';

%% Obtain parameters, frequencies, and fixed calibration distance

% Number of channels for the packet type and bandwidth
num_ch = length(sve.fc);

% Inputs which depend on the packet type and bandwidth
[delta_f,chN,chIdx,chIdxInterp,p_idx,...
    delta_f_LLTF,chN_LLTF,chIdx_LLTF,chIdxInterp_LLTF,p_idx_LLTF] = get_impparams(sve);

% Frequencies of the channels (VHT-LTF or HE-LTF)
f = nan(length(chIdx),num_ch); % Non-interpolated
fInterp = nan(chN,num_ch); % Interpolated for null-subcarriers within each channel
for i = 1:num_ch
    f(:,i) = chIdx.*delta_f + sve.fc(i);
    fInterp(:,i) = chIdxInterp.*delta_f + sve.fc(i);
end

% Frequencies of the channels (L-LTF)
f_LLTF = nan(length(chIdx_LLTF),num_ch);
fInterp_LLTF = nan(chN_LLTF,num_ch);
for i = 1:num_ch
    f_LLTF(:,i) = chIdx_LLTF.*delta_f_LLTF + sve.fc(i);
    fInterp_LLTF(:,i) = chIdxInterp_LLTF.*delta_f_LLTF + sve.fc(i);
end

% Get fixed calibration distance given the 
switch c
    case 3e8
        cal_dist = 23.34; % For wireless measurements
    case 3e8 / sqrt(1.4)
        cal_dist = 19.5662; % For wired measurements (dielectric_const = 1.4)
    otherwise
        error('No calibration distance for given propagation speed c')
end

%% Compensate for the TX and RX gains per channel

% Leakage HE/VHT
chanEst_L_gcal = nan(size(sve_L.chanEst));
chanEst_L_gcal(:,:,1,:,:) = sve_L.chanEst(:,:,1,:,:) ./ sqrt(10.^((channel_txgain_user1 + channel_rxgain_user1)/10));
chanEst_L_gcal(:,:,2,:,:) = sve_L.chanEst(:,:,2,:,:) ./ sqrt(10.^((channel_txgain_user2 + channel_rxgain_user2)/10));

% Leakage Legacy
chanEst_LLTF_L_gcal = nan(size(sve_L.chanEst_LLTF));
chanEst_LLTF_L_gcal(:,:,1,:,:) = sve_L.chanEst_LLTF(:,:,1,:,:) ./ sqrt(10.^((channel_txgain_user1 + channel_rxgain_user1)/10));
chanEst_LLTF_L_gcal(:,:,2,:,:) = sve_L.chanEst_LLTF(:,:,2,:,:) ./ sqrt(10.^((channel_txgain_user2 + channel_rxgain_user2)/10));

% Calibration HE/VHT
chanEst_cal_gcal = nan(size(sve_cal.chanEst));
chanEst_cal_gcal(:,:,1,:) = sve_cal.chanEst(:,:,1,:) ./ sqrt(10.^((channel_cal_txgain_user1 + channel_cal_rxgain_user1)/10));
chanEst_cal_gcal(:,:,2,:) = sve_cal.chanEst(:,:,2,:) ./ sqrt(10.^((channel_cal_txgain_user2 + channel_cal_rxgain_user2)/10));

% Calibration Legacy
chanEst_LLTF_cal_gcal = nan(size(sve_cal.chanEst_LLTF));
chanEst_LLTF_cal_gcal(:,:,1,:) = sve_cal.chanEst_LLTF(:,:,1,:) ./ sqrt(10.^((channel_cal_txgain_user1 + channel_cal_rxgain_user1)/10));
chanEst_LLTF_cal_gcal(:,:,2,:) = sve_cal.chanEst_LLTF(:,:,2,:) ./ sqrt(10.^((channel_cal_txgain_user2 + channel_cal_rxgain_user2)/10));

% HE/VHT
chanEst_gcal = nan(size(sve.chanEst));
chanEst_gcal(:,:,1,:,:) = sve.chanEst(:,:,1,:,:) ./ sqrt(10.^((channel_txgain_user1 + channel_rxgain_user1)/10));
chanEst_gcal(:,:,2,:,:) = sve.chanEst(:,:,2,:,:) ./ sqrt(10.^((channel_txgain_user2 + channel_rxgain_user2)/10));

% Legacy
chanEst_LLTF_gcal = nan(size(sve.chanEst_LLTF));
chanEst_LLTF_gcal(:,:,1,:,:) = sve.chanEst_LLTF(:,:,1,:,:) ./ sqrt(10.^((channel_txgain_user1 + channel_rxgain_user1)/10));
chanEst_LLTF_gcal(:,:,2,:,:) = sve.chanEst_LLTF(:,:,2,:,:) ./ sqrt(10.^((channel_txgain_user2 + channel_rxgain_user2)/10));

%% Magnitude calibration

% Average the calibration data across time
chanEst_cal_gcal = mean(chanEst_cal_gcal,4);
chanEst_LLTF_cal_gcal = mean(chanEst_LLTF_cal_gcal,4);

% Magnitude calibration for leakage and wireless measurements
chanEst_L_magcal_gcal = chanEst_L_gcal ./ abs(chanEst_cal_gcal);
chanEst_LLTF_L_magcal_gcal = chanEst_LLTF_L_gcal ./ abs(chanEst_LLTF_cal_gcal);
chanEst_magcal_gcal = chanEst_gcal ./ abs(chanEst_cal_gcal);
chanEst_LLTF_magcal_gcal = chanEst_LLTF_gcal ./ abs(chanEst_LLTF_cal_gcal);

%% Phase Calibration (Leakage and Low-Pass Filter)

if sve.fs == 160e6

    % Only need one channel measurement from a single USRP for calibration
    % as the hardware is the same in both USRPs
    chanEst_cal_mean = chanEst_cal_gcal(:,1,1);
    chanEst_LLTF_cal_mean = chanEst_LLTF_cal_gcal(:,1,1);

    % Get linear fit of the phase indices of interest and extrapolate
    chanEst_cal_phase = unwrap(angle(chanEst_cal_mean(p_idx)));
    p = polyfit(chIdx(p_idx),chanEst_cal_phase,1);
    p_linfit = p(1)*chIdx.' + p(2);

    chanEst_LLTF_cal_phase = unwrap(angle(chanEst_LLTF_cal_mean(p_idx_LLTF)));
    p = polyfit(chIdx_LLTF(p_idx_LLTF),chanEst_LLTF_cal_phase,1);
    p_linfit_LLTF = p(1)*chIdx_LLTF.' + p(2);

    % Get the nonlinear phase offset for the phase indices at the edges
    p_offset = unwrap(angle(chanEst_cal_mean)) - p_linfit;
    p_offset_LLTF = unwrap(angle(chanEst_LLTF_cal_mean)) - p_linfit_LLTF;

end

chanEst_gap = nan(size(sve.chanEst));
chanEst_LLTF_gap = nan(size(sve.chanEst_LLTF));
for i = 1:num_ch % Loop for each channel
    for k = 1:2 % Loop for reflector and initiator
        for j = 1:num_ant_pairs % Loop for each antenna pair
            for u = 1:td % Loop for each measurement in time

                % Get the CFRs for the capture
                chanEst_in = chanEst_magcal_gcal(:,i,k,j,u);
                chanEst_L = chanEst_L_magcal_gcal(:,i,k,1,u);

                chanEst_LLTF_in = chanEst_LLTF_magcal_gcal(:,i,k,j,u);
                chanEst_LLTF_L = chanEst_LLTF_L_magcal_gcal(:,i,k,1,u);

                if sve.fs == 160e6

                    % 160 MHz low-pass filter phase calibration
                    chanEst_L = abs(chanEst_L) .* exp(1i*(unwrap(angle(chanEst_L)) - p_offset));
                    chanEst_LLTF_L = abs(chanEst_LLTF_L) .* exp(1i*(unwrap(angle(chanEst_LLTF_L)) - p_offset_LLTF));
                    chanEst_160cal = abs(chanEst_in) .* exp(1i*(unwrap(angle(chanEst_in)) - p_offset));
                    chanEst_LLTF_160cal = abs(chanEst_LLTF_in) .* exp(1i*(unwrap(angle(chanEst_LLTF_in)) - p_offset_LLTF));
        
                    % Leakage phase calibration
                    chanEst_LCal = chanEst_160cal ./ exp(1i*angle(chanEst_L));
                    chanEst_LLTF_LCal = chanEst_LLTF_160cal ./ exp(1i*angle(chanEst_LLTF_L));
            
                else
        
                    % Leakage phase calibration
                    chanEst_LCal = (chanEst_in ./ exp(1i*angle(chanEst_L)));
                    chanEst_LLTF_LCal = (chanEst_LLTF_in ./ exp(1i*angle(chanEst_LLTF_L)));
 
                end

                % Save the interpolated channel estimate
                chanEst_gap(:,i,k,j,u) = chanEst_LCal;
                chanEst_LLTF_gap(:,i,k,j,u) = chanEst_LLTF_LCal;

            end
        end
    end
end

%% Get the two-way CFR measurements after magnitude and phase calibrations
% Note: The integer delay is added to the initiator (:,:,2,:,:) one-way CFR
% phase before saving measurements on the testbed.
chanEst_2W_gap = chanEst_gap(:,:,1,:,:) .* chanEst_gap(:,:,2,:,:);
chanEst_LLTF_2W_gap = chanEst_LLTF_gap(:,:,1,:,:) .* chanEst_LLTF_gap(:,:,2,:,:);

%% Reshape to remove extra dimension
chanEst_2W_gap = reshape(chanEst_2W_gap,length(chIdx),num_ch,num_ant_pairs,td);
chanEst_LLTF_2W_gap = reshape(chanEst_LLTF_2W_gap,length(chIdx_LLTF),num_ch,num_ant_pairs,td);

%% Fixed distance calibration due to leakage delay and hardware delays
chanEst_2W_gap = chanEst_2W_gap .* exp(1i*4*pi*f*(cal_dist/c));
chanEst_LLTF_2W_gap = chanEst_LLTF_2W_gap .* exp(1i*4*pi*f_LLTF*(cal_dist/c));

%% Interpolation for null subcarriers within each channel
chanEst_2W = nan(length(chIdxInterp),num_ch,num_ant_pairs,td);
chanEst_LLTF_2W = nan(length(chIdxInterp_LLTF),num_ch,num_ant_pairs,td);
for i = 1:num_ch % Loop for each channel
    for j = 1:num_ant_pairs % Loop for each antenna pair
        for u = 1:td % Loop for each measurement in time

            % Get the CFRs for the capture
            chanEst_in = chanEst_2W_gap(:,i,j,u);
            chanEst_LLTF_in = chanEst_LLTF_2W_gap(:,i,j,u);

            % Interpolation for magnitude and phase separately
            interpMag = interp1(chIdx,abs(chanEst_in),chIdxInterp,'linear').';
            interpPhase = interp1(chIdx,unwrap(angle(chanEst_in)),chIdxInterp,'linear').';
            chanEst_interp = interpMag .* exp(1i*interpPhase);

            interpMag = interp1(chIdx_LLTF,abs(chanEst_LLTF_in),chIdxInterp_LLTF,'linear').';
            interpPhase = interp1(chIdx_LLTF,unwrap(angle(chanEst_LLTF_in)),chIdxInterp_LLTF,'linear').';
            chanEst_LLTF_interp = interpMag .* exp(1i*interpPhase);

            % Output the interpolated two-way measurements
            chanEst_2W(:,i,j,u) = chanEst_interp;
            chanEst_LLTF_2W(:,i,j,u) = chanEst_LLTF_interp;

        end
    end
end