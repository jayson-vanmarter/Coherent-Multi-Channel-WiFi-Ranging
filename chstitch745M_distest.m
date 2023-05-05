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
% This function uses all 20 MHz channel measurements in the 5 GHz band for
% a multi-channel stitch of 745 MHz in bandwidth. AR-model-based
% interpolation is used to interpolate across the large 120 MHz gap between
% the lower and upper 5 GHz bands. MUSIC is utilized for distance
% estimation using the stitched and downsampled CFR.
%
% Inputs:
%   f : Subcarrier channel frequencies corresponding to the two-way CFR
%   measurements
%
%   chanEst_2W : Two-way CFR measurements across channels and antennas
%
%   Delta_f : Frequency measurement spacing
%
%   M : Number of sources (MPCs) specified for subspace decomposition
%
%   c : Propagation speed. This changes the scan window and resolution
%   depending on if we are looking at wireless (3e8) or wired
%   (3e8/sqrt(1.4) measurements
%
%   Dr : Downsampling rate
%
% Outputs:
%   dist_est : Estimated distance
%
%   MUSIC_spec : Computed MUSIC pseudospectrum


function [dist_est,MUSIC_spec] = chstitch745M_distest(f,chanEst_2W,Delta_f,M,c,Dr)

%% Inputs

% Number of antenna measurements
P = size(chanEst_2W,3);

% Index the channel frequencies
bidx = [1,11; % 2.4 GHz channels
    12,21; % 5 GHz lower channels
    22,42]; % 5 GHz upper channels

% Get the subcarrier frequencies of all channels in each band
f_ch_l = f(:,bidx(2,1):bidx(2,2)); % For the 5 GHz lower channels
f_ch_h = f(:,bidx(3,1):bidx(3,2)); % For the 5 GHz upper channels
f_l = f_ch_l(:); % Concatenate to single dimension
f_h = f_ch_h(:); % Concatenate to single dimension
f_l_interp1 = (f_l(1):Delta_f:f_l(end)).'; % Subcarrier frequencies before downsampling for 5 GHz lower
f_h_interp1 = (f_h(1):Delta_f:f_h(end)).'; % Subcarrier frequencies before downsampling for 5 GHz upper
K_B1_l = length(f_l_interp1); % Number of CFR measurements before downsampling for 5 GHz lower
K_B1_h = length(f_h_interp1); % Number of CFR measurements before downsampling for 5 GHz upper

% Get minimum downsampling lengths
K_Dr_l = floor((K_B1_l - Dr)/Dr) + 1; % For 5 GHz lower
K_Dr_h = floor((K_B1_h - Dr)/Dr) + 1; % For 5 GHz upper

% Downsample and truncate frequency indices
f_l_interp = f_l_interp1(1:Dr:end); % Downsample for 5 GHz lower
f_l_interp = f_l_interp(1:K_Dr_l); % Truncate for 5 GHz lower
f_h_interp = f_h_interp1(1:Dr:end); % Downsample for 5 GHz upper
f_h_interp = f_h_interp(1:K_Dr_h); % Truncate for 5 GHz upper
f_5_interp = f_l_interp(1):Dr*Delta_f:f_h_interp(end); % For the full 745 MHz

% Get the missing indices between the lower and upper CFRs to be
% interpolated across
k_l = (f_l_interp - f_l_interp(1)) / (Dr*Delta_f) + 1; % Indices for 5 GHz lower
k_h = (f_h_interp - f_l_interp(1)) / (Dr*Delta_f) + 1; % Indices for 5 GHz upper
k_miss = k_l(end)+1:k_h(1)-1; % Missing measurement indices to be interpolated

% Define the number of measurements after stitching and smoothing subarray
% size
K_B_5 = length(f_5_interp); % Number of measurements for the full 745 MHz
N_5 = round(0.5*K_B_5); % Smoothing subarray size

%% MUSIC struct of parameters for Multi-Band
music.threshold = -15; % MUSIC peak threshold in dB down from the max
music.delta_f = Dr*Delta_f; % Subcarrier spacing after downsampling
music.c = c; % Propagation speed
music.k = 0:K_B_5-1; % CFR measurement indices after stitching
music.N = N_5; % Smoothing subarray size

switch c
    case 3e8
        music.dist_min = 0;
        music.dist_max = 15; % Max distance considered in meters. Larger distance estimations will be shifted to this value. 
        music.delta = 0.02; % MUSIC scan resolution in m
    case 3e8 / sqrt(1.4)
        music.dist_min = -1;
        music.dist_max = 15; % Max distance considered in meters. Larger distance estimations will be shifted to this value. 
        music.delta = 0.001; % MUSIC scan resolution in m
end

music.dists = (music.dist_min:music.delta:music.dist_max); % Corresponding steering matrix distances considered

% Compute steering matrix for the set of possible delays we are considering
music.S = zeros(music.N,length(music.dists));
for m = 1:length(music.dists)
    music.S(:,m) = exp(-1i*4*pi*music.delta_f*(music.dists(m)/music.c)*music.k(1:music.N));
end

%% AR model parameters
order = round(K_B_5/6); % AR model order
K_g = length(k_miss); % Number of subcarriers to predict for extrapolation and subsequent interpolation

%% Interpolate across small gaps between channels
% Index the measurements for the lower and upper 5 GHz bands
H_l_gaps = chanEst_2W(:,bidx(2,1):bidx(2,2),:);
H_h_gaps = chanEst_2W(:,bidx(3,1):bidx(3,2),:);

% Concatenate CFR measurements across multiple channels
H_l_gaps = reshape(H_l_gaps,[],P);
H_h_gaps = reshape(H_h_gaps,[],P);

% Interpolate
H_l = nan(K_B1_l,P);
H_h = nan(K_B1_h,P);
for p = 1:P
    % Interpolate magnitude and phase separately
    H_l_interpmag = interp1(f_l,abs(H_l_gaps(:,p)),f_l_interp1,'linear');
    H_l_interpmphase = interp1(f_l,unwrap(angle(H_l_gaps(:,p))),f_l_interp1,'linear');
    H_h_interpmag = interp1(f_h,abs(H_h_gaps(:,p)),f_h_interp1,'linear');
    H_h_interpmphase = interp1(f_h,unwrap(angle(H_h_gaps(:,p))),f_h_interp1,'linear');

    H_l(:,p) = H_l_interpmag .* exp(1i*H_l_interpmphase);
    H_h(:,p) = H_h_interpmag .* exp(1i*H_h_interpmphase);
end

%% Downsample
H_l_deci = nan(K_Dr_l,P,Dr);
H_h_deci = nan(K_Dr_h,P,Dr);
for i = 1:Dr
    H_l_deci_full = H_l(i:Dr:end,:);
    H_l_deci(:,:,i) = H_l_deci_full(1:K_Dr_l,:);
    H_h_deci_full = H_h(i:Dr:end,:);
    H_h_deci(:,:,i) = H_h_deci_full(1:K_Dr_h,:);
end

%% Interpolate across the large gap bewtween the lower and upper 5 GHz bands
% Use Burg's algorithm to estimate the AR parameters. Then, use the
% estimated AR model to extrapolate across the 120 MHz gap in the 5 GHz
% band
H_5 = nan(K_B_5,P,Dr);
for p = 1:P
    for i = 1:Dr
        H_5(:,p,i) = burgs_interp(H_l_deci(:,p,i),H_h_deci(:,p,i),K_g,order);
    end
end

%% Distance estimation
[MUSIC_spec,dist_est] = music1_est(music,squeeze(H_5),M);