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
% This function is used to get the parameters for the packet type and
% bandwidth obtained from a loaded file in sve. sve contains a set of
% two-way measurements for a given location, packet type, and bandwidth.
% These measurements span all channels for the given bandwidth and all
% antennas.
%
% Outputs:
%   Delta_f : HE-LTF or VHT-LTF frequency spacing
%
%   K : Number of frequency measurements after interpolation for null
%   subcarriers on the given channel using the HE-LTF or VHT-LTF.
% 
%   chIdx : Frequency measurement indices for the HE-LTF or VHT-LTF.
% 
%   chIdxInterp : Frequency measurement indices after interpolation across
%   null subcarrier gaps for the HE-LTF or VHT-LTF.
%
%   p_idx : Frequency indices to use for linear fitting the 160 MHz phase
%   response. The fitted response is used for mitigation of nonlinear
%   low-pass filter effects in the phase. This is for the HE-LTF or
%   VHT-LTF.
%
%   Delta_f_LLTF : L-LTF frequency spacing.
%
%   K : Number of frequency measurements after interpolation for null
%   subcarriers on the given channel using the L-LTF.
%
%   chIdx_LLTF : Frequency measurement indices for the L-LTF.
% 
%   chIdxInterp_LLTF : Frequency measurement indices after interpolation
%   across null subcarrier gaps for the L-LTF.
%
%   p_idx_LLTF : Frequency indices to use for linear fitting the 160 MHz
%   phase response. The fitted response is used for mitigation of nonlinear
%   low-pass filter effects in the phase. This is for the L-LTF.


function [Delta_f,K,chIdx,chIdxInterp,p_idx,...
    Delta_f_LLTF,K_LLTF,chIdx_LLTF,chIdxInterp_LLTF,p_idx_LLTF] = get_impparams(sve)

if sve.implementation == 0
        
    Delta_f = 78.125e3; % Frequency spacing

    if sve.fs == 20e6
        K = 245; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-122:-2,2:122];
        chIdxInterp = -122:122;
    elseif sve.fs == 40e6
        K = 489; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-244:-3,3:244];
        chIdxInterp = -244:244;
    elseif sve.fs == 80e6
        K = 1001; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-500:-3,3:500];
        chIdxInterp = -500:500;
    elseif sve.fs == 160e6
        K = 2025; % Number of CFR frequencies for MUSIC
        p_idx = (-500:500) + 1013; % Indexes for 160 MHz phase calibration
        chIdx = [-1012:-515,-509:-12,...
            12:509,515:1012];
        chIdxInterp = -1012:1012;
    end
elseif sve.implementation == 2

    Delta_f = 312.5e3; % Frequency spacing

    if sve.fs == 20e6
        K = 57; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-28:-1,1:28];    
        chIdxInterp = -28:28;
    elseif sve.fs == 40e6
        K = 117; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-58:-2,2:58];
        chIdxInterp = -58:58;
    elseif sve.fs == 80e6
        K = 245; % Number of CFR frequencies for MUSIC
        p_idx = nan;
        chIdx = [-122:-2,2:122];
        chIdxInterp = -122:122;
    elseif sve.fs == 160e6
        K = 501; % Number of CFR frequencies for MUSIC
        p_idx = [-126:-6,...
            6:126] + 251; % Indexes for 160 MHz phase calibration without interpolation
        chIdx = [-250:-130,-126:-6,...
            6:126,130:250];
        chIdxInterp = -250:250;
    end
end

% Info for legacy LTF
Delta_f_LLTF = 312.5e3;
if sve.fs == 20e6
    K_LLTF = 53; % Number of CFR frequencies for MUSIC
    p_idx_LLTF = nan;
    chIdx_LLTF = [-26:-1,1:26];    
    chIdxInterp_LLTF = -26:26;
elseif sve.fs == 40e6
    K_LLTF = 117; % Number of CFR frequencies for MUSIC
    p_idx_LLTF = nan;
    chIdx_LLTF = [-58:-33,-31:-6,...
        6:31,33:58];
    chIdxInterp_LLTF = -58:58;
elseif sve.fs == 80e6
    K_LLTF = 245; % Number of CFR frequencies for MUSIC
    p_idx_LLTF = nan;
    chIdx_LLTF = [-122:-97,-95:-70,...
        -58:-33,-31:-6,...
        6:31,33:58,...
        70:95,97:122];
    chIdxInterp_LLTF = -122:122;
elseif sve.fs == 160e6
    K_LLTF = 501; % Number of CFR frequencies for MUSIC
    p_idx_LLTF = [-122:-97,-95:-70,...
        -58:-33,-31:-6,...
        6:31,33:58,...
        70:95,97:122] + 251; % Indexes for 160 MHz phase calibration without interpolation
    chIdx_LLTF = [-250:-225,-223:-198,...
        -186:-161,-159:-134,...
        -122:-97,-95:-70,...
        -58:-33,-31:-6,...
        6:31,33:58,...
        70:95,97:122,...
        134:159,161:186,...
        198:223,225:250];
    chIdxInterp_LLTF = -250:250;
end