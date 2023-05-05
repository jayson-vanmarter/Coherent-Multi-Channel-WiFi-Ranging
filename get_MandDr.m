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
% This function is used to return the signal subspace size for MUSIC and
% the downsampling rates given the bandwidth.
%
% Inputs:
%   BW : Bandwidth in MHz used to set outputs
%
% Outputs:
%   M : Number of sources specified for single-channel distance esstimation
%
%   M_320 : Number of sources specified for the 320 MHz multi-channel
%   stitch. Stitching is only done for the 20 MHz single-channel bandwidth
%   case. Output is nan if BW is not 20 MHz.
%
%   M_745 : Number of sources specified for the 745 MHz multi-channel
%   stitch. Stitching is only done for the 20 MHz single-channel bandwidth
%   case. Output is nan if BW is not 20 MHz.
%
%   Dr_VHT : Downsampling rate for the VHT-LTF
%
%   Dr_HE : Downsampling rate for the HE-LTF
%
%   Dr_Legacy : Downsampling rate for the L-LTF
%
%   Dr_VHT_320 : Downsampling rate for the 320 MHz multi-channel stitch
%   using the VHT-LTF. Output is nan if BW is not 20 MHz.
%
%   Dr_HE_320 : Downsampling rate for the 320 MHz multi-channel stitch
%   using the HE-LTF. Output is nan if BW is not 20 MHz.
%
%   Dr_Legacy_320 : Downsampling rate for the 320 MHz multi-channel stitch
%   using the L-LTF. Output is nan if BW is not 20 MHz.
%
%   Dr_VHT_745 : Downsampling rate for the 745 MHz multi-channel stitch
%   using the VHT-LTF. Output is nan if BW is not 20 MHz.
%
%   Dr_HE_745 : Downsampling rate for the 745 MHz multi-channel stitch
%   using the HE-LTF. Output is nan if BW is not 20 MHz.
%
%   Dr_Legacy_745 : Downsampling rate for the 745 MHz multi-channel stitch
%   using the L-LTF. Output is nan if BW is not 20 MHz.

function [M,M_320,M_745,...
    Dr_VHT,Dr_HE,Dr_Legacy,...
    Dr_VHT_320,Dr_HE_320,Dr_Legacy_320,...
    Dr_VHT_745,Dr_HE_745,Dr_Legacy_745] = ...
    get_MandDr(BW)

switch BW
    case 20
        M = 3; % For single-channel case
        M_320 = 40; % For 320 MHz stitch (if using 20 MHz bandwidth)
        M_745 = 60; % For 320 MHz stitch (if using 20 MHz bandwidth)
        Dr_VHT = 2; % For single-channel VHT-LTF case
        Dr_VHT_320 = 7; % For VHT-LTF 320 MHz multi-channel stitch
        Dr_VHT_745 = 10; % For VHT-LTF 745 MHz multi-channel stitch
    case 40
        M = 5;
        M_320 = nan;
        M_745 = nan;
        Dr_VHT = 3;
        Dr_VHT_320 = nan;
        Dr_VHT_745 = nan;
    case 80
        M = 10;
        M_320 = nan;
        M_745 = nan;
        Dr_VHT = 4;
        Dr_VHT_320 = nan;
        Dr_VHT_745 = nan;
    case 160
        M = 20;
        M_320 = nan;
        M_745 = nan;
        Dr_VHT = 5;
        Dr_VHT_320 = nan;
        Dr_VHT_745 = nan;
end

Dr_HE = 4*Dr_VHT; % For single-channel HE-LTF case
Dr_Legacy = Dr_VHT; % For L-LTF case
Dr_HE_320 = 4*Dr_VHT_320; % For HE-LTF 320 MHz multi-channel stitch
Dr_Legacy_320 = Dr_VHT_320; % For L-LTF 320 MHz multi-channel stitch
Dr_HE_745 = 4*Dr_VHT_745; % For HE-LTF 320 MHz multi-channel stitch
Dr_Legacy_745 = Dr_VHT_745; % For L-LTF 320 MHz multi-channel stitch