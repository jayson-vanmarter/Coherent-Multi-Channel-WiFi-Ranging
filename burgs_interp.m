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
% This function estimates AR model parameters using Burg's method for the
% lower and upper 5 GHz bands. Then, the AR model is used to extrapolate
% the CFR across the CFR gap. The extrapolated CFRs using the lower and
% upper CFR AR model parameters is linearly cross-faded for the
% interpolated CFR estimate.
%
% Inputs:
%   r_l : CFR for 5 GHz lower
%
%   r_h : CFR for 5 GHz upper
%
%   K_g : Size of the CFR gap in number of subcarriers
%
%   order : AR model order
%
% Output:
%   r_comb : Stitched CFR output of 745 MHz BW

function r_comb = burgs_interp(H_l,H_h,K_g,order)

K_l = length(H_l); % Length of the lower set of measurements
K_h = length(H_h); % Length of the upper set of measurements

% Linearly decaying averaging filters
w_f = linspace(1,0,K_g).'; 
w_b = linspace(0,1,K_g).';

% Use Burg's algorithm to estimate the AR model parameters
a_51 = arburg(H_l,order);
a_52 = arburg(H_h,order);

% Forward prediction
H_f = zeros(1,K_l+K_g);
H_f(1:K_l) = H_l;
for j = K_l+1:K_l+K_g
    H_f(j) = -sum(a_51(2:end).*H_f((j-1):-1:(j-order)));
end
H_f_exp = H_f(K_l+1:K_l+K_g).';

% Backward prediction
H_b = zeros(1,K_h+K_g);
H_b(K_g+1:end) = H_h;
for j = K_g:-1:1
    H_b(j) = -sum(conj(a_52(2:end)).*H_b((j+1):(j+order)));
end
H_b_exp = H_b(1:K_g).';

% Average and stitch
H_exp = w_f.*H_f_exp + w_b.*H_b_exp; % Average by linear cross-fading
r_comb = [H_l;H_exp;H_h]; % Stitch