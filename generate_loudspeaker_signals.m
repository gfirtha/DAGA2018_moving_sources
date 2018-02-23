%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is a supplementary material for the paper:                  %
% G.Firtha and P.Fiala, Theory and Implementation of 2.5D WFS of moving   %
%    sources with arbitrary trajectory, in proceedings of DAGA2018        %
%                                                                         %
%  The script generates driving functions for a circular SSD from an      %
%  arbitrary input sound file, being the excitation signal of a moving    %
%  point source on an arbitrary trajectory                                %
%  Outputs:                                                               %
%    d_wfs: i.-th column contains the i.-th loudspeaker's driving signal  %
%    if output = 'wav' is set: out.wav is created                         %
%                                                                         %
% (c) 2018 by Gergely Firtha                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
addpath('Files')
c = 343.1;
%% User defined parameters
% Virtual source properties
v = 20;                                       % Virtual source velocity [m/s]
[in,fs] = audioread('Samples/fire.mp3');      % Source excitation signal

% Anchor points for source trajectory
x_a = [  -3   -3  -3  -2.5   -1.5    0  100;  %x_coordinates [m]
        -100  -1   0   1.5    2.5    3   3 ]; %y_coordinates [m]

% SSD properties
N_ssd =   50;                                 % Number of loudspeakers
Rssd  =   2;                                  % Radius of circular SSD
Rref  =   0.1;                                % Radius of reference circle
AA_filter = 'on';                             % Antialiasing filtering enabled (on/off)
AA_type   = 'freq_domain';                    % AA filter type: freq_domain/time_domain

output = 'wav';                               % saving output signals as a multichannel wav file
subsamp = 200;                                % subsampling parameter of source trajectory
%% Create circular SSD
fi = (0:2*pi/N_ssd:2*pi-2*pi/N_ssd)';
x0 =  [ cos(fi)  sin(fi) ]*Rssd;
n0 = -[ cos(fi)  sin(fi) ];
v0 =  [ sin(fi) -cos(fi) ];
% Get actual driving function lengths
[ p,xp,yp ] = make_path( x_a(1,:), x_a(2,:), 150 );
Tsim = p(end)/v;
Nt = floor(Tsim*fs/subsamp)*subsamp;
t = (0:Nt-1)'/fs;
ts = t(1:subsamp:end);
in = sum(in(1:Nt),1)';
% Get source trajectory as the function of time t
xs = get_trajectory( p,xp,yp, ts, v );
f = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
p1 = plot(xs(:,1),xs(:,2));
hold on
draw_ssd( p1, x0(1:1:end,:), n0(1:1:end,:), 0.03 );
axis equal tight
grid on
xlim( [-Rssd-2,Rssd+2] );
ylim( [-Rssd-2,Rssd+2] );
drawnow
clear xp yp p
% Calculate initial source position/propagation time delay for each SSD element at t = 0
Tau0  = get_initial_position( v,c, x_a, x0 );
% Get amplitudes and delays at time instants ts
[ A, Tau, wc ] = get_amps_and_taus( ts, x0,n0,v0, xs, Tau0,c, Rref );
% Filter input signal with ideal WFS prefilters and apply amplitude and delays
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)'/(Nt)*fs );
s_wfs = ifft( sqrt(1i*w/(c*2*pi)).*fft(in) );

d_wfs = zeros(subsamp*size(Tau,1),length(x0));
wb = waitbar(0,'Calculating driving functions');
for n = 1 : length(x0)
    waitbar(n/length(x0),wb);
    d_wfs(:,n) = interp1( t,in, t-interp1(ts,Tau(:,n),t), 'linear','extrap' ).*interp1(ts,A(:,n),t, 'linear','extrap');
end
close(wb);
clear A Tau w
%%
if strcmp(AA_filter,'on')
    if strcmp(AA_type,'freq_domain')
        % Frequency domain antialiasing filtering
        wlen = 2048;
        hop = wlen/4;
        nfft = wlen;
        d_wfs  = anti_aliasing_fd(fs,x0,d_wfs,ts,wc,wlen,hop,nfft);
    elseif strcmp(AA_type,'time_domain')
        % Time domain antialiasing filtering
        Nfil = 200;
        step = 4096;
        d_wfs = anti_aliasing_td( d_wfs, t,ts, wc, Nfil, step );
    end
elseif strcmp(AA_filter,'off')
end
d_wfs = real(d_wfs);
%%
if strcmp(output,'wav')
    audiowrite('out.wav',d_wfs,fs);
end

subplot(1,2,2)
plot(t,sum(d_wfs,2))
xlim([t(1),t(end)]);
grid on
xlabel('t -> [s]')
ylabel('Sum( s(t))');
title('Synthesized field at the center of SSD')
