%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is a supplementary material for the paper:                  %
% G.Firtha and P.Fiala, Theory and Implementation of 2.5D WFS of moving   %
%    sources with arbitrary trajectory, in proceedings of DAGA2018        %
%                                                                         %
%  The script generates driving functions for a circular SSD from an      %
%  arbitrary input sound file, being the excitation signal of a moving    %
%  point source on an arbitrary trajectory                                %
%  Outputs:                                                               %
%   d_wfs: i.-th column contains the i.-th loudspeaker's driving signal   %
%   output = 'wav' is set: i channel out.wav is created                   %
%   output = 'SSR' is set: i number of wave and an asdf file are generated%
% (c) 2018 by Gergely Firtha                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
addpath('Files')
c = 343.1;
%% User defined parameters
% Virtual source properties
v = 15;                                       % Virtual source velocity [m/s]
[input,fs] = audioread('Samples/fire.mp3');   % Source excitation signal
output = 'SSR';                               % wav / SSR
% if wav: saves output signals as a multichannel wav file
% if SSR: saves i number of wav files and generates an asdf file for the SSR

% Anchor points for source trajectory
x_a = [  -3   -3  -3  -2.5   -1.5    0  100;  %x_coordinates [m]
    -100  -1   0   1.5    2.5    3   3 ]; %y_coordinates [m]

% SSD properties
N_ssd =   50;                                 % Number of loudspeakers
R_ssd  =   2;                                 % Radius of circular SSD
R_ref  =   0.1;                               % Radius of reference circle
AA_filter = 'on';                             % Antialiasing filtering enabled (on/off)
AA_type   = 'freq_domain';                    % AA filter type: freq_domain/time_domain

subsamp = 200;                                % subsampling parameter of source trajectory
%% Create circular SSD
fi = (0:2*pi/N_ssd:2*pi-2*pi/N_ssd)';
x0 =  [ cos(fi)  sin(fi) ]*R_ssd;
n0 = -[ cos(fi)  sin(fi) ];
v0 =  [ sin(fi) -cos(fi) ];
% Get actual driving function lengths
[ p,xp,yp ] = make_path( x_a(1,:), x_a(2,:), 150 );
T_sim = p(end)/v;
Nt = floor(T_sim*fs/subsamp)*subsamp;
t = (0:Nt-1)'/fs;
ts = t(1:subsamp:end);
input = sum(input(1:Nt),1)';
% Get source trajectory as the function of time t
xs = get_trajectory( p,xp,yp, ts, v );

% Plot trajectory
f = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
p1 = plot(xs(:,1),xs(:,2));
hold on
draw_ssd( p1, x0(1:1:end,:), n0(1:1:end,:), 0.04 );
axis equal tight
grid on
xlim( [-R_ssd-2,R_ssd+2] );
ylim( [-R_ssd-2,R_ssd+2] );
plot(x_a(1,:),x_a(2,:),'ok')
title('SSD geometry and source trajectory')
xlabel('x -> [m]')
ylabel('y -> [m]')
drawnow
clear xp yp p
% Calculate initial source position/propagation time delay for each SSD element at t = 0
Tau0  = get_initial_position( v,c, x_a, x0 );
% Get amplitudes and delays at time instants ts
[ A, Tau, wc ] = get_amps_and_taus( ts, x0,n0,v0, xs, Tau0,c, R_ref );
% Filter input signal with ideal WFS prefilters and apply amplitude and delays
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)'/(Nt)*fs );
s_wfs = ifft( sqrt(1i*w/(c*2*pi)).*fft(input) );
%%
d_wfs = zeros(subsamp*size(Tau,1),length(x0));
wb = waitbar(0,'Calculating driving functions');
for n = 1 : length(x0)
    waitbar(n/length(x0),wb);
    d_wfs(:,n) = interp1( t,input, t-interp1(ts,Tau(:,n),t), 'linear','extrap' ).*interp1(ts,A(:,n),t, 'linear','extrap');
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
subplot(1,2,2)
plot(t,sum(d_wfs,2))
xlim([t(1),t(end)]);
grid on
xlabel('t -> [s]')
ylabel('Sum( s(t))');
title('Synthesized field at the center of SSD')
%%
if strcmp(output,'wav')
    audiowrite('out.wav',d_wfs,fs);
elseif strcmp(output,'SSR')
    mkdir Audio
    for n = 1 : size(x0,1)
        audiowrite(sprintf('Audio/out_%03i.wav',n),d_wfs(:,n),fs);
    end
    
    fileID = fopen(sprintf('WFS_%d_ls.asd',N_ssd),'w');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n<asdf>\n  <header>\n    <name>WFS system simulation</name>\n');
    fprintf(fileID,'    <description></description>\n  </header>\n\n');
    fprintf(fileID,'  <scene_setup>\n\n');
    for n = 1 : N_ssd
        fprintf(fileID,'<source name="ls%03i" model="point"><file>Audio/out_%03i.wav</file><position x="%f" y="%f"/><orientation azimuth="%2d"/></source>\n',...
            n,n, x0(n,1),x0(n,2),round(atand(-n0(n,1)./n0(n,2))));
    end   
    fprintf(fileID,'\n  </scene_setup>\n');
    fprintf(fileID,'</asdf>');
    fclose(fileID);
else
end
