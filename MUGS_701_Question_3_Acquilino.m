% McGill University - Schulich School of Music - Comprehensive Exam - MUGS 701
% Matlab code Question 3
% Student: Alberto Acquilino

clear all
close all
clc

filename = 'C:\Users\Alberto\Desktop\McGill\Comp exam\Q3\sounds\Good1.wav';
[y_stereo, Fs] = audioread(filename);           % import audio file
A4 = 442;                                       % reference pitch (Hz)

y1 = y_stereo(:, 1);                            % consider only one channel
N = length(y1);                                 % length of the signal
t = (0 : N-1) / Fs;                             % time vector (in seconds)

y  = detrend(y1, 'constant');                   % removes DC component
y_abs = abs(y);                                 % calculate the absolute value of the audio signal

% ENVELOPE TRACKING
n_env = 250;                                    % window length for the linear peak envelope algorithm
local_maxima = zeros(1, floor(N/n_env));        % linear peak envelope vector initialization
time_ref = zeros(1, floor(N/n_env));            % linear envelope time vector initialization

for i = n_env : n_env : N
  time_ref(i/n_env) = i / Fs;
  local_maxima(i/n_env) = max(y_abs(i-n_env+1:i)); % linear peak envelope
end

RMS_env = envelope(y,500,'rms');                % root-mean-square envelope
[up,lo] = envelope(y_abs,500,'peak');           % spline peak envelope


% Plot options
lwidth = 1.5;

% Text options
fsize = 12;
tsize = 16;
font  = 'Segoe UI';

% Figure 1
figure(1)
plot(t, y_abs, 'b', 'linewidth', .5)
hold on
plot(t, RMS_env, 'g', 'linestyle', '--', 'linewidth', lwidth)
plot(t, up, 'r', 'linestyle', '--', 'linewidth', lwidth)
plot(time_ref, local_maxima, 'c', 'linestyle', '--', 'linewidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[0 0 1 1])
set(gca, 'fontname', font, 'fontsize', fsize);
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Normalized amplitude', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
title('Amplitude envelope methods', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
legend('Audio signal', 'RMS envelope', 'Spline peak envelope', 'Linear peak envelope')
box on;

% The considered trumpet range is among E3 = 164.8 Hz to Bb5 = 932 Hz
f_min = 150;    % minimum frequency considered for the pitch function
f_max = 1000;   % maximum frequency considered for the pitch function

% Given the interest in the transient, a short window with high overlap is
% chosen
windowLength  = round(Fs * 3 / f_min);
overlapLength = round(windowLength * .9);


[f0,idx]     = pitch(y, Fs, 'WindowLength', windowLength, 'OverlapLength', overlapLength, 'Range', [f_min,f_max], 'Method', 'NCF');
[f0_SRH,idx] = pitch(y, Fs, 'WindowLength', windowLength, 'OverlapLength', overlapLength, 'Range', [f_min,f_max], 'Method', 'SRH');
[f0_CEP,idx] = pitch(y, Fs, 'WindowLength', windowLength, 'OverlapLength', overlapLength, 'Range', [f_min,f_max], 'Method', 'CEP');
time_f0 = idx/Fs;

% Notes frequency calculation
n = -17 : 13;
f_notes = A4 * (2^(1/12)).^n;


% Figure 2
figure(2)
subplot(4, 1, 1)
window = hamming(512); % window with a size of 512 points
noverlap = 256; % number of points for repeating the window
nfft = 1024; % size of the fft
[S, F, T, P] = spectrogram(y, window, noverlap, nfft, Fs, 'yaxis');
surf(T, F, 10*log10(P), 'edgecolor', 'none'); axis tight; view(0,90);
title('Spectrogram', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);

subplot(4, 1, 2)
plot(time_f0, 12*log2(f0/A4), 'b', 'linewidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[0 0 1 1])
a = get(gca,'YTickLabel');
set(gca, 'fontname', font, 'fontsize', fsize);
set(gca,'YTickLabel', a, 'FontName', 'Times', 'fontsize', 5, 'xlim', [time_f0(1), time_f0(end)], 'ylim', [min(n), max(12*log2(f0/A4))])
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize)
yticks(n)
yticklabels({'', '', '', 'G_3', '', '', '', '', 'C_4', '', '', '', '', '', '', 'G_4', '', '', '', '', 'C_5', '', '', '', '', '', '', 'G_5', '', '', ''})
% To have a detailed note grid on the y axis comment the previous line and uncomment the following one
% yticklabels({'E_3', 'F_3', 'F#_3', 'G_3', 'Ab_3', 'A_3', 'Bb_3', 'B_3', 'C_4', 'C#_4', 'D_4', 'Eb_4', 'E_4', 'F_4', 'F#_4', 'G_4', 'Ab_4', 'A_4', 'Bb_4', 'B_4', 'C_5', 'C#_5', 'D_5', 'Eb_5', 'E_5', 'F_5', 'F#_5', 'G_5', 'Ab_5', 'A_5', 'Bb_5'})
title('Fundamental frequency - NCF method', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
box on;

subplot(4, 1, 3)
plot(time_f0, 12*log2(f0_SRH/A4), 'b', 'linewidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[0 0 1 1])
a = get(gca,'YTickLabel');
set(gca, 'fontname', font, 'fontsize', fsize);
set(gca,'YTickLabel', a, 'FontName', 'Times', 'fontsize', 5, 'xlim', [time_f0(1), time_f0(end)])
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize)
yticks(n)
yticklabels({'', '', '', 'G_3', '', '', '', '', 'C_4', '', '', '', '', '', '', 'G_4', '', '', '', '', 'C_5', '', '', '', '', '', '', 'G_5', '', '', ''})
% To have a detailed note grid on the y axis comment the previous line and uncomment the following one
% yticklabels({'E_3', 'F_3', 'F#_3', 'G_3', 'Ab_3', 'A_3', 'Bb_3', 'B_3', 'C_4', 'C#_4', 'D_4', 'Eb_4', 'E_4', 'F_4', 'F#_4', 'G_4', 'Ab_4', 'A_4', 'Bb_4', 'B_4', 'C_5', 'C#_5', 'D_5', 'Eb_5', 'E_5', 'F_5', 'F#_5', 'G_5', 'Ab_5', 'A_5', 'Bb_5'})
title('Fundamental frequency - SRH method', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
box on;

subplot(4, 1, 4)
plot(time_f0, 12*log2(f0_CEP/A4), 'b', 'linewidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[0 0 1 1])
a = get(gca,'YTickLabel');
set(gca, 'fontname', font, 'fontsize', fsize);
set(gca,'YTickLabel', a, 'FontName', 'Times', 'fontsize', 5, 'xlim', [time_f0(1), time_f0(end)])
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize)
yticks(n)
yticklabels({'', '', '', 'G_3', '', '', '', '', 'C_4', '', '', '', '', '', '', 'G_4', '', '', '', '', 'C_5', '', '', '', '', '', '', 'G_5', '', '', ''})
% To have a detailed note grid on the y axis comment the previous line and uncomment the following one
% yticklabels({'E_3', 'F_3', 'F#_3', 'G_3', 'Ab_3', 'A_3', 'Bb_3', 'B_3', 'C_4', 'C#_4', 'D_4', 'Eb_4', 'E_4', 'F_4', 'F#_4', 'G_4', 'Ab_4', 'A_4', 'Bb_4', 'B_4', 'C_5', 'C#_5', 'D_5', 'Eb_5', 'E_5', 'F_5', 'F#_5', 'G_5', 'Ab_5', 'A_5', 'Bb_5'})
title('Fundamental frequency - CEP method', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
box on;

% Figure 3
figure(3)
subplot(3, 1, 1)
window = hamming(512); % window with a size of 512 points
noverlap = 256; % number of points for repeating the window
nfft = 1024; % size of the fft
[S, F, T, P] = spectrogram(y, window, noverlap, nfft, Fs, 'yaxis');
surf(T, F, 10*log10(P), 'edgecolor', 'none'); axis tight; view(0,90);
title('Spectrogram', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);

subplot(3, 1, 2)
plot(time_f0, 12*log2(f0/A4), 'b', 'linewidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[0 0 1 1])
a = get(gca,'YTickLabel');
set(gca, 'fontname', font, 'fontsize', fsize);
set(gca,'YTickLabel', a, 'FontName', 'Times', 'fontsize', 5, 'xlim', [time_f0(1), time_f0(end)], 'ylim', [min(n), max(12*log2(f0/A4))])
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize)
yticks(n)
yticklabels({'', '', '', 'G_3', '', '', '', '', 'C_4', '', '', '', '', '', '', 'G_4', '', '', '', '', 'C_5', '', '', '', '', '', '', 'G_5', '', '', ''})
% To have a detailed note grid on the y axis comment the previous line and uncomment the following one
% yticklabels({'E_3', 'F_3', 'F#_3', 'G_3', 'Ab_3', 'A_3', 'Bb_3', 'B_3', 'C_4', 'C#_4', 'D_4', 'Eb_4', 'E_4', 'F_4', 'F#_4', 'G_4', 'Ab_4', 'A_4', 'Bb_4', 'B_4', 'C_5', 'C#_5', 'D_5', 'Eb_5', 'E_5', 'F_5', 'F#_5', 'G_5', 'Ab_5', 'A_5', 'Bb_5'})
title('Fundamental frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
box on;

% ATTACK ESTIMATION/DETECTION ALGORITHM

amp_max = max(y_abs);          % maximum amplitude of the signal
thr_att = .15 * amp_max;       % threshold to define the start attack time is a percentage of the maximum value of the ampitude envelope
thr_rel = .25 * amp_max;       % threshold to define the release time is a percentage of the maximum value of the amplitude envelope

subplot(3, 1, 3)
plot(time_ref, local_maxima, 'b', 'linestyle', '-', 'linewidth', lwidth)
set(gca, 'xlim', [time_ref(1), time_ref(end)])
hold on 
plot(time_ref, thr_att * ones(size(time_ref)), 'r', 'linestyle', '--', 'linewidth', lwidth)
plot(time_ref, thr_rel * ones(size(time_ref)), 'g', 'linestyle', '--', 'linewidth', lwidth)
xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
ylabel('Normalized amplitude', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
title('Amplitude envelope', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
legend('Amplitude envelope', 'Attack threshold', 'Release threshold')

%  START OF THE ATTACK DEFINITION
% I define the start attack time as the first sample in which the amplitude
% reaches the defined threshold
t_att_start = 1;

while y_abs(t_att_start) < thr_att
  t_att_start = t_att_start + 1;
end

% END OF THE ATTACK DEFINITION
% I firstly set the end of the attack equal to the start of the attack on
% the same time resolution of the pitch function
t_att_end = 1;
while idx(t_att_end) < t_att_start
  t_att_end = t_att_end + 1;
end
% t_att_end = t_att_end - 1;

% I then define the end of the attack as the point in which the fundamental
% frequency remains within a neighborhood of its value for at least 0.3 seconds

t_min = Fs * .3;  % duration in seconds of the pitch duration stability multiplied by the sampling rate
f0_thres = 1;     % frequency wide band within which a stable pitch is expected. 1 = 100 cents = 1 semitone

t_min_idx = 1;
while idx(t_min_idx) < t_min
  t_min_idx = t_min_idx + 1;
end

while max(12*log2(f0(t_att_end:t_att_end+t_min_idx)/A4)) - min(12*log2(f0(t_att_end:t_att_end+t_min_idx)/A4)) > f0_thres
  t_att_end = t_att_end + 1;
end

t_att_end = idx(t_att_end) / Fs; % end attack in seconds
t_att_start = t_att_start / Fs;  % start attack in seconds

% Attack duration in milliseconds
att_dur = (t_att_end - t_att_start) * 10^3;

% Array of obtained results (hard-coded)
Index_attack_results = 1 : 12;
Label_attack_results = {'Good 1', 'Good 2', 'Good 3', 'Good 4', 'Good 5', 'Good 6', 'Bad 1', 'Bad 2', 'Bad 3', 'Bad 4', 'Bad 5', 'Bad 6'};

Results    = [6.6042, 3.6875, 3.0208, 9.9167, 1.1667, 1.5833, 89.2708, 138.1042, 304.4375, 169.7917, 57.2292, 206.5417]; % attack time measures [ms]
LAT_TT     = [-2.0846, -0.7435, -1.0003, -1.4181, -0.9854, -0.8840, -0.2239, -1.0282, -1.4082, -1.1751, -0.8518, -1.9231]; % log-attack-time results from the Timbre Toolbox
Results_TT = 10.^LAT_TT * 10^3; % Timbre Toolbox attack time measures [ms]

% Figure 4
figure(4)
scatter(Index_attack_results, Results_TT, 'r', 'filled','d')
hold on
scatter(Index_attack_results, Results, 'MarkerEdgeColor', [0 .5 .5], 'MarkerFaceColor', [0 .7 .7], 'LineWidth', lwidth)
% Figure options
set(gcf, 'color', 'w', 'units','normalized','outerposition',[.2 .2 .7 .7])
a = get(gca,'XTickLabel');
set(gca, 'fontname', font, 'fontsize', fsize, 'xlim', [0, length(Index_attack_results) + 1]);
set(gca,'XTickLabel', a, 'FontName', 'Times', 'fontsize', fsize)
ylabel('Attack duration [ms]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize)
xticks(Index_attack_results)
xticklabels(Label_attack_results)
xtickangle(45)
title('Attack duration results', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);
legend('Weakest-effort method', 'Proposed method')
box on;


% RELEASE ESTIMATION/DETECTION ALGORITHM
% The release time is firstly imposed as the time in which the maximum
% of the envelope is reached
[max_env, idx_rel] = max(local_maxima);

% The release time is then defined as the first sample in which the amplitude
% envelope reaches the specified threshold
while local_maxima(idx_rel) > thr_rel
  idx_rel = idx_rel + 1;
end
idx_rel = idx_rel - 1;

t_rel = time_ref(idx_rel); % release time in seconds

% f0_der1 = gradient(f0,idx(2)-idx(1));
% 
% figure
% plot(idx/Fs,f0_der1, 'b', 'linewidth', lwidth)
% xlabel('time [s]', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
% ylabel('f_{0} I^{st} derivative', 'fontweight', 'bold', 'fontname', font, 'fontsize', fsize);
% set(gca, 'xlim', [time_f0(1), time_f0(end)])
% title('First derivative of fundamental frequency', 'fontweight', 'bold', 'fontname', font, 'fontsize', tsize);

% centroid   = spectralCentroid(y, Fs);
% t_centroid = linspace(0, N/Fs, size(centroid, 1));
% 
% plot(t_centroid, centroid)
% set(gca, 'xlim', [t_centroid(1), t_centroid(end)])
% 
% windowWidth = 1;
% B = 1/windowWidth*ones(windowWidth, 1);
% out = filter(B,1,f0);
% close all
