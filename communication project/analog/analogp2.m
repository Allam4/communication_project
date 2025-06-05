fs = 100;
df = 0.01;
T = 1/df;
ts = 1/fs;
N = ceil(T/ts);

t = linspace(-T/2, T/2 - ts, N);

% Baseband x(t)
x = zeros(1, N);
x(t >= -4 & t <= 4) = 1;

X_fft = fftshift(fft(x));
if rem(N, 2) == 0
    f = -(0.5*fs) : df : (0.5*fs - df);
else
    f = -(0.5*fs - 0.5*df) : df : (0.5*fs - 0.5*df);
end

BW_lpf = 1;
H_lpf = zeros(1, N);
H_lpf(abs(f) <= BW_lpf) = 1;
Y_fft = X_fft .* H_lpf;
y_filter = real(ifft(ifftshift(Y_fft)));

% Modulate x(t) with fc1
fc1 = 20;
c1 = cos(2 * pi * fc1 * t);
s1 = y_filter .* c1;

figure;
subplot(3,1,1);
plot(t, y_filter);
title("Filtered x(t)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-1 1]);

subplot(3,1,2);
plot(t, c1);
title("Carrier c1(t) (fc = 20 Hz)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-0.2 0.2]);

subplot(3,1,3);
plot(t, s1);
title("DSB-SC Modulated Signal s1(t)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-1 1]);

% Define m(t)
m = zeros(1, N);
condition = (t > 0 & t < 4);
m(condition) = cos(2 * pi * 0.5 * t(condition));

% Hilbert Transform
M_fft = fft(m);
freq_in = -floor(N/2):ceil(N/2)-1;
H_hilbert = -1j * sign(freq_in);
H_hilbert(ceil(N/2)) = 0;  % DC component
Mh_fft = M_fft .* fftshift(H_hilbert);
m_h = real(ifft(Mh_fft));

% SSB Modulation using fc2
fc2 = 23;
c2_cos = cos(2 * pi * fc2 * t);
c2_sin = sin(2 * pi * fc2 * t);
s2 = m .* c2_cos - m_h .* c2_sin;  % USB signal

figure;
subplot(3,1,1);
plot(t, m);
title("m(t)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-1 1]);

subplot(3,1,2);
plot(t, c2_cos);
title("Carrier c₂(t) = cos(2πfc₂t), fc₂ = 23 Hz");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-0.2 0.2]);

subplot(3,1,3);
plot(t, s2);
title("SSB Modulated Signal s2(t)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-1 1]);

% Q11: Combined signal s(t)
s = s1 + s2;
figure;
subplot(2,1,1);
plot(t, s);
title("Combined Modulated Signal s(t) = s1(t) + s2(t)");
xlabel("Time (s)");
ylabel("Amplitude");
grid on;
xlim([-0.5 0.5]);

S_fft = fftshift(fft(s)) * ts;
subplot(2,1,2);
plot(f, abs(S_fft));
title("Magnitude Spectrum |S(f)| of Combined Signal");
xlabel("Frequency (Hz)");
ylabel("|S(f)|");
grid on;
xlim([-30 30]);
hold on;
maxY = max(abs(S_fft(f > 15 & f < 30)));
plot([19 19], [0 maxY], 'r--');
plot([21 21], [0 maxY], 'g--');
plot([23 23], [0 maxY], 'r--');
plot([25 25], [0 maxY], 'g--');
legend("Spectrum", "s1 edge 1", "s1 edge 2", "s2 edge 1", "s2 edge 2");

% Q12: Demodulation
c1_demod = cos(2 * pi * fc1 * t);
demod_s1_mix = s .* c1_demod;
DEMOD_LPF1 = zeros(1, N);
DEMOD_LPF1(abs(f) <= BW_lpf) = 1;
DEMOD_S1_MIX_FFT = fftshift(fft(demod_s1_mix));
R_X_FFT = DEMOD_S1_MIX_FFT .* DEMOD_LPF1;
r_x = real(ifft(ifftshift(R_X_FFT))) * 2;

c2_demod = cos(2 * pi * fc2 * t);
demod_s2_mix = s .* c2_demod;
DEMOD_LPF2 = zeros(1, N);
DEMOD_LPF2(abs(f) <= BW_lpf) = 1;
DEMOD_S2_MIX_FFT = fftshift(fft(demod_s2_mix));
R_M_FFT = DEMOD_S2_MIX_FFT .* DEMOD_LPF2;
r_m = real(ifft(ifftshift(R_M_FFT)));

figure;
subplot(2,1,1);
plot(t, y_filter, 'b');
hold on;
plot(t, r_x, 'r--');
title("Channel 1: Received x(t) vs Original Filtered x(t)");
xlabel("Time (s)");
ylabel("Amplitude");
legend("Original x(t)", "Received x(t)");
grid on;
xlim([-5 5]);

subplot(2,1,2);
plot(t, m, 'b');
hold on;
plot(t, r_m, 'r--');
title("Channel 2: Received m(t) vs Original m(t)");
xlabel("Time (s)");
ylabel("Amplitude");
legend("Original m(t)", "Received m(t)");
grid on;
xlim([-1 5]);

