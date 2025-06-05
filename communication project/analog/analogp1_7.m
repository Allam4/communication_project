% Q1
t_plot = linspace(-5, 5, 1000);
x_t = zeros(size(t_plot));
idx1 = find(t_plot >= -4 & t_plot <= 0);
x_t(idx1) = t_plot(idx1) + 5;
idx2 = find(t_plot > 0 & t_plot <= 4);
x_t(idx2) = -t_plot(idx2) + 5;
figure(1);
plot(t_plot, x_t);
hold on;
plot([-5 5], [0 0], 'k--');
plot([0 0], [0 5.5], 'k--');
grid on;
title('Function x(t)');
xlabel('Time (s)');
ylabel('Amplitude');
axis([-5 5 -0.5 5.5]);
legend('x(t)');
set(gca, 'FontSize', 12);
drawnow;

% Q3
fs = 100;
df = 0.01;
T = 1/df;
ts = 1/fs;
N = ceil(T/ts);
t = linspace(-T/2, T/2 - ts, N);
xr = zeros(size(t));
for i = 1:length(t)
    if t(i) >= -4 && t(i) <= 4
        xr(i) = 5 - abs(t(i));
    else
        xr(i) = 0;
    end
end
NU_Xf = fftshift(fft(xr)) * ts;
if (rem(N,2)==0)
    f = -(0.5*fs):df:(0.5*fs-df);
else
    f = -(0.5*fs-0.5*df):df:(0.5*fs-0.5*df);
end

function y = safe_sinc(x)
    y = ones(size(x));
    in= (x ~= 0);
    y(in) = sin(pi*x(in)) ./ (pi*x(in));
end

A_Xf = 8 * safe_sinc(8*f) + 16 * (safe_sinc(4*f)).^2;

figure;
plot(f, abs(NU_Xf));
hold on;
plot(f, abs(A_Xf));
title("Comparison of Numerical FFT and Analytical FT for x(t)");
xlabel("Frequency (f) [Hz]");
ylabel("|X(f)|");
legend("Num_FFT", "An_FT");
grid on;
xlim([-5 5]);

% Q4
Power = abs(NU_Xf).^2;
max_p = max(Power);
threshold_p = 0.05 * max_p;
p = find(Power >= threshold_p);
minfreq_in = min(p);
maxfreq_in = max(p);
freqmin = f(minfreq_in);
freqmax = f(maxfreq_in);
pfreq_in = find(f >= 0);
power_p = Power(pfreq_in);
freq_p = f(pfreq_in);
bp_in = find(power_p < threshold_p);
if isempty(bp_in)
    f_bw_positive = freq_p(end);
    warning("Power did not drop below 5%% threshold within the calculated frequency range.");
else
    first_below_index = bp_in(1);
    f_bw_positive = freq_p(first_below_index);
end
BW = 2 * f_bw_positive;
result_str = sprintf("Estimated 5%% Power Bandwidth (BW) for rev x(t): %.4f Hz\n", BW);
disp(result_str);
fileID = fopen("bw_xr.txt", "w");
fprintf(fileID, "%s", result_str);
fclose(fileID);
disp("Bandwidth estimation saved to bw_xr.txt");

figure;
plot(f, 10*log10(Power / max_p + eps));
hold on;
plot(f, 10*log10(0.05) * ones(size(f)), "r--");
plot([-f_bw_positive, f_bw_positive], [10*log10(0.05), 10*log10(0.05)], "g-x", "LineWidth", 2, "MarkerSize", 10);
title("Normalized PSD of Revised x(t) and 5% Bandwidth Estimation");
xlabel("Frequency (f) [Hz]");
ylabel("Normalized PSD (dB)");
legend("Normalized PSD", "5% Threshold (-13 dB)", "BW Limits");
grid on;
xlim([-2 2]);
ylim([-40 5]);

% Q5
Xf_input = fftshift(fft(xr));
BW_lpf = 1;
Hf_lpf = zeros(size(f));
Hf_lpf(abs(f) <= BW_lpf) = 1;
Xf_output = Xf_input .* Hf_lpf;
yfilter = real(ifft(ifftshift(Xf_output)));

figure;
plot(t, xr, "b-", "LineWidth", 1.5);
hold on;
plot(t, yfilter, "r--", "LineWidth", 1.5);
title("Input vs. Output of Perfect LPF (BW=1Hz) for Revised x(t)");
xlabel("Time (t) [s]");
ylabel("Amplitude");
legend("Input x(t)", "Output y(t)");
grid on;
xlim([-10 10]);

% Q6
Xf_input = fftshift(fft(xr));
BW_lpf = 0.3;
Hf_lpf = zeros(size(f));
Hf_lpf(abs(f) <= BW_lpf) = 1;
Xf_output = Xf_input .* Hf_lpf;
yfilter = real(ifft(ifftshift(Xf_output)));

figure;
plot(t, xr, "b-", "LineWidth", 1.5);
hold on;
plot(t, yfilter, "r--", "LineWidth", 1.5);
title("Input vs. Output of Perfect LPF (BW=0.3Hz) for Revised x(t)");
xlabel("Time (t) [s]");
ylabel("Amplitude");
legend("Input x(t)", "Output y(t)");
grid on;
xlim([-10 10]);

% Q7-1
ts = 0.01;
t = -1:ts:5;
m_r= zeros(size(t));
for i = 1:length(t)
    if t(i) > 0 && t(i) < 4
        m_r(i) = cos(pi * t(i));
    else
        m_r(i) = 0;
    end
end
figure;
plot(t, m_r);
title("Plot of m(t) = cos(pi*t) for 0 < t < 4");
xlabel("Time (t)");
ylabel("m(t)");
grid on;
axis([-1 5 -1.2 1.2]);

% Q7-3
fs = 100;
df = 0.01;
T = 1/df;
ts = 1/fs;
N = ceil(T/ts);
t = linspace(-T/2, T/2 - ts, N);
m_r= zeros(size(t));
for i = 1:length(t)
    if t(i) > 0 && t(i) < 4
        m_r(i) = cos(pi * t(i));
    else
        m_r(i) = 0;
    end
end
Nu_Mf = fftshift(fft(m_r)) * ts;
if (rem(N,2)==0)
   f=-(0.5*fs) : df : (0.5*fs-df);
else
   f=-(0.5*fs-0.5*df) : df : (0.5*fs-0.5*df);
end
function y = safe_sinc(x)
    y = ones(size(x));
    in = (x ~= 0);
    y(in) = sin(pi*x(in)) ./ (pi*x(in));
end
An_Mf = 2 .* exp(-1j*4*pi*f) .* (safe_sinc(4.*f - 2) + safe_sinc(4.*f + 2));
figure;
plot(f, abs(Nu_Mf));
hold on;
plot(f, abs(An_Mf));
title("Comparison of Numerical FFT and Analytical FT for m(t)");
xlabel("Frequency (f) [Hz]");
ylabel("|M(f)|");
legend("Numerical FFT", "Analytical FT");
grid on;
xlim([-5 5]);

% Q7-4
PowerSD = abs(Nu_Mf).^2;
max_po = max(PowerSD);
thr_power = 0.05 * max_po;
abv_thr_in = find(PowerSD >= thr_power);
minfreq_in = min(abv_thr_in);
maxfreq_in = max(abv_thr_in);
freq_min = f(minfreq_in);
freq_max = f(maxfreq_in);
BW= freq_max - freq_min;
peak_freq = 0.5;
[~, peak_in] = min(abs(f - peak_freq));
in_upper = peak_in;
while in_upper <= N && PowerSD(in_upper) >= thr_power
    in_upper = in_upper + 1;
end
if in_upper > N
    f_upper = f(N);
else
    f_upper = f(in_upper);
end
in_lower = peak_in;
while in_lower >= 1 && PowerSD(in_lower) >=thr_power
    in_lower = in_lower - 1;
end
if in_lower < 1
    f_lower = f(1);
else
    f_lower = f(in_lower);
end
BW_total = (f_upper - f_lower) * 2;
result_str = sprintf("Estimated 5%% Power Bandwidth (BW) for m(t): %.4f Hz\n", BW);
disp(result_str);
fileID = fopen("bw_m_r.txt", "w");
fprintf(fileID, "%s", result_str);
fclose(fileID);
disp("Bandwidth estimation saved to bw_m_r.txt");

figure;
plot(f, 10*log10(PowerSD/max_po + eps));
hold on;
plot(f, 10*log10(0.05) * ones(size(f)), "r--");
plot([freq_min, freq_max], [10*log10(0.05), 10*log10(0.05)], "g-x","MarkerSize", 10);
title("Normalized PowerSD of m(t) and 5% Bandwidth");
xlabel("Frequency (f) [Hz]");
ylabel("Normalized PowerSD (dB)");
legend("Normalized PowerSD", "5% Threshold (-13 dB)", "BW Limits");
grid on;
xlim([-2 2]);
ylim([-40 5]);


