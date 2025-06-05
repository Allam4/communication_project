%% Parameters
fc = 1e3;
num_bits = 256;
samples_per_bit = 100;
bitstream = randi([0 1], 1, num_bits);
N = num_bits * samples_per_bit;
t = linspace(0, num_bits/fc, N);


unipolar_nrz = zeros(1, N);
bipolar_nrz = zeros(1, N);

last_bit = -1;


for i = 1:num_bits
    idx_start = (i-1)*samples_per_bit + 1;
    idx_end = i * samples_per_bit;

    % Unipolar NRZ
    if bitstream(i) == 1
        unipolar_nrz(idx_start:idx_end) = 1;
    else
        unipolar_nrz(idx_start:idx_end) = 0;
    end

    % Bipolar NRZ (AMI)
    if bitstream(i) == 1
        last_bit = -last_bit;
        bipolar_nrz(idx_start:idx_end) = last_bit;
    else
        bipolar_nrz(idx_start:idx_end) = 0;
    end
end


figure(1);
plot(t, unipolar_nrz, 'b'); xlabel('Time (s)'); ylabel('Amplitude');
title('Unipolar NRZ - Time Domain'); grid on;

figure(2);
plot(t, bipolar_nrz, 'r'); xlabel('Time (s)'); ylabel('Amplitude');
title('Bipolar NRZ (AMI) - Time Domain'); grid on;

%% Frequency Domain Setup
ts = t(2) - t(1);
fs = 1 / ts;
df = fs / N;


if rem(N,2) == 0
    f = -fs/2 : df : fs/2 - df;
else
    f = -(fs/2 - df/2) : df : fs/2 - df/2;
end


UNIPOLAR_F = fftshift(fft(unipolar_nrz) / N);
BIPOLAR_F = fftshift(fft(bipolar_nrz) / N);


figure(3);
plot(f, abs(UNIPOLAR_F), 'b'); xlabel('Frequency (Hz)'); ylabel('|Unipolar NRZ(f)|');
title('Unipolar NRZ - Frequency Domain'); grid on;

figure(4);
plot(f, abs(BIPOLAR_F), 'r'); xlabel('Frequency (Hz)'); ylabel('|Bipolar NRZ(f)|');
title('Bipolar NRZ (AMI) - Frequency Domain'); grid on;
