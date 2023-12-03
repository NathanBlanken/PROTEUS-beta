function pulse = getPulse(f,Ncy,PA,Fs,Tresp,dispFig, T)
% Simulate the pulse generated by a transducer with transfer function T,
% resulting from a three-level driving signal.


    pulse.f = f;                % Centre frequency(Hz)
    pulse.w = 2*pi*pulse.f;     % Angular centre frequency (Hz)

    pulse.Nc = Ncy;             % Number of cycles
    pulse.A = PA;               % Acoustic pressure amplitude (Pa)
    pulse.fs = Fs;              % Sample frequency (Hz)

    
    % Compute bandwidth of transducer (Hz):
    Df = compute_bandwidth(T,Fs);
    
    % Time domain vector (s)
    T_IR = 2/Df;                % Approximate imulse response duration (s)
    T_drive = pulse.Nc/pulse.f; % Drive signal duration (s)
    pulse.t = 0:1/pulse.fs:(T_drive + T_IR + Tresp);

    % Compute the trilevel drive signal
    V = getDriveSignal(Ncy,f,Fs,length(T));
    
    % Compute the acoustic pressure waveform of the incident wave (maximum
    % amplitude normalised to 1):
 
    M       = length(pulse.t);
    [p,dp]  = getPressure(V,T,Fs,M);
    
    % Make sure values at start of signal and end of signal are zero:
    p = removeOffsets(p);

    pulse.p     = p*PA;
    pulse.dp    = dp*PA;
 
    plotWaveform(pulse,dispFig)
    
end

function Df = compute_bandwidth(Tfit,Fs)
% Compute bandwidth of transducer with transfer function Tfit.

% Transfer function on dB scale:
Y_dB = 20*log10(abs(Tfit)/max(abs(Tfit)));
N = length(Y_dB);

% Find the -6 dB bandwidth of the transducer:
I1 = find(Y_dB(1:round(N/2))>-6,1);         % Left  -6 dB point
I2 = find(Y_dB(1:round(N/2))>-6,1,'last');  % Right -6 dB point
Df = (I2 - I1)/N*Fs;                        % Bandwidth (-6 dB) (Hz)

end

function p = removeOffsets(p)
% Make sure values at start of signal and end of signal are zero.

a = p(1);       % Start offset
b = p(end);     % End offset

N = length(p);
n = (0:(N-1));

p = p - a - (b-a)/(N-1)*n;
end

function [p,dp] = getPressure(V,Tfit,Fs,M)
% Compute transmit pressure of a transducer with transfer function Tfit due
% to a drive voltage V.

% Make sure sampling frequency is the same as for Tfit
if Fs ~= 250e6
    error(['Sampling rate of pulse should be the same as sampling'...
        ' rate of transfer function.'])
end

% Apply transfer function Tfit
pfft = Tfit.*fft(V);        % Fourier transform of pressure signal
p = real(ifft(pfft));

% Apply transfer function to compute derivative
N = length(pfft);
f = (0:(N-1))/N*Fs;         % Frequency vector
omega = 2*pi*f;             % Angular frequency vector
omega(ceil(N/2+1):N) = -omega(floor(1+N/2):-1:2);
dpfft = 1i*omega.*pfft;     % Fourier transform of derivative
dp = real(ifft(dpfft));

% Truncate the signal to desired length
p = p(1:M);
dp = dp(1:M);

% Normalise the signal
p0 = max(abs(p));
p = p/p0;
dp = dp/p0;
end

function V = getDriveSignal(Ncy,f,Fs,N)
% Get a pulse train of alternating positive and negative block pulses, with
% Ncy cycles and frequency f at sampling rate F.

ON_Frac = 0.67;             % Fraction of half cycle with high level

NT = Fs/f;                  % Number of sample points per cycle
NT_half = round(Fs/(2*f));  % Number of sample points per half cycle
NT_ON = round(ON_Frac*Fs/(2*f));

V = zeros(1,N);

for k = 1:Ncy
    % Positive pulse
    Nstart = 1 + round((k-1)*NT);
    V(Nstart:Nstart+NT_ON) = 1;
    % Negative pulse
    Nstart = Nstart + NT_half;
    V(Nstart:Nstart+NT_ON) = -1;
end


end

function plotWaveform(pulse,dispFig)

    ti = pulse.t;
    Pacc = pulse.p;
    dPacc = pulse.dp;

    if dispFig == false
        return
    end
    figure(1);
    subplot(1,2,1);
    plot(ti*10^6,Pacc);
    xlabel('$t$ ($\mu$s)','interpreter','latex')
    ylabel('$P_{ac}$ (Pa)','interpreter','latex')
    grid on
    % plot(ti(1:length(input_signal)),abs(input_signal));
    subplot(1,2,2);
    plot(ti*10^6,dPacc);
    xlabel('$t$ ($\mu$s)','interpreter','latex')
    ylabel('d$P$/d$t$ (Pa/s)','interpreter','latex')
    % plot(ti(1:length(input_signal)),abs(input_signal));
    grid on
    %%%%%%%%%%%%%%
end