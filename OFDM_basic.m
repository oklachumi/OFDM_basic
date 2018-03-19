% Simulate the effect of ISI as the length of a guard interval(CP, CS, or ZP) 
% It considers the BER performance of an OFDM system with 64-point FFT and 
% 16 virtual carriers for 16-QAM signaling in the AWGN or multipath 
% Rayleigh fading channel(with the maximum delay of 15 samples).
clear,close,clc all

function y = guard_interval(Ng,Nfft,NgType,ofdmSym)
if NgType == 1 % CP
    y = [ofdmSym(Nfft-Ng+1:Nfft) ofdmSym(1:Nfft)];
elseif NgType == 2 % ZP
    y = [zeros(1,Ng) ofdmSym(1:Nfft)];
endif
endfunction

function y = remove_GI(Ng,Lsym,NgType,ofdmSym)
if Ng ~= 0
    if NgType == 1 % CP  
        y = ofdmSym(Ng+1:Lsym);
    elseif NgType == 2 % CS
        y = ofdmSym(1:Lsym-Ng) + [ofdmSym(Lsym-Ng+1:Lsym) zeros(1,Lsym-2*Ng)];
    endif
else
    y = ofdmSym;
endif
endfunction

function y = remove_CP(x,Ncp,Noff)
if nargin < 3
    Noff = 0;
endif
y = x(:,Ncp+1-Noff:end-Noff);
endfunction

function y = Q(x)
y = erfc(x/sqrt(2))/2;
endfunction

function plot_ber(file_name,Nbps)
EbN0dB = [0:1:30];
M = 2^Nbps;
ber_AWGN = ber_QAM(EbN0dB,M,'AWGN');
ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
semilogy(EbN0dB,ber_AWGN,'r:');
hold on;
semilogy(EbN0dB,ber_Rayleigh,'r-');
a = load(file_name);
semilogy(a(:,1),a(:,2),'b--s');
grid on
legend('AWGN analytic','Rayleigh fading analytic','Simulation');
xlabel('EbN0[dB]'); ylabel('BER'); axis([a(1,1) a(end,1) 1e-5 1])
endfunction

function ber = ber_QAM(EbN0dB,M,AWGN_or_Rayleigh)
N = length(EbN0dB);
sqM = sqrt(M); 
a = 2*(1-power(sqM,-1))/log2(sqM);
b = 6*log2(sqM)/(M-1);
if nargin < 3
    AWGN_or_Rayleigh = 'AWGN';
endif
if lower(AWGN_or_Rayleigh(1)) == 'a'
    ber = a*Q(sqrt(b*10.^(EbN0dB/10)));
else 
    rn = b*10.^(EbN0dB/10)/2;
    ber = 0.5*a*(1-sqrt(rn./(rn+1)));
endif
endfunction


NgType = 1; % NgType=1/2 for cyclic prefix/zero padding
if NgType == 1
    nt = 'CP';
elseif NgType == 2
    nt = 'ZP';
end
Ch=1; % Ch=0/1 for AWGN/multipath channel
if Ch == 0
    chType = 'AWGN';
    Target_neb = 100;
else 
    chType = 'CH';
    Target_neb = 500;
end
%figure(Ch+1), clf
PowerdB = [0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
Delay = [0 3 5 6 8];          % Channel delay 'sample'
Power = 10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Ntap = length(PowerdB);       % Chanel tap number
Lch = Delay(end)+1;           % Channel length
Nbps = 4;
M = 2^Nbps; % Modulation order=2/4/6 for QPSK/16QAM/64QAM
Nfft = 64;  % FFT size
Ng = 3;%Nfft/4;  % Ng=0: Guard interval length
% Ng=Nfft/4;
Nsym = Nfft + Ng;      % Symbol duration
Nvc = Nfft/4;        % Nvc=0: no virtual carrier 16
Nused = Nfft-Nvc;    % 48 = 64 - 16
 
EbN0 = [0:3:30];    % EbN0
N_iter = 1e2;       % Number of iterations for each EbN0
Nframe = 3;         % Number of symbols per frame
sigPow = 0;         % Signal power initialization
file_name = ['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat']; % OFDM_BER_AWGN_CP_GL16
fid = fopen(file_name, 'w+');
norms = [1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM
for i = 0:length(EbN0)
    randn('state',0); rand('state',0); %rand狀態都一樣
    Ber = 0; % BER initialization
    Neb = 0; Ntb = 0; % Initialize the number of error/total bits
    for m = 1:N_iter
        % Tx
        X = randint(1,Nused*Nframe,M); % bit: integer vector 1 x 48*3, M=16
        Xmod = qammod(X,M)/norms(Nbps); % qammod(X,M,0,'gray')/norms(Nbps);
        if NgType ~= 2
            x_GI = zeros(1,Nframe*Nsym); % 3*80
        elseif NgType == 2
            x_GI = zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros
        end
        kk1 = [1:Nused/2]; % Nused = 48
        kk2 = [Nused/2+1:Nused];
        kk3 = 1:Nfft; % 64
        kk4 = 1:Nsym; % 80
        for k = 1:Nframe
            if Nvc ~= 0 % 16
                X_shift = [0 Xmod(kk2) zeros(1,Nvc-1) Xmod(kk1)];
            else
                X_shift = [Xmod(kk2) Xmod(kk1)];
            end
            x = ifft(X_shift);
            x_GI(kk4) = guard_interval(Ng,Nfft,NgType,x);
            kk1 = kk1 + Nused;
            kk2 = kk2 + Nused;
            kk3 = kk3 + Nfft;
            kk4 = kk4 + Nsym;
        end
        if Ch == 0 % AWGN
            y=x_GI; % No channel
        else % Multipath fading channel
            channel = (randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);
            h = zeros(1,Lch);
            h(Delay+1) = channel; % cir: channel impulse response
            y = conv(x_GI,h); 
        endif
        if i == 0 % Only to measure the signal power for adding AWGN noise
            y1 = y(1:Nframe*Nsym);
            sigPow = sigPow + y1*y1';
            continue;
        endif
        % Add AWGN noise
        snr = EbN0(i) + 10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
        noise_mag = sqrt((10.^(-snr/10))*sigPow/2); % N0=Eb/SNR
        y_GI = y + noise_mag*(randn(size(y)) + j*randn(size(y)));
        % Rx
        kk1 = (NgType == 2)*Ng + [1:Nsym];
        kk2 = 1:Nfft;
        kk3 = 1:Nused;
        kk4 = Nused/2 + Nvc + 1:Nfft;
        kk5 = (Nvc~=0) + [1:Nused/2];
        if Ch == 1
            H = fft([h zeros(1,Nfft-Lch)]); % Channel frequency response
            H_shift(kk3) = [H(kk4) H(kk5)]; 
        end
        for k = 1:Nframe
            Y(kk2) = fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
            Y_shift = [Y(kk4) Y(kk5)];
            if Ch == 0
                Xmod_r(kk3) = Y_shift;
            else 
                Xmod_r(kk3) = Y_shift./H_shift;  % Equalizer - channel compensation
            end
            kk1 = kk1 + Nsym;
            kk2 = kk2 + Nfft;
            kk3 = kk3 + Nused;
            kk4 = kk4 + Nfft;
            kk5 = kk5 + Nfft;
        end
        X_r = qamdemod(Xmod_r*norms(Nbps),M);
        Neb = Neb + sum(sum(de2bi(X_r,Nbps)~=de2bi(X,Nbps)));
        Ntb = Ntb + Nused*Nframe*Nbps;  %[Ber,Neb,Ntb]=ber(bit_Rx,bit,Nbps); 
        if Neb > Target_neb
            break;
        end
    end
    if i == 0
        sigPow = sigPow/Nsym/Nframe/N_iter;
        fprintf('Signal power = %11.3e\n', sigPow);
        fprintf(fid,'%%Signal power= %11.3e\n%%EbN0[dB]       BER\n', sigPow);
    else
        Ber = Neb/Ntb;     
        fprintf('EbN0=%3d[dB], BER=%4d/%8d =%11.3e\n', EbN0(i), Neb, Ntb, Ber)
        fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
        if Ber < 1e-6
            break;
        end
    end
end

if (fid ~= 0)
    fclose(fid);
endif
disp('Simulation is finished');
plot_ber(file_name,Nbps);

