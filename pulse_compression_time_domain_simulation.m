% This script was written by Mark D. Hickle of Missouri University of
% Science and Technology, March 2023


clearvars;
clc;

j=1i;

create_gif = 1;%0 = no, 1 = yes
gif_name = 'lfm_mf_3targets_20p0us_1e6_500m_2.gif';
gif_delay = 1/30; %Frame period
iter = 1;

waveform_selection = 1; %0 = simple pulse, 1 = LFM, 2 = phase-coded
windowing = 0;          %0 = no windowing, 1 = hamming, 2 = blackman
apply_matched_filter = 1;
fs = 20e6;              % sampling frequency
dt = 1/fs;              % sampling period

rmax = 40000;            % km
tmax = rmax/3e8;
tau = 20e-6;
t_wave = 0:dt:tau;
t = -tau:dt:(tmax-dt); % time vector
range = t*3e8;
range_rec = range/2;
npoints = length(t);
npoints_pulse = length(t_wave)-1;
plot_period = 5; %Generates gif frame every plot_period time steps

f0 = 1e9;
target_ranges = [10000 10500 11000];
target_amplitudes = [1, 1, 1];

targets = target_ranges/3e8;
target_phases = ones(1,length(targets));
targets_idx = round((targets+tau)/dt);



B = 1e6;

disp(['This simulation will generate ',num2str(floor(npoints/plot_period)),' gif frames'])


waveform = zeros(1,npoints_pulse);
switch waveform_selection
    case 0
        waveform = rectpuls(t_wave-tau/2,tau);
    case 1
        waveform = 1i*(rectpuls(t_wave-tau/2,tau).*exp(1i*2*pi*B/tau.*(t_wave).^2));
    case 2
        % length-5 biphase barker
%         phase_sequence = [1 1 1 -1 1]; 

        % length-13 biphase barker
%         phase_sequence = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; 

        %28-length polyphase barker        
        phase_sequence = [0 0 4.553 4.087 2.215 1.894 3.946 3.771 .0922 .968 2.931 1.003 2.452 3.461 1.330 3.745 1.917 4.431 2.297 4.417]; phase_sequence = exp(1i*phase_sequence);  

        nchips = length(phase_sequence);
        npoints_chip = ceil(npoints_pulse/nchips);
        t_chip = (0:(npoints_chip-1))*dt;
        for k=1:nchips
            waveform((1+(k-1)*npoints_chip):(k*npoints_chip))=phase_sequence(k).*rectpuls(t_chip-tau/nchips/2,tau/nchips);
        end
end
npoints_pulse = length(waveform);

mf = conj(waveform(npoints_pulse:-1:1))/(tau*fs);

if((waveform_selection==1)&&(windowing == 1))
    window = hamming(length(mf))';
    mf = mf.*window/mean(window);
end
if((waveform_selection==1)&&(windowing == 2))
    window = blackman(length(mf))';
    mf = mf.*window/mean(window);
end



%Initialize signals before entering for loop
x_ref = zeros(1,npoints);
x_t = zeros(1,npoints);
x_rec = NaN*zeros(1,npoints);
x_rec_filt = NaN*zeros(1,npoints);
x_rec_filt_temp = zeros(1,npoints);

for ts = 1:1:npoints
    if(ts<npoints_pulse)
        x_t(1) = waveform(ts);
    else
        x_t(1)=0;
    end

    x_t(2:(npoints)) = x_t(1:(npoints-1)); %at each time step, shift Tx signal samples to the right by 1 (forward propagating wave)
    x_ref(targets_idx) = x_ref(targets_idx) + x_t(targets_idx).*target_amplitudes.*target_phases;
    x_ref(1:(npoints-1)) = x_ref(2:npoints); 
    x_ref_filt = fftfilt(mf,x_ref);
    x_rec(ts) = x_ref(npoints_pulse);
    x_rec_filt_temp(ts) = x_ref_filt(npoints_pulse)+0.000000001*randn();
    if(ts>npoints_pulse)
        x_rec_filt(ts-npoints_pulse)=x_rec_filt_temp(ts);
    end
    if((mod(ts,plot_period)==0)&&(ts>=npoints_pulse))
        figure(1);
        subplot(2,1,1);
        plot(range/1e3,real(x_t),'k','LineWidth',1); %Transmitted signal
        hold on;
        plot(range/1e3,real(x_ref),'r','LineWidth',1);%Reflected signal
        
        %Plot dashed target lines
        plot([target_ranges(1)/1e3, target_ranges(1)/1e3+.01],[0,1]*abs(target_amplitudes(1)),'k--','LineWidth',1);
        plot([target_ranges(2)/1e3, target_ranges(2)/1e3+.01],[0,1]*abs(target_amplitudes(2)),'k--','LineWidth',1);
        plot([target_ranges(3)/1e3, target_ranges(3)/1e3+.01],[0,1]*abs(target_amplitudes(3)),'k--','LineWidth',1);
        hold off;
        ylim([-3 3])
        xlim([0 30]);
        xlabel('Actual range (km)');
        ylabel({'Signal Amplitude'})
        legend('Real part of transmitted signal','Real aprt of reflected signal','Targets');
        title('Transmitted & reflected signals');
        grid on;

        subplot(2,1,2);
        
        if(apply_matched_filter)
            plot((range_rec+tau)/1e3,(abs(x_rec_filt)),'r','LineWidth',1);
        else
            plot((range_rec+tau)/1e3,(abs(x_rec)),'r','LineWidth',1);
        end
        xlim([0 15]);
        xticks(0:1:15);
        ylim([0 1.1]);
        yticks(0:.2:1);
        xlabel('Measured range (km)');
        ylabel('Signal Amplitude')
        title('Received signal after matched filtering');
        grid on;
        set(gcf,'Color','w')
        drawnow;
        if(create_gif)
            if(iter==1)
                gif(gif_name,'DelayTime',gif_delay);
            else
                gif;
            end
        end
        iter = iter+1;
    end
end

