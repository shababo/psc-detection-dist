function [filtered_trace, event_times, event_sizes] = wiener_filter(trace,template,ar_noise_params,nfft,dt,threshold,min_window)
%
% ex = wienerFilter(y,h,sigma,gamma,alpha);
%
% Generalized Wiener filter using parameter alpha. When
% alpha = 1, it is the Wiener filter. It is also called
% Regularized inverse filter.
%
% Reference: Richb's paper
% Created: Tue May 4 16:24:06 CDT 1999, Huipin Zhang

Fs = 1/dt;

N = length(trace);
trace_f = fft(trace,nfft); % replace with multitaper
template_f = fft(template,nfft);
[trace_P, freq_mt] = pmtm(trace,9,nfft*2,Fs); 
% min(freq_mt)
% max(freq_mt)
% size(trace_P)
% trace_P_period = abs(trace_f).^2;
trace_P = trace_P(1:end-1)';
trace_P = trace_P/sum(trace_P);

% get parameterized ar noise psd
[noise_P, freqs_noise_P] = freqz(sqrt(ar_noise_params.sigma_sq),ar_noise_params.phi,nfft,Fs);

% min(freqs_noise_P)
% max(freqs_noise_P)

% direct implementation of the regularized inverse filter, 
% when alpha = 1, it is the Wiener filter
% Gf = conj(template_f).*Pxf./(abs(template_f.^2).*Pxf+alpha*noise_P);
%
% Since we don't know Pxf, the following 
% handle singular case (zero case)
% template_f_clean = template_f.*(abs(template_f)>0)+1/gamma*(abs(template_f)==0);
template_f_clean = template_f;
inverse_template_f = 1./template_f_clean;
% inverse_template_f = inverse_template_f.*(abs(template_f)*gamma>1)+gamma*abs(template_f_clean).*inverse_template_f.*(abs(template_f_clean)*gamma<=1);



noise_P = noise_P'/sum(noise_P);
% size(trace_P)
% size(noise_P)
first_i = find(trace_P<=noise_P,1,'first');
% trace_P = trace_P.*(trace_P>noise_P)+noise_P.*(trace_P<=noise_P);
trace_P = trace_P.*((1:nfft)<first_i)+noise_P.*((1:nfft)>=first_i);
wien_filter = inverse_template_f.*(trace_P-noise_P)./(trace_P); %in denom: -(1-alpha)*noise_P


% loglog(freqs_noise_P,abs(trace_P),'b'); hold on;
% % loglog(freqs_noise_P,abs(trace_P_period),'g'); hold on;
% loglog(freqs_noise_P,abs(noise_P),'r'); hold on;



% max(max(abs(Gf).^2)) % should be equal to gamma^2
% Restorated image without denoising
filtered_trace = wien_filter.*trace_f;
filtered_trace = real(ifft2(filtered_trace));
filtered_trace = filtered_trace(1:length(trace));


[~, event_times] = findpeaks(filtered_trace,'MinPeakHeight',threshold*std(filtered_trace),'MinPeakDistance',min_window);
event_sizes = zeros(size(event_times));
for j = 1:length(event_times)
    baseline = min(trace(max(1,event_times(j)-40):event_times(j)));
    event_sizes(j) = max(trace(event_times(j):min(event_times(j)+200,length(trace)))) - 10 - baseline;
end

event_times(event_sizes < 5) = [];
event_sizes(event_sizes < 5) = [];
