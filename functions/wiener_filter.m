function [filtered_trace, event_times] = wiener_filter(trace,template,ar_noise_params,gamma,nfft,dt,threshold,min_window)
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
trace_P = abs(trace_f).^2;

% get parameterized ar noise psd
[noise_P, freqs_noise_P] = freqz(ar_noise_params.sigma_sq,ar_noise_params.phi,nfft,Fs);

% direct implementation of the regularized inverse filter, 
% when alpha = 1, it is the Wiener filter
% Gf = conj(template_f).*Pxf./(abs(template_f.^2).*Pxf+alpha*noise_P);
%
% Since we don't know Pxf, the following 
% handle singular case (zero case)
template_f_clean = template_f.*(abs(template_f)>0)+1/gamma*(abs(template_f)==0);
inverse_template_f = 1./template_f_clean;
inverse_template_f = inverse_template_f.*(abs(template_f)*gamma>1)+gamma*abs(template_f_clean).*inverse_template_f.*(abs(template_f_clean)*gamma<=1);

% THIS NEEDS ELEMENTWISE LOGIC OPERATIONS... OR AT LEAST CHECKING THAT IT
% IS

noise_P = noise_P';
trace_P = trace_P.*(trace_P>noise_P)+noise_P.*(trace_P<=noise_P);
wien_filter = inverse_template_f.*(trace_P-noise_P)./(trace_P-noise_P); %(1-alpha)*noise_P

figure;
loglog(noise_P,'r'); hold on;
loglog(trace_P,'b'); hold on;


% max(max(abs(Gf).^2)) % should be equal to gamma^2
% Restorated image without denoising
filtered_trace = wien_filter.*trace_f;
filtered_trace = real(ifft2(filtered_trace));

event_times = findpeaks(filtered_trace,'MinPeakHeight',threshold,'MinPeakDistance',min_window);

