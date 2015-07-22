function event_onsets = find_pscs(trace, dt, tau, amp_thesh, conv_thresh, low_passed, plot_figs)

% tau = .001
t_vector = (0:length(trace)-1)*dt;

% convolution = zeros(length(trace),1);

noise_sd = std(trace - smooth(trace,100,'sgolay',4)');

% mean(trace - smooth(trace,50,'sgolay',4)')

prev_samples = floor(3*tau/dt);
template = alpha_synapse(t_vector,0,tau,-1);
template = [zeros(1,prev_samples) template];

if low_passed
    template = template - smooth(template,500,'sgolay',2)';
end

template(floor(13*tau/dt):end) = [];

template = template *8;

% template = fliplr(template);

% figure;
% plot(trace)
% hold on
% plot(950:(950 + length(template) -1),template-29)




colored_noise = conv(normrnd(0,noise_sd,size(trace)),template,'same');



convolution = conv([trace zeros(1,length(template))],template,'same');
% convolution(1:ceil(length(template)/2)-1) = [];
if length(convolution) > length(trace)
    convolution(length(trace)+1:end) = [];
end

window_size = floor(2*tau/dt);
if mod(window_size,2) == 0
    window_size = window_size + 1;
end


thresh_vec = smooth([(trace - median(trace)) zeros(1,window_size)],window_size)';
% thresh_vec(1:ceil(window_size/2)-1) = [];
thresh_vec(1:ceil(tau/dt)-1) = [];
if length(thresh_vec) > length(trace)
    thresh_vec(length(trace)+1:end) = [];
end

% figure; plot(trace); hold on; plot(thresh_vec); title('test')
thresh_vec = thresh_vec < -amp_thesh;



% if plot_figs
%     figure;
%     plot(t_vector,trace,'b',t_vector,convolution/max(convolution)*10,'g',t_vector,thresh_vec*10,'r')
% 
% 
% 
%     figure; 
%     plot(trace,'g')
%     hold on
%     plot(trace - smooth(trace,100,'sgolay',4)','b')
%     hold on
%     plot(smooth(trace,100,'sgolay',4),'r')
% end
% thresh_vec = zeros(size(trace));
% 
% for t_idx = 1:length(trace)
%     
%     t = t_vector(t_idx);
%     
% %     template_max = min(-5,-abs(mean(trace(t_idx:min(length(trace),t_idx+floor(2*tau/dt))))));
%     
%     thresh_vec(t_idx) = mean(trace(t_idx:min(length(trace),t_idx+floor(tau/dt)))) < -2.5*noise_sd;
%     
%     template = alpha_synapse(t_vector,t,tau,-1);
%     
%     
%     convolution(t_idx) = template*(trace)';
%     
% end

% convolution(convolution < 0) = 0;
% convolution = convolution/max(convolution);
convolution = smooth(convolution,floor(tau/dt))';
% convolution = [zeros(1,prev_samples) convolution(1:end-prev_samples)];

% colored_noise(colored_noise < 0) = 0;
% colored_noise = colored_noise/max(colored_noise);

% mean(convolution)
% 
% mean(colored_noise)
% 
% std(colored_noise)


conv_thresh = conv_thresh*std(colored_noise);

event_onsets = find(convolution(2:end-1) > convolution(1:end-2) &...
                    convolution(2:end-1) > convolution(3:end)   &...
                    thresh_vec(2:end-1));
                
event_onsets = event_onsets + 1;

% bad_event_idx = [];
% 
% for i = 1:length(event_onsets)
%     
%     diff_idx = min(event_onsets(i)+floor(.5*tau/dt),length(trace));
%     
%     slope = (trace(diff_idx) - trace(event_onsets(i)))/((diff_idx-event_onsets(i))*dt);
%     
%     if slope > -2000
%         bad_event_idx = [bad_event_idx i];
%     end
%     
% end

% event_onsets(bad_event_idx) = [];
        

% if plot_figs
%     norm_conv = convolution/max(convolution)*max(trace);
%     norm_conv_thresh = conv_thresh/max(convolution)*max(trace);
%     figure;
%     plot(t_vector,trace,'b',t_vector,norm_conv,'g');
%     hold on
%     plot(t_vector,norm_conv_thresh*ones(size(t_vector)),'r','LineWidth',3);
%     hold on
%     plot(t_vector,thresh_vec*max(trace),'m');
%     hold on
%     scatter(t_vector(event_onsets),max(trace)*ones(length(event_onsets),1),'r*')
% end