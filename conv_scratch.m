t = 0:.01:100;
tau_decay = 10; decay = exp(-t/tau_decay);
tau_rise = 1; rise = -exp(-t/tau_rise);
figure; plot(tau_decay); hold on; plot(tau_rise)
figure; plot(decay); hold on; plot(rise)
hold on; plot(decay + rise)
title('difference of exponentials')
help conv
kernel = (decay + rise)/max(decay + rise);
hold on; plot(kernel)
input = [zeros(1,2000) ones(1,4000) zeros(1,4001)];
output = conv(input,kernel);
figure; plot(output)
input = [zeros(1,2000) ones(1,7000) zeros(1,1001)];
output = conv(input,kernel);
figure; plot(output)
kernel = (decay + rise)/sum(decay + rise);
output = conv(input,kernel);
figure; plot(output)