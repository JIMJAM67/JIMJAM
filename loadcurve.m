data=[0 6 20;
      6 10 25;
      10 12 30;
      12 16 25;
      16 20 35;
      20 24 20];

power = data(:,3);
Dt = data(:,2) - data(:,1);

Total_power = sum(power .* Dt);

Average_load = Total_power / sum(Dt);
peak_load = max(power);
Daily_LF = Average_load/peak_load*100;

% Display only required results
fprintf('Average Load = %.2f MW\n', Average_load);
fprintf('Peak Load = %.2f MW\n', peak_load);
fprintf('Daily Load Factor = %.2f %%\n', Daily_LF);

% Plot load curve
L = length(power);
timeinterval = data(:,1:2);
t = sort(reshape(timeinterval,1,2*L));

P = zeros(1,2*L);
for n = 1:L
    P(2*n-1) = power(n);
    P(2*n) = power(n);
end

figure;
plot(t,P,'LineWidth',1.5);
xlabel('Time, Hr');
ylabel('Power, MW');
title('Generation Power, MW versus Time, Hour');
grid on;