s_choice = [1, (crit.n_s + 1) / 2, crit.n_s];

for i = 1:3
    figure(i);
    subplot(2, 1, 1);
    plot(kgrid_dense, V(s_choice(i), :));
    title('value function');
    subplot(2, 1, 2);
    bar(para.moment_kgrid, Dist((s_choice(i) - 1)*crit.m_g(2)+1:s_choice(i)*crit.m_g(2)));
    title('distribution');
end

c = [0.25, 0.5, 1];
figure(4);
for i = 1:3
    subplot(3, 1, i);
    ss = (crit.n_s + 1) / 2;
    tmp = (c0_hat(ss, :) < c(i));
    kp = ka(ss, :) .* (1-tmp) + kn(ss, :) .* tmp;
    kp = kp - para.moment_kgrid';
    plot(para.moment_kgrid, kp);
    title(sprintf('Policy for realized shock c%d', i));
end