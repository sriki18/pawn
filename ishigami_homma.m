%{
Test PAWN algorithm by making call to runPAWN.m using the analytical
Ishingami-Homma function discussed in [1]. 

## Refernces
[1]: Pianosi, F., Wagener, T., 2015. A simple and efficient method for 
global sensitivity analysis based on cumulative distribution functions. 
Environ. Model. Softw. 67, 1–11. doi:10.1016/j.envsoft.2015.01.004
%}

clear;

% Ishigami-Homma function
f = @(x,p) sin(x(1)) + p.a * sin(x(2)^2) + p.b * x(3)^4 * sin(x(1));
p.a = 2; p.b = 1;
ih = @(x) f(x, p);

% Bounds
lb = ones(1,3) * -pi; ub = ones(1,3) * pi;

% Parameters from Figure 4's caption
n=15;
Nu = 100;
Nc = 100; % Is actually 50 in the paper

% Other parameters
npts = 100; seed = 4;

[KS,xvals,y_u, y_c, par_u, par_c, ft] = PAWN(f, p, lb, ub, Nu, n, Nc, npts, seed);

m1 = min([y_c, y_u']);
m2 = max([y_c, y_u']);
[f,ci] = ksdensity(y_u, linspace(m1,m2,npts), 'Function', 'cdf');

% Begin Plotting
subplot(331); ylabel('y'); hold on;box on;
subplot(334); ylabel('CDF'); hold on;  box on;
subplot(337); ylabel('KS'); hold on; xlabel('x_1'); box on;
subplot(338); xlabel('x_2'); hold on;box on;
subplot(339); xlabel('x_3'); hold on;box on;
for ind=1:length(lb)
    subplot(330+ind)
    plot(par_c(:,ind),y_c, 'xk', 'markersize', 2);
end

crit_c = [1.22,1.36,1.48,1.63,1.73,1.95]; % 0.1, 0.05, 0.025, 0.010, 0.005, 0.001
critval = crit_c(2) * sqrt((Nu+Nc)/(Nu*Nc));

for ind=1:length(lb)
    subplot(333+ind)
    plot([-100, 100], [1, 1], 'k--'); hold on
    plot(ci,ft((ind-1)*n+1:ind*n,:), 'color', [0.8,0.8,0.8]); ylim([0,1]);
    plot(ci,f, 'r:','linewidth',4); ylim([0,1.2])
    hold off
end

colData = linspace(0,n,n)'/n; colData = [colData colData colData];

for ind=1:length(lb)
    subplot(336+ind)
    [xtoplot, indices] = sort(xvals(ind,:));
    ytoplot = KS(ind, :); ytoplot = ytoplot(indices);
    plot([min(xvals(:)), max(xvals(:))], [critval,critval], 'k--'); hold on;
    plot(xtoplot, ytoplot, 'k');
    scatter(xtoplot, ytoplot, [],  colData, 'filled', 'markeredgecolor', 'k'); 
    ylim([0,1])
    hold off
end

disp(strcat('Median of KS : ', num2str(median(KS,2)')))
disp(strcat('Max of KS : ', num2str(max(KS'))))

% saveas(gcf, 'figures/fig4.png')
print('figures/fig4.png', '-dpng', '-r300')