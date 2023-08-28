change phase res = 1000; %resolution
b = 0.05;
k = 2;
a = 1;
theta = pi/6;
x = linspace(-12,12,res);
y = pi/2;
E = exp(-b*k*(x.*cos(theta) + y.*cos(theta))).*(cosh(x./a)).^(1i*k*a*cos(theta)).*exp(-1i*k*sin(theta).*y);

figure()
plot(x, real(E), 'LineWidth', 2)
grid on
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 22)
ylabel('$E(x,y=\frac{\pi}{2})$', 'Interpreter', 'latex', 'FontSize', 22)
xlim([min(x) max(x)]);
xline(0)
yline(0)
yticks([-3,-2,-1,0,1,2,3]);
yticklabels({'', '-2','-1','0','1','2',''})
ylim([-2.5, 2.5])
%export_fig test.png
saveas(ax,'yfixed.png')