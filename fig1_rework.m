%% USER PARAMETERS
clear variables

SAVETOFILE = 0; % CARE: Will overwrite files with same name!!!

res = 500; %resolution, #datapoints in each x,y-direction
%betakappa = 0.05;
k = 2;
x0 = 1;
theta = pi/6;

%c_0 = 3e8;
eps_R = 2;
mu_R = 2;
mu_I1 = 0.01;
mu_I2 = 0.04;
eps_I1 = 0.01;
eps_I2 = 0.04;

kappa = 0.5*k*sqrt(eps_R/mu_R)*(2*mu_R - 1i*(mu_I1-mu_I2));
%kappa = k + 1i*0.5*k*sqrt(eps_R/mu_R)*(mu_I2-mu_I1);
beta_mu = (mu_I1+mu_I2)/(2*mu_R - 1i*(mu_I1-mu_I2));
beta_eps = (eps_I1+eps_I2)/(2*eps_R - 1i*(eps_I1-eps_I2));
if beta_mu ~= beta_eps
    error("Impedance matching criteria not fulfilled!")
end
kappabeta_mu = kappa*beta_mu;
kappabeta_eps = kappa*beta_eps;

%% CALCULATION
x = linspace(-10,10,res);
y = transpose(linspace(-10,10,res));
E = exp(-kappabeta_mu*(x.*cos(theta) + y.*sin(theta))).*(cosh(x./x0)).^(1i*kappa*x0*cos(theta)).*exp(-1i*kappa*sin(theta).*y);

reE = real(E);
max_magn_E = max( abs(max(max(reE))), abs(min(min(reE))));

fprintf("max(|E|) = %f\n", max_magn_E)

% Convert to grid data
[X,Y] = meshgrid(x,y);
xi = min(x):(1/res):max(x);
yi = min(y):(1/res):max(y);
[X1, Y1] = meshgrid(xi,yi);
Z = griddata(X, Y, reE, X, Y);

%% FIGURE
figure()
imagesc(xi, yi, Z);
xlabel('$x$ (m)', 'Interpreter', 'latex')
ylabel('$y$ (m)', 'Interpreter', 'latex')
ax = gca;
set(ax, 'YDir', 'normal');
c = colorbar;%('XTick',;
caxis([-max_magn_E, max_magn_E])
set(get(c,'title'),'string','$E(x,y)$', 'Interpreter', 'latex');
ax.XAxis.FontSize = 14; % for presentations: 18
ax.YAxis.FontSize = 14;
c.FontSize = 14;


% Draw arrows
arrowlength = sqrt((.3*max(x))^2 + (.3*max(y))^2);
drawArrow([-arrowlength*2, -arrowlength], [-arrowlength*2 + arrowlength*cos(theta), -arrowlength + arrowlength*sin(theta)],'r'); %RHS
drawArrow([arrowlength*2, -arrowlength], [arrowlength*2 - arrowlength*cos(theta), -arrowlength + arrowlength*sin(theta)], 'b'); %LHS
xlim([min(x) max(x)])
ylim([min(y) max(y)])
xline(0)
yline(0)

colormap(brewermap([],'*RdBu')); % *-prefix reverses the color scheme

% SAVE FIGURE TO FILE:
if SAVETOFILE
    exportgraphics(ax,"analytic_grid_theta="+theta*180/pi+"_x0="+x0+"_beta="+betakappa+".png", 'Resolution', 600);
end

%% MATERIAL FUNCTION PLOT
x0 = [0.1; 1];
%eps = -eps_R*(tanh(x./x0) + 1i*beta_eps);
eps = -(eps_I1+eps_I2)/(2*beta_eps)*(tanh(x./x0) + 1i*beta_eps);

lwidth = 2.3;
legendfsize = 12;
labelXsize = 16;
labelYsize = 18;
axesfsize = 12;
titlefsize = 18;

figure(); 
sgtitle("$$\varepsilon_R = $$"+eps_R+"; $$\kappa\beta = $$"+kappabeta_eps, 'FontSize',titlefsize, 'Interpreter', 'latex');
grid on; hold on;
ax = gca;
ax.FontSize = axesfsize;
h = zeros(1,6);
h(5) = xline(0); h(6) = yline(0);
h(1:2) = plot(x, real(eps), 'LineWidth', lwidth);
h(3:4) = plot(x, imag(eps), 'k:', 'LineWidth', lwidth);

xlabel("$$x$$ (m)",'FontSize',labelXsize,'Interpreter','latex')
ylabel("$$\varepsilon/\varepsilon_0$$",'FontSize',labelYsize,'Interpreter','latex')
legend(h(1:2),"$$\varepsilon_R, x_0 = $$"+x0(1),"$$\varepsilon_R, x_0 = $$"+x0(2),'Interpreter','latex', 'FontSize', legendfsize)
yticks([-2 -1 0 1 2])
hold off
newcolors = {'#0072BD','#D95319'};%,'#77AC30'};
colororder(newcolors)
