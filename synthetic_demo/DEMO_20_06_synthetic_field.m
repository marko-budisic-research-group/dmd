%% Synthetic field for DMD analysis
% Each of the constituents will be given its own dynamics, in order to present
% a clear target for DMD.
%
% Explanation screencast: https://youtu.be/TcKInGLe3M0
%
% For ease of visualization, the domain will be 1D, which shouldn't really affect
% the analysis.
%
% The list of constituents
%%
% # Oscillating component with weak growth/decay (2 states)
% # Slowly decaying monotonic (non-oscillating) component (1 state)
% # Oscillating component with linear (non-normal) growth (4 states)
% # Non-sinusoidal oscillating component (2 states)
% # Oscillating component undergoing a Hopf bifurcation (2 states)
% # Background mode undergoing a pitchfork bifurcation (1 state)
% # Additive (non-gaussian) noise
%%
% Total: 12 states+noise
%
% Dynamics are simulated from a 12-component initial condition which can then
% be used to observe the pattern in each of the components. Passing non-zero value
% to a component initial condition makes it possible to observe its pattern.


%% Parameters:
TimeScale = 10;
dt = 1e-1; % timestep
dx = 1e-3; % spatial step
TimeDomain = 0:dt:3*TimeScale;
SpaceDomain = -1:dx:1;

NoiseMagnitudeGain = .5; % 50% of noise


x0 = nan(12,1); % initial condition
component_names = strings(12,1); % names of components

% # Oscillating component with weak growth/decay (2 states)

p.Comp1.DoublingTime = -TimeScale*1.56;
p.Comp1.Period = TimeScale/5;
x0(1:2) = [1;1];
component_names(1:2) = "Weak growth/decay osc.";

% # Slowly decaying monotonic (non-oscillating) component (1 state)

p.Comp2.DoublingTime = -TimeScale;
x0(3) = 1;
component_names(3) = "Monotonic";

% # Oscillating component with linear (non-normal) growth (4 states)

p.Comp3.DoublingTime = inf;
p.Comp3.Period = TimeScale/7;
x0(4:5) = 0.1;
component_names(4:5) = "Resonating/Periodic";
x0(6:7) = 0.1;
component_names(6:7) = "Res. driver (periodic)";

% # Non-sinusoidal oscillating component (2 states)

p.Comp4.DoublingTime = inf;
p.Comp4.Period = TimeScale/1.8;
x0(8:9) = 0.5;
component_names(8:9) = "Non-sin. osc";

% # Oscillating component undergoing a Hopf bifurcation (2 states)

p.hopf = 0.5;
p.hopf_period = TimeScale/6.28;
x0(10:11) = 0.1;
component_names(10:11)= "Limit cycle (Hopf)";

% # Background mode undergoing a pitchfork bifurcation (1 state)

p.pitchfork = 0.25;
x0(12) = 1;
component_names(12) = "Bkg. pfork";




% NoiseMagnitude


%% Simulate coefficients
[t,x] = ode45( @(t,x)rhs_synthetic(t,x,p), TimeDomain, x0 );

% add harmonics to signals 8 and 9
clip = @(x,c) x - (x-sign(x).*c).*(abs(x) > c); % symmetric clip of the signal at value c
x(:, 8:9) = clip( x(:,8:9), max(abs(x0(8:9)))/2); % clip at half-value of initial condition


%% Spatial field functions
% All spatial fields will be evaluated on [-1,1] domain, but they are
% specified using formulas (to allow for finer/coarser sampling)
% We need 12 such fields (as we have 12 coefficients)

fieldfn = {};
fieldfn{1} = @(x)cos(2*pi*x-0.5);
fieldfn{2} = @(x)sin(2*pi*x-0.5);

fieldfn{3} = @(x)pdf('norm',x, 0,0.5);

fieldfn{4} = @(x)triangularPulse( -1,-0.8, -0.6,x )*(-1);
fieldfn{5} = @(x)triangularPulse( -0.8, -0.6, -0.4, x );
fieldfn{6} = @(x)triangularPulse( 0.4,0.6, 0.8,x )*(-1);
fieldfn{7} = @(x)triangularPulse( 0.6, 0.8, 1, x );

fieldfn{8} = @(x)pdf('norm',x, 0.4,0.1)*(-1);
fieldfn{9} = @(x)pdf('norm',x, 0.6,0.1);

fieldfn{10} = @(x)pdf('lognorm',x+1, 0,0.25);
fieldfn{11} = @(x)pdf('lognorm',-x+2, 0,0.25);

fieldfn{12} = @(x)pdf('norm',x,0,0.5).*cos(2*pi*x);

SpatialFields = nan( numel(SpaceDomain), numel(fieldfn) );
SpatiotemporalFields = nan( numel(SpaceDomain), numel(TimeDomain), numel(fieldfn) );

for k = 1:numel(fieldfn)
    % evaluate the field function on a grid
    SpatialFields(:,k) = fieldfn{k}(SpaceDomain);
    SpatiotemporalFields(:,:,k) = SpatialFields(:,k) * x(:,k)';
end


%% Assembling the data set
DataWoNoise = sum(SpatiotemporalFields,3);

%% Add noise
% Gaussian noise with variance scaled so that the L2 norm of the total
% noise integrated over the spatiotemporal domain is exactly ==
% NoiseMagnitudeGain * norm(data w/o noise)

NoiseComponent = NoiseMagnitudeGain*norm(DataWoNoise,'fro')...
    *randn( numel(SpaceDomain), numel(TimeDomain) )*sqrt(dt/range(TimeDomain))*sqrt(dx/range(SpaceDomain));


DataAssembled = DataWoNoise + NoiseComponent;

DEMO_20_06_synthetic_field_COMPLETE = true;

try
    assert(any(x0 ~= 0));
    IncludedComponents =  "Included: " + join( unique( component_names( x0 ~= 0 ) ), ', ' );
catch
    IncludedComponents = "Included: ";
end
try
    assert(any(x0 == 0));
    ExcludedComponents =  "Not included: " + join( unique( component_names( x0 == 0 ) ), ', ' );
catch
    ExcludedComponents = "Not included: ";
end



%% VISUALIZATION

set(0,'DefaultFigureWindowStyle','docked') % dock all figures

%% Visualize Spectrum
figure(97); clf;
[~,M] = rhs_synthetic(0,x0,p);

E = eig(M);
scatter( real(E), imag(E), 64, 'filled' );

xline( log(2)/p.Comp1.DoublingTime, 'b--', 'C1' );
xline( log(2)/p.Comp2.DoublingTime, 'b--', 'C2' );
xline( log(2)/p.Comp3.DoublingTime, 'b--', 'C3' );
xline( log(2)/p.Comp4.DoublingTime, 'b--', 'C4' );
xline( p.hopf, 'b--', 'Hopf' );
xline( p.pitchfork, 'b--', 'PF' );

yline(2*pi/p.Comp1.Period, 'r--', 'C1' );
yline(0, 'r--', 'C2' );
yline(2*pi/p.Comp3.Period, 'r--', 'C3' );
yline(2*pi/p.Comp4.Period, 'r--', 'C4' );
yline(2*pi/p.hopf_period, 'r--', 'Hopf' );

yline(-2*pi/p.Comp1.Period, 'r--', 'C1' );
yline(0, 'r--', 'PF' );
yline(-2*pi/p.Comp3.Period, 'r--', 'C3' );
yline(-2*pi/p.Comp4.Period, 'r--', 'C4' );
yline(-2*pi/p.hopf_period, 'r--', 'Hopf' );

title("Spectrum of the linearization of the coefficient system Jacobian")
xlabel("Re(\omega)")
ylabel("Im(\omega)")
%% Visualize coefficients

figure(98); clf;
h_coeff = tiledlayout('flow');
for k = 1:12
    nexttile;
    plot(t, x(:,k), 'DisplayName',"C "+k);
    title("C " + k);
end
title(h_coeff, "Time evolution of combination coefficients");

%% Visualize individual components
cmap = [winter(8);flipud(autumn(8)) ];

for k = 1:numel(fieldfn)
    figure(k); clf;
    subplot(2,2,1);
    plot(SpaceDomain, SpatialFields(:,k) );
    ylim([-1,1]*max(abs(SpatialFields(:,k))));
    view([90 -90]);
    xlabel('Space');
    title('Spatial profile');
    
    subplot(2,2,2);
    contourf(t,SpaceDomain, SpatiotemporalFields(:,:,k) );
    colormap(cmap);shading flat;
    hcb = colorbar('east');
    xlabel('Time'); ylabel('Space');
    title('Spatiotemporal field');
    caxis([-1,1]*max(abs( SpatiotemporalFields(:,:,k)),[],'all'));
    
    
    subplot(2,2,3);
    [power, phase, F] = spectral_fft( x(:,k), dt );
    stem( F,power );
    xlabel('Frequency [1/time]');
    ylabel('Power');
    set(gca,'yscale','log');
    title('Power Spectrum of the time coefficient');
    
    subplot(2,2,4);
    plot(t, x(:,k), 'DisplayName',"C "+k);
    xlabel('Time');
    title('Time coefficient');
    
    sgtitle(component_names{k} + "; x(0) = " + x0(k));
end

%% Visualize assembled field
clim = max(abs(bounds(DataWoNoise,'all')))*[-1,1];

figure(99);

subplot(2,2,1);
contourf(t,SpaceDomain, DataWoNoise );
colormap(cmap);shading flat;
hcb = colorbar('east'); set(hcb,'xcolor','r','ycolor','r');
xlabel('Time'); ylabel('Space');
title('No noise')
caxis(clim);

subplot(2,2,2);
pcolor(t,SpaceDomain, NoiseComponent );
colormap(cmap);shading flat;
hcb = colorbar('east'); set(hcb,'xcolor','r','ycolor','r');
xlabel('Time'); ylabel('Space');
title('Noise');
caxis(clim);

subplot(2,2,3);
DataMean = repmat( mean(DataAssembled,2), size(TimeDomain) );
pcolor(t,SpaceDomain, DataMean );
colormap(cmap);shading flat;
hcb = colorbar('east'); set(hcb,'xcolor','r','ycolor','r');
xlabel('Time'); ylabel('Space');
title('Time-mean');
caxis(clim);


subplot(2,2,4);
pcolor(t,SpaceDomain, DataAssembled );
colormap(cmap);shading flat;
hcb = colorbar('east'); set(hcb,'xcolor','r','ycolor','r');
xlabel('Time'); ylabel('Space');
title('Fully assembled field');
caxis(clim);

sgtitle({"Input data; TimeScale = " + TimeScale + " Duration = "+max(TimeDomain); ...
    IncludedComponents; ExcludedComponents;
    },'horizontalAlignment', 'left' );

