%% Using DMD on Lorenz'96

%% Running the model
% This section runs the Lorenz'96 and stores the result in a data matrix
% suitable for input into DMD.

% Time grid
dt = 0.05;
T=0:dt:120;

% state vector dimension
N=100;

% initial condition
IC = 1./(1:N).';

% run the simulation
[t,y] = ode45( @(t,y)l96(t,y,3), T, IC );
y = y'; % transpose it - columns are state vectors

% % remove the transient
% y = y(:, T>60);
% t = t(T>60);


% spatiotemporal matrix
figure(1)
tiledlayout('flow');
nexttile;
pcolor(t, 1:N, y); shading flat
colormap(gca, 'parula')
colorbar;
caxis([-1.1,1.1]*max(abs(caxis)));
colorrange=caxis;
title("Input data");

%% DMD analysis
% We run the DMD here to compute DMD modes. Frequencies are converted into
% decay times and periods (this is completely analogous to how you'd
% interpret eigenvalues of x' = Ax linear system).

addpath('./dmd') % this is my code folder

% compute DMD without and with removal of frequencies
desired_rank = 50;
fitsteps = ceil( linspace(1, size(y,2), 5) ); % five equispaced steps to fit over

out = dmd(y, dt, ...
    desired_rank, ... % dimension of reduced DMD matrix
    'rom_type', 'tlsq', ...  % standard or total DMD
    'sortbyb',false,... % sort modes by |b| or mean L2 (default:false)
    'step', fitsteps ... % snapshots over which to optimize value of b
    );

% out.meanL2norm = [out.AvgMeanL2norm; out.meanL2norm];
% out.b = [out.AvgB; out.b];
% out.omega = [out.AvgOmega; out.omega];
% out.lambda = [out.AvgLambda; out.lambda];
% out.Phi = [out.AvgPhi, out.Phi];


out.periodDecay = log(2)./real( out.omega);
out.periodOsc   = 2*pi./imag( out.omega );

%%
% Label modes by their decay/doubling times and oscillation periods
modelabels = {};
for k = 1:numel(out.omega)
modelabels{k} = sprintf(...
    "TC %.2f -- P %.2f",  ...
    out.periodDecay(k), ... % computing doubling time from real part
    out.periodOsc(k) );    % computing period from imaginary part
end


%% Reconstruction
% Here we use DMD modes to produce reconstruction of the input data.
romrank = 12;

% In this case, we're using only top modes (most energetic), but any subset
% of indices can be used.
mode_selection = 1:romrank;

try
ROM = reduce_order( out.Phi, out.omega, out.b, t, mode_selection );
assert ( max(imag(ROM),[],'all') < 1e-8, ...
    'UNPAIRED',...
    "Bad reconstruction - did you unpair conjugate modes by accident?")
catch e
    if strcmpi( e.identifier, 'UNPAIRED' )
        disp("Reducing rank by one to make sure not to unpair conjugate modes");
        romrank = romrank-1;
        mode_selection = 1:romrank;
        ROM = reduce_order( out.Phi, out.omega, out.b, t, mode_selection);
        assert ( max(imag(ROM),[],'all') < 1e-8, "Bad reconstruction - did you unpair conjugate modes by accident?")
    end
end

ROM = real(ROM);

%% Visualizations
% First, let's compare spatiotemporal fields.
figure(1)

nexttile;
pcolor(t, 1:N, ROM ); shading flat
colorbar;
caxis(colorrange);
title("ROM rank = " + romrank );

ax2=nexttile;
pcolor(t, 1:N, abs(y-ROM)); shading flat
colorbar;
caxis([0,3]);
colormap(ax2,parula(3));
title("Error of ROM rank = " + romrank );



%% Then, let's look at eigenvalue plots.
figure(2);

tiledlayout('flow');
nexttile;
scatter( real( out.omega), imag( out.omega ), 2*(1 + out.meanL2norm).^2, sign( real(out.omega) ),'filled');
colormap( gca, flipud( eye(3) ) );
title( colorbar, 'Sign of re(evalue)')
title('Continuous-time eigenvalues');
ylabel('Im(omega)  - ang. freq')
xlabel('Re(omega)  - decay coeff');

nexttile;
scatter( out.periodDecay, out.periodOsc, 2*(1 + out.meanL2norm).^2, sign( real(out.omega) ),'filled');
colormap( flipud( eye(3) ) );
title( colorbar, 'Sign of re(evalue)')
title('Continuous-time eigenvalues');
ylabel('Oscillation period')
xlabel('Doubling time (+), Half-life time (-)');

nexttile;
stem(imag( out.omega ), out.meanL2norm );
title('Continuous-time eigenvalues');
ylabel('Im(omega)  - ang. freq')
ylabel('Mean L2 norm (in time)');

nexttile;
stem(real( out.omega ), out.meanL2norm );
title('Continuous-time eigenvalues');
xlabel('Re(omega)  - decay coeff');
ylabel('Mean L2 norm (in time)');

%% Finally, spatial profiles associated with modes and their individual time evolution.
linecolors = lines(5);

figure(3);
tiledlayout('flow');
nexttile;
ROM1 = reduce_order( out.Phi, out.omega, out.b, t, 1 );
pcolor(t, 1:N, real( ROM1 ) ); shading flat
caxis([0,2]);
colorbar;
title("Modes 1 " + modelabels{1});


nexttile;
plot(1:N, real( out.b(1)*out.Phi(:,1) ),'-', 'Color',linecolors(1,:), 'DisplayName',modelabels{1}, 'LineWidth',3 );
title('First mode (non-oscillatory, mean flow)');
xlabel('State number')
ylim([-1.1,1.1]*max(abs(ylim)));

legend;

figure(4);
tiledlayout('flow');
nexttile;
ROM23 = reduce_order( out.Phi, out.omega, out.b, t, 2:3 );
pcolor(t, 1:N, real( ROM23 ) ); shading flat
colorbar;
title("Modes 2,3 " + modelabels{2});

nexttile;
plot(1:N, real( out.b(2)*out.Phi(:,2) ),'-', 'Color', linecolors(2,:), 'DisplayName',"RE " + modelabels{2}, 'LineWidth',3 ); hold on;
plot(1:N, imag( out.b(2)*out.Phi(:,2) ),'-', 'Color',linecolors(4,:), 'DisplayName',"IM " + modelabels{2}, 'LineWidth',3 ); hold off;
title('Second mode pair (oscillatory)');
xlabel('State number')
ylim([-1.1,1.1]*max(abs(ylim)));

legend;
hold off;




%% AUXILIARY - L96 model

function dx = l96(t, x, F)
% Evaluates the right hand side of the Lorenz '96 system
% \frac{dx_i}{dt} = -x_{i-2}x_{i-1} + x_{i-1}x_{i+1} - x_i + F

DIM = size(x,1);

dx = zeros(DIM,1);

for j=3:DIM-1
    dx(j) = -x(j)+x(j-1)*(x(j+1)-x(j-2)) + F;
end
dx(1) = -x(1)+x(DIM)*(x(2)-x(DIM-1)) + F;
dx(2) = -x(2)+x(1)*(x(3)-x(DIM)) + F;
dx(DIM) = -x(DIM)+x(DIM-1)*(x(1)-x(DIM-2)) + F;

end
