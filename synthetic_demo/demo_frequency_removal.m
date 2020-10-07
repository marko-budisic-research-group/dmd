%% DEMO OF FREQUENCY REMOVAL BY HARMONIC AVERAGING
% This is a generalization of mean removal. See paper by Hirsh (2019)

%% create the dataset
linear_combination

%% add the noise
noise = randn(size(X0));
noise = noise / norm(noise(:,1));
X0 = X0 + 10;
X = X0 + 0.01*noise;

%% choose which frequencies to remove by harmonic averaging
toremove = [[1j,-1j]*2,0];


%% PCT of modes to keep in SVD step (precise number not important)
pct = 0.45;
svd_rank_truncation = floor( rank(X)*pct/2 )*2;

%% select a random subset of steps to fit (so we don't have to use all)
fitsteps = unique([1, randi([1,size(X,2)], 1,10), size(X,2)] );

%% compute DMD without and with removal of frequencies
out1 = dmd(X, dt, svd_rank_truncation, 'rom_type', 'tlsq','sortbyb',false, 'step', fitsteps);
out2 = dmd(X, dt, svd_rank_truncation, 'rom_type', 'tlsq','sortbyb',false, 'step', fitsteps, ...
    'removefrequencies', toremove );


%% COMPARE FREQUENCY SPECTRA
figure(1); clf;

% match DMD eigenvalues between two calculations
OmegaDistance = abs( out1.lambda - out2.lambda.' );
[M,uOm1, uOm2] = matchpairs( OmegaDistance, 0.05 );

tiledlayout('flow');

% PLOT DISCRETE-TIME EIGENVALUES
nexttile 
plot(out1.lambda,'ro','MarkerSize',9,'DisplayName','Pure DMD'); hold on; 
plot(out2.lambda,'bx','MarkerSize',9,'DisplayName','DMD w/ mode removal');
plot(exp(dt*frequencies),'*', 'Color',repmat(0.35,[1,3]),'DisplayName','True freq');

% plot connection lines
for k = 1:size(M,1)    
    line( real([out1.lambda(M(k,1)); out2.lambda(M(k,2))]),...
        imag([out1.lambda(M(k,1)); out2.lambda(M(k,2))]),'HandleVisibility','off' );    
end
axis([-1,1,-1,1]);
legend('location','best');
title("Discrete-time (DMD matrix) eigenvalues")


nexttile;
plot(out1.omega,'ro','MarkerSize',9,'DisplayName','Pure DMD'); hold on; 
plot(out2.omega,'bx','MarkerSize',9,'DisplayName','DMD w/ mode removal');
plot(frequencies,'*','Color',repmat(0.35,[1,3]),'DisplayName','True freq');

% plot mode index
for k = 1:20
    text(real(out1.omega(k))+0.01,imag(out1.omega(k))+0.1, " " + k, 'Color','r','clipping',true );
    text(real(out2.omega(k))-0.01,imag(out2.omega(k))-0.1, " " + k, 'Color','b','clipping',true );
    
end
% plot connection lines
for k = 1:size(M,1)    
    line( real([out1.omega(M(k,1)); out2.omega(M(k,2))]),...
        imag([out1.omega(M(k,1)); out2.omega(M(k,2))]),'HandleVisibility','off' );    
end

xlim([-1,0.5])
ylim([-5,5]);
legend('location','best');
title("Continuous-time (flow) eigenvalues")
grid minor;

%% SHOW ROM MADE WITH PURE DMD
colorrange = [-1,1]*max(abs(X),[],'all');

%% reconstruct using all DMD modes (to convince ourselves that full DMD is roughly correct)
DataMatrixROM1 = reduce_order( out1.Phi, out1.omega, out1.b, t);

figure(3); 
tiledlayout('flow'); 

nexttile; pcolor(X); colorbar; shading interp; caxis(colorrange);
title("Input data")

nexttile; pcolor(X0); colorbar; shading interp; caxis(colorrange);
title("Input data w/o noise")


nexttile; pcolor(real( DataMatrixROM1 ) ); colorbar; shading interp; caxis(colorrange);
title("ROM");

nexttile; pcolor( ( real( X0 ) - real( DataMatrixROM1 ) ) ); shading interp; colorbar;caxis(colorrange);
title("Difference between ROM and data w/o noise");


%% SHOW ROM MADE USING DMD with MODE REMOVAL

%% reconstruct using all DMD modes (to convince ourselves that full DMD is roughly correct)
%% notice that we manually add modes that were computed by harmonic averaging
DataMatrixROM2 = reduce_order( [out2.Phi, out2.AvgPhi], [out2.omega;toremove(:)], [out2.b;out2.AvgB(:)], t);

figure(4); 

tiledlayout('flow'); 

nexttile; pcolor(X); colorbar; shading interp; caxis(colorrange);
title("Input data")

nexttile; pcolor(X0); colorbar; shading interp; caxis(colorrange);
title("Input data w/o noise")


nexttile; pcolor(real( DataMatrixROM2 ) ); colorbar; shading interp; caxis(colorrange);
title("ROM");

nexttile; pcolor( ( real( X0 ) - real( DataMatrixROM2 ) ) ); shading interp; colorbar;caxis(colorrange);
title("Difference between ROM and data w/o noise");

%% MODE PROFILE

figure(5);
tiledlayout('flow');

for k = 1:numel(toremove)
    
    [~,idx1] = min( abs(out1.omega - toremove(k)))
    [~,idx2] = min( abs(out2.AvgOmega - toremove(k)))
    
    nexttile
    plot(x, abs( out1.Phi(:,idx1) ), 'r--', 'DisplayName','DMD-computed mode'); hold all;
    plot(x, abs( out2.AvgPhi(:,idx2) ), 'b:', 'DisplayName','Harmonic averaged mode');
    legend('location','best')
    
    title("Amplitude of removed mode " + k)
    nexttile
    plot(x, unwrap(angle( out1.Phi(:,idx1) )), 'r--','DisplayName','DMD-computed mode'); hold all;
    plot(x, unwrap(angle( out2.AvgPhi(:,idx2))), 'b:', 'DisplayName','Harmonic averaged mode');
    legend('location','best')
    title("Phase of removed mode " + k)
    
end
