% demo code.

%% create the dataset
linear_combination

%% add the noise
noise = randn(size(X0));
noise = noise / norm(noise(:,1));
X = X0 + 0.01*noise;


%% PCT of modes to keep in SVD step (precise number not important)
pct = 0.45;
svd_rank_truncation = floor( rank(X)*pct/2 )*2;

%% select a random subset of steps to fit (so we don't have to use all)
fitsteps = unique([1, randi([1,size(X,2)], 1,10), size(X,2)] );

%% compute DMD without and with removal of frequencies
out = dmd(X, dt, svd_rank_truncation, 'rom_type', 'tlsq', 'step', fitsteps);


figure(1);

tiledlayout('flow');

% PLOT DISCRETE-TIME EIGENVALUES
nexttile 
plot(out.lambda,'ro','MarkerSize',9,'DisplayName','Pure DMD'); hold on; 
plot(exp(dt*frequencies),'+', 'Color',repmat(0.35,[1,3]),'DisplayName','True freq');

axis([-1,1,-1,1]);
legend('location','best');
title("Discrete-time (DMD matrix) eigenvalues")


nexttile;
plot(out.omega,'ro','MarkerSize',9,'DisplayName','Pure DMD'); hold on; 
plot(frequencies,'+','Color',repmat(0.35,[1,3]),'DisplayName','True freq');

% plot mode index
for k = 1:min(20, svd_rank_truncation)
    text(real(out.omega(k))+0.01,imag(out.omega(k))+0.1, " " + k, 'Color','r','clipping',true );    
end

xlim([-1,0.5])
ylim([-5,5]);
legend('location','best');
title("Continuous-time (flow) eigenvalues")
grid minor;


figure(2);
tiledlayout('flow');

for k = 1:min(20, svd_rank_truncation)
    
    nexttile
    plot(x, abs( out.Phi(:,k) )); hold all;    
    title({"Amplitude of DMD mode " + k; "Omega = " + out.omega(k)})
    nexttile
    plot(x, unwrap(angle( out.Phi(:,k) ))); hold all;    
    title({"Phase of DMD mode " + k; "Omega = " + out.omega(k)})
    
end

%% SHOW ROM MADE WITH PURE DMD
colorrange = [-1,1]*max(abs(X),[],'all');

%% reconstruct using all DMD modes (to convince ourselves that full DMD is roughly correct)
ROM_N = 10;

DataMatrixROM = reduce_order( out.Phi, out.omega, out.b, t, 1:ROM_N);

if norm(imag(DataMatrixROM)) > 0.1
    disp('Make sure to use both complex modes')
    ROM_N = ROM_N-1;
    DataMatrixROM = reduce_order( out.Phi, out.omega, out.b, t, 1:ROM_N);
end

figure(3); 
tiledlayout('flow'); 

nexttile; pcolor(X); colorbar; shading interp; caxis(colorrange);
title("Input data")

nexttile; pcolor(X0); colorbar; shading interp; caxis(colorrange);
title("Input data w/o noise")


nexttile; pcolor(real( DataMatrixROM ) ); colorbar; shading interp; caxis(colorrange);
title("ROM");

nexttile; pcolor( ( real( X0 ) - real( DataMatrixROM ) ) ); shading interp; colorbar;caxis(colorrange);
title("Difference between data w/o noise and ROM of order " + ROM_N );
