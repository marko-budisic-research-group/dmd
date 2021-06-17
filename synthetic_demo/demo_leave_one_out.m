%% CREATE SYNTHETIC DATASET
addpath('..');
% create the dataset
linear_combination

% add the noise
noise = randn(size(X0));
noise = norm(X0(:,1))* noise / norm(noise(:,1));
X = X0 + 0.01*noise;


%% COMPUTE DMD

% percentage of modes to keep in SVD step (precise number not important)

% select a random subset of steps to fit (so we don't have to use all)
fitsteps = unique([1, randi([1,size(X,2)], 1,10), size(X,2)] );

%% LEAVE ONE OUT VALIDATION
figure(1); tiledlayout('flow');

for pct = [0.1:0.2:0.9, 0.95]

svd_rank_truncation = floor( rank(X)*pct )

Ntests = 10;
DMD_eigenvalues = nan( svd_rank_truncation, Ntests);

for k = 1:Ntests
    
    % select a random snapshot pair to drop
    sel = randi(size(X,2)-1);
    
    
    % compute DMD without and with removal of frequencies
    out = dmd(X, dt, ...
        svd_rank_truncation, ... % dimension of reduced DMD matrix
        'rom_type', 'tlsq', ...  % standard or total DMD
        'dmd_type','exact',...
        'removefrequencies',0,...
        'sortby','l2',... % sort modes by |b| or mean L2
        'step', fitsteps, ... % snapshots over which to optimize value of b
        'excludepairs', sel...
        );
    
    
    DMD_eigenvalues(:,k) = out.lambda(:);
    
end

nexttile;

% visualize eigenvalues as 'heat plot'
[Xg,Yg] = meshgrid(-1:0.005:1,[-1:0.005:1]);
V = ksdensity( [real(DMD_eigenvalues(:)), imag(DMD_eigenvalues(:))],...
    [Xg(:),Yg(:)], ...
    'Bandwidth',0.025 );
pcolor(Xg,Yg,reshape(V,size(Xg))); shading flat
colormap(flipud(hot));
% unit circle
hold on; 
rectangle("Position",[-1,-1,2,2],"Curvature",1, "EdgeColor",'k'); 

% scattered eigenvalue data
hold on; h = scatter( real(DMD_eigenvalues(:)), imag(DMD_eigenvalues(:)),'filled');
alpha(h,0.1)
h.CData = [0,0.9,0];
h.SizeData = 16;
hold off;
axis equal square
title({"Leave-one-out test of robustness of eigenvalues; SVD dim=" +svd_rank_truncation})

end
