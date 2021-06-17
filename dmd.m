function out = dmd( DataMatrix, dt, rom_dim, varargin )
%% DMD Dynamic Mode Decomposition (exact version)
%
% out = dmd( DataMatrix, dt, rom_dim, ... )
%
% Inputs:
%     DataMatrix (:,:) double {mustBeNumeric, mustBeReal, mustBeFinite} -- INPUT SNAPSHOT MATRIX (columns are snapshots)
%     dt (1,1) double {mustBePositive, mustBeFinite} -- INPUT SNAPSHOT TIMESTEP (columns are snapshots)
%     rom_dim (1,1) double {mustBePositive,mustBeInteger,mustBeFinite} -- REQUESTED DIMENSION OF REDUCED-ORDER MODEL
%
% Additional optional name-value pairs, specified as
% out = dmd( DataMatrix, dt, rom_dim, NAME, VALUE, NAME, VALUE,... )
%
%     out = dmd( DataMatrix, dt, rom_dim, 'rom_type', VALUE)
%         VALUE = 'lsq' - standard least-squares truncation {default}
%         VALUE = 'tlsq' - total least squares truncation (Hemati 2017)
%
%     out = dmd( DataMatrix, dt, rom_dim, 'dmd_type', VALUE)
%         VALUE = 'exact' Exact DMD style eigenvalue calculation
%         VALUE = 'rrr' (Drmac et al. 2018) DDMD_RRR optimization {default}
%
%    out = dmd( DataMatrix, dt, rom_dim, 'sortby', VALUE)
%         VALUE = 'initcond' - sort DMD modes by |b|
%         VALUE = 'l2' - sort DMD modes by average L2 norm of time
%                        coefficient
%         VALUE = 'residual' - sort DMD modes by optimal residual {default}
%
%
%    out = dmd( DataMatrix, dt, rom_dim, 'svdcode', VALUE)
%         VALUE = 'QR' - use QR-SVD from Lapack for SVD algorithm {default}
%                      - requires (automatic) mex compilation on first run
%         VALUE = 'DD' - use DD-SVD from MATLAB for SVD algorithm
%
%    out = dmd( DataMatrix, dt, rom_dim, 'normalize', VALUE)
%         VALUE = true - normalize snapshot matrices by column L2 norms,
%                        to avoid scaling effects on POD
%
%    out = dmd( DataMatrix, dt, rom_dim, 'step', VALUE)
%         VALUE is a row-vector of snapshot indices at which optimization is
%         performed to determine b coefficients in the expansion.
%         If VALUE = -1, all snapshots are used.
%
%    out = dmd( DataMatrix, dt, rom_dim, 'numericalRankTolerance', VALUE)
%         VALUE is the threshold for ratios sigma_n/sigma_1 of singular values
%               used to determine the numerical rank/"significance" of the POD subspace
%
%    out = dmd( DataMatrix, dt, rom_dim, 'excludepairs', VALUE)
%         VALUE is either a logical index vector or the list of numerical
%               index of snapshot pairs that will be excluded in DMD
%               calculation.  This option allows including (where VALUE=false)
%               or dropping (where VALUE=true, or an integer) pairs of snapshots in X1 and
%               X2 matrices, without the user needing to modify the input
%               matrix, which can be useful, e.g., for leave-one-out testing of variance.
%
%    out = dmd( DataMatrix, dt, rom_dim, 'ritzMaxIteration', IT, ...
%                 'ritzATOL', ATOL, 'ritzRTOL', RTOL )
%
%                 These parameters
%                 control the stopping criterion for the residual optimization
%                 loop. The loop will perform at most IT iterates, it stops
%                 sooner if subsequent Ritz values are within RTOL/ATOL range.
%
%    out = dmd( DataMatrix, dt, rom_dim, 'step', VALUE)
%         VALUE is a row-vector of snapshot indices at which optimization is
%         performed to determine b coefficients in the expansion.
%         If VALUE = -1, all snapshots are used.
%
%
%    out = dmd( DataMatrix, dt, rom_dim, 'removefrequencies', VALUE)
%         VALUE = vector of continuous-time complex-valued frequencies to
%                 remove from data by harmonic averaging before DMD is
%                 applied.
%
%                 If the passed vector is non-empty, output structure will
%                 contain vectors out.AvgOmega, out.AvgModes, etc. that
%                 correspond to modes that were computed by harmonic
%                 averaging. They are normalized and scaled in the same way
%                 as DMD modes, so they can be used in the same manner.
%
%
% Outputs:
% out.meanL2norm - time-average of L2 norm of mode
% out.b - contribution of the mode to 0th timestep
% out.Phi - DMD mode vector
% out.omega - continuous time DMD eigenvalue omega = log( lambda ) / dt
% out.lambda - discrete time DMD eigenvalue lambda = exp( omega * dt )
% out.model_rank - rank of the model (= input r parameter)
% out.optimalResiduals - residual after adjustment (if DDMD-RRR was NOT used, it will be populated by NaN)
%

p = inputParser;

p.addRequired('DataMatrix',@(M)validateattributes(M, {'numeric'},{'2d','real','finite'})  );
p.addRequired('dt',@(s)validateattributes(s, {'numeric'},{'scalar','positive','finite'}) );
p.addRequired('rom_dim', @(s)validateattributes(s, {'numeric'},{'scalar','positive','finite','integer'}));

p.addParameter('rom_type', 'lsq', @(x)assert(ismember(x,{'lsq','tlsq'}),'Input not recognized: %s',x) );
p.addParameter('dmd_type', 'exact', @(x)assert(ismember(x,{'exact','rrr'}),'Input not recognized: %s',x) );

p.addParameter('removefrequencies', [], @(s)validateattributes(s,  {'numeric'},{'row'}));
p.addParameter('sortby', 'residual', @(x)assert(ismember(x,{'initcond','l2','residual'}),'Input not recognized: %s',x) );

p.addParameter('step',1,@(s)validateattributes(s, {'numeric'}, {'row','integer','finite'}) );
p.addParameter('normalize',true,@(s)validateattributes(s, {'logical'}, {'scalar'}))

p.addParameter('excludepairs',[],@(s)validateattributes(s, {'logical','numeric'}, {'row', 'finite'}))

p.addParameter('ritzMaxIteration',1,@(s)validateattributes(s, {'numeric'}, {'scalar','integer','nonnegative'}) );

p.addParameter('ritzATOL',1e-4,@(s)validateattributes(s, {'numeric'}, {'scalar','nonnegative'}) );
p.addParameter('ritzRTOL',1e-2,@(s)validateattributes(s, {'numeric'}, {'scalar','nonnegative'}) );

p.addParameter('numericalRankTolerance',0.0,@(s)validateattributes(s, {'numeric'}, {'scalar','nonnegative','finite'}) );

p.addParameter('svdcode', 'QR', @(x)assert(ismember(x,{'QR','DD'}),'Input not recognized: %s',x) );

parse(p,DataMatrix, dt, rom_dim, varargin{:});
options = p.Results;


% use QR-SVD if available (see Drmac et al 2018)
if exist('svd_lapack') && strcmpi( options.svdcode,'QR')
    disp('Using QR SVD (Lapack GESVD)')
    mysvd = @(x)svd_lapack(x, 0,'gesvd');
else
    disp('Using DD SVD (built-in SVD)')
    mysvd = @(x)svd(x, 0);
end

if options.step == -1
    options.step = 1:size(DataMatrix,2);
end

%% Preprocessing
Nsnapshots = size(DataMatrix,2);
T = (Nsnapshots-1)*dt;

X1=DataMatrix(:,1:end-1); % likely a very slow part of the code
X2=DataMatrix(:,2:end);


%% Mean/frequency removal
if ~isempty(options.removefrequencies)
    disp('Removing predetermined frequencies')

    % row vector of continuous frequencies
    OmegaR = reshape( options.removefrequencies, 1, [] );

    % discrete frequencies row vector
    LambdaR = exp( OmegaR * dt );

    % Vandermonde
    Vandermonde = LambdaR .^ transpose( 0:( Nsnapshots-1 ) );

    % truncated vandermonde (for dealing with X1 and X2)
    VR_o = Vandermonde(2:end,:);

    % projectors
    PI = pinv(transpose(Vandermonde))*transpose(Vandermonde);
    PI_o = pinv(transpose(VR_o))*transpose(VR_o);

    X1 = X1 - X1*PI_o;
    X2 = X2 - X2*PI_o;


    % the following two calculations are in-principle equivalent, but not
    % exactly when number of time steps is finite
    HarmonicAverage = DataMatrix*pinv(transpose(Vandermonde));
    %HarmonicAverage = DataMatrix*conj(Vandermonde)/Nsnapshots;


    [out.AvgPhi, out.AvgB] = normalizeModes(HarmonicAverage);

    out.AvgOmega = options.removefrequencies(:);
    out.AvgLambda = LambdaR(:);
    out.AvgB = out.AvgB(:);

    meanL2 = abs(out.AvgB) .* sqrt( (exp(2*real(out.AvgOmega)*T)-1)./(2*real(out.AvgOmega)*T) );
    meanL2(isnan(meanL2)) = abs(out.AvgB(isnan(meanL2)));

    out.AvgMeanL2norm = meanL2;

    clear PI_o;
    clear PI;

end

%% Subselect pairs
if ~isempty(options.excludepairs)
    X1(:,options.excludepairs) = [];
    X2(:,options.excludepairs) = [];
end

% pre-normalization
if options.normalize
    X1 = normalize(X1, 'norm',2);
    X2 = normalize(X2, 'norm',2);
end

% total-DMD Hemati debias
switch(options.rom_type)
    case 'tlsq'
        disp("Total least-squares truncation to order " + rom_dim)
        Z = [X1;X2];
        [~,~,Q] = mysvd(Z);
        Qr = Q(:,1:rom_dim);
        X1 = X1*Qr;
        X2 = X2*Qr;
        clear Z Qr;
    case 'lsq'
        disp("Standard least-squares truncation to order " + rom_dim)
end

%% SVD subspace selection
[U,Sigma,V] = mysvd( X1 );
sigma = diag(Sigma);
numericalRank = find( sigma >= sigma(1)*options.numericalRankTolerance, 1, 'last' );
subspaceSize = min(rom_dim, numericalRank);
if subspaceSize < rom_dim
    warning("Numerical rank = " + numericalRank + " is smaller than requested ROM dimension = " + rom_dim);
end
Ur = U(:, 1:subspaceSize);
Vr = V(:,1:subspaceSize);
Sigmar = Sigma(1:subspaceSize, 1:subspaceSize);

clear U V Sigma;
clear X1;

%%

switch(options.dmd_type)
    case 'exact'
        disp('Exact DMD');
        % Build Atilde
        singular_values = diag(Sigmar);

        Atilde = transpose(Ur)*X2*Vr*diag(1./singular_values);

        assert( all(~isnan(Atilde),'all') && all(~isinf(Atilde),'all'), ...
            "Atilde shouldn't contain NaN or Inf - check invertibility of Sr")

        assert( length(Atilde) == rom_dim, "Requested model rank not achieved")

        % Compute DMD Modes
        [W, Lambda] = eig(Atilde);
        Phi = X2*Vr*diag(1./singular_values)*W;
        clear X2;
        %% Compute continuous-time eigenvalues
        lambda = diag(Lambda);
        optimalResiduals = nan(size(lambda));

    case 'rrr'
        disp('DDMD_RRR');

        % Algorithm 2 from:
        % Drmač, Zlatko, Igor Mezić, and Ryan Mohr. 2018.
        % “Data Driven Modal Decompositions: Analysis and Enhancements.”
        % SIAM Journal on Scientific Computing 40 (4): A2253–85. https://doi.org/10.1137/17M1144155.


        AUr = X2*(Vr*diag(1./sigma(1:subspaceSize)));
        clear X2;

        R = triu( qr( [Ur, AUr], 0 ) );
        numRankPrime = min( size(Ur,1) - subspaceSize, subspaceSize );
        Rblocks = mat2cell( R, [subspaceSize, size(Ur,1)-subspaceSize], [subspaceSize, subspaceSize] );
        [R11,~,R12,R22] = deal(Rblocks{:});

        % Rayleigh quotient
        Sk = diag( conj( diag(R11) ) )*R12;
        ritz = eig(Sk); % Ritz values - can be used as DMD values

        W = nan(subspaceSize);
        optimalResiduals = nan(subspaceSize,1);


        %%
        % Ritz value adjustment loop
        %
        ritzOld = zeros(size(ritz));
        count = 0;
        discrepancy = nan(options.ritzMaxIteration,1);
        ritzes = ritz;
        Mleft = [R12; R22];
        Mright = [R11; zeros(size(R22))];


        while not( allclose(ritzOld, ritz, options.ritzATOL, options.ritzRTOL) ) && ...
                count < options.ritzMaxIteration
            count = count+1;
            [~,discrepancy(count)] = allclose(ritzOld, ritz, options.ritzATOL, options.ritzRTOL);
            disp("Ritz correction try " + count + " Discrepancy = " + discrepancy(count) );
            ritzOld = ritz;
            for i = 1:subspaceSize

                M = Mleft - ritz(i)*Mright;

                % compute smallest singular value
                [~,SS,VV] = mysvd( M );
                optimalResiduals(i) = SS(end,end);

                svec_min = VV(:,end);
                ritz(i) = svec_min' * Sk * svec_min;
                W(:,i) = svec_min;
            end
            ritzes(:,end+1) = ritz;
        end
        clear Mleft Mright M
        %%
        out.discrepancy = discrepancy;
        out.ritzes = ritzes;

        Phi = Ur*W;
        lambda = ritzes(:,1);
end

Lambda = diag(lambda);
omega = log(lambda)/dt;

% normalize each mode by its 2-norm (division by a constant)
[Phi,~] = normalizeModes( Phi ); % see appendix

%% compute combination coefficients
% by L2-fit to a sequence of snapshots

disp("Computing b by L2 fit at snapshots " + ...
    num2str(options.step, "%d") );
disp( "time = " + ...
    num2str( (options.step-1)*dt, "%.2f ") );

% LHS: modes advanced to required steps
LHS = zeros([size(Phi), numel(options.step)]);

% advance mode vectors by steps and set them into layers of LHS tensor
for k = 1:numel(options.step)
    LHS(:,:,k) = Phi*Lambda^(options.step(k)-1);
end

% reshape LHS into a 2D matrix so that the layers are stacked on top of
% each other
LHS = permute(LHS, [1,3,2]);
LHS = reshape(LHS,[], size(Phi,2) );

% basically
% LHS = [ Phi * D^0 ; Phi * D^1 ; ...]


% RHS: corresponding snapshots
RHS = DataMatrix(:, options.step);
RHS = reshape(RHS,[], 1);

% basically
% RHS = [Datamatrix(:,1); Datamatrix(:,2);, ...]

% not always equivalent to LHS\RHS
% see: https://www.mathworks.com/help/matlab/ref/lsqminnorm.html#mw_e9cfa6e3-ccc3-4830-802c-df496b9e452b
b = lsqminnorm(LHS, RHS);


%% Compute mean L2 contribution of each mode
meanL2norm = abs(b) .* sqrt( (exp(2*real(omega)*T)-1)./(2*real(omega)*T) );
meanL2norm(isnan(meanL2norm)) = abs(b(isnan(meanL2norm))); % L2 norm of non-exponentials

% when real(omega) is zero, the formula above doesn't work as magnitude of
% the mode is constant - so we compute it manually
idx_small = abs(real(omega)*T) < 1e-12;
meanL2norm(idx_small) = abs(b(idx_small));

% catch any remaining errors
assert(all( ~isnan(meanL2norm) ), "Problems: L2 norm computation failed");


%% Sorting

switch(options.sortby)
    case 'initcond'
    [out.ranking,idx] = sort( abs(b), 'descend' );
    case 'l2'
    [out.ranking,idx] = sort( meanL2norm, 'descend' );
    case 'residual'
    [out.ranking,idx] = sort( optimalResiduals, 'ascend' );
    otherwise
        error('Unknown sorting option')
end


out.optimalResiduals = optimalResiduals(idx);
out.meanL2norm = meanL2norm(idx);
out.b = b(idx);
out.Phi = Phi(:,idx);
out.omega = omega(idx);
out.lambda = lambda(idx);
out.model_rank = length(lambda);
end

function [NormMatrix, Coefficients ] = normalizeModes( ModeMatrix )

NormalizationPhase = median( angle(ModeMatrix), 1 );
Coefficients = vecnorm( ModeMatrix ) .* exp( 1j * NormalizationPhase );
NormMatrix = ModeMatrix ./ Coefficients;
end

function [out,maxd] = allclose(a,b,atol,rtol)
% check that two arrays are within numerical tolerance to each other
d = abs(a(:)-b(:));
out = all( d - atol+rtol*abs(b(:)) <= 0 );
maxd= max(d - atol+rtol*abs(b(:)));

end
