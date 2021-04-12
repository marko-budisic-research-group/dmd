 function out = dmd( DataMatrix, dt, r, options )
%% DMD Dynamic Mode Decomposition (exact version)
%
% out = dmd( DataMatrix, dt, r, ... )
%
% Inputs:
%     DataMatrix (:,:) double {mustBeNumeric, mustBeReal, mustBeFinite} -- INPUT SNAPSHOT MATRIX (columns are snapshots)
%     dt (1,1) double {mustBePositive, mustBeFinite} -- INPUT SNAPSHOT TIMESTEP (columns are snapshots)
%     r (1,1) double {mustBePositive,mustBeInteger,mustBeFinite} -- REQUESTED DIMENSION OF REDUCED-ORDER MODEL
%
% Additional optional name-value pairs, specified as
% out = dmd( DataMatrix, dt, r, NAME, VALUE, NAME, VALUE,... )
%
%     out = dmd( DataMatrix, dt, r, 'rom_type', VALUE)
%         VALUE = 'lsq' - standard least-squares truncation {default}
%         VALUE = 'tlsq' - total least squares truncation (Hemati 2017)
%
%     out = dmd( DataMatrix, dt, r, 'rrr', VALUE)
 %        Drmac et al. DDMD_RRR optimization
%         VALUE = 'none' Exact DMD style eigenvalue calculation
%         VALUE = 'compute' - Compute residual errors for Az = \lambda z
%                             eigenvalue computation (stored in out.residual)
%         VALUE = 'optimize' - Optimize residual error for Az = \lambda z
%                             eigenvalue computation (stored in
%                             out.residual)
%
%    out = dmd( DataMatrix, dt, r, 'sortbyb', VALUE)
%         VALUE = true - sort DMD modes by |b|
%         VALUE = false - sort DMD modes by average L2 norm of time
%                          coefficient {default}
%
%    out = dmd( DataMatrix, dt, r, 'sortbyb', VALUE)
%         VALUE = true - sort DMD modes by |b|
%         VALUE = false - sort DMD modes by average L2 norm of time
%                          coefficient {default}
%
%    out = dmd( DataMatrix, dt, r, 'step', VALUE)
%         VALUE is a row-vector of snapshot indices at which optimization is
%         performed to determine b coefficients in the expansion.
%         If VALUE = -1, all snapshots are used.
%
%    out = dmd( DataMatrix, dt, r, 'removefrequencies', VALUE)
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
%
%
% Outputs:
% out.meanL2norm - time-average of L2 norm of mode
% out.b - contribution of the mode to 0th timestep
% out.Phi - DMD mode vector
% out.omega - continuous time DMD eigenvalue omega = log( lambda ) / dt
% out.lambda - discrete time DMD eigenvalue lambda = exp( omega * dt )
% out.model_rank - rank of the model (= input r parameter)

arguments

    DataMatrix (:,:) double {mustBeNumeric, mustBeReal, mustBeFinite}
    dt (1,1) double {mustBePositive, mustBeFinite}
    r (1,1) double {mustBePositive,mustBeInteger,mustBeFinite}
    options.rom_type {mustBeMember(options.rom_type,{'lsq','tlsq'})} = ...
        'lsq'
    options.rrr {mustBeMember(options.rrr,{'none','compute','optimize'})} = ...
        'none'
    options.removefrequencies (1,:) double = []
    options.sortbyb (1,1) logical = false
    options.step (1,:) double {mustBeInteger,mustBeFinite} = 1
end

if options.step == -1
    options.step = 1:size(DataMatrix,2);
end

assert( r<= min(size(DataMatrix)), ...
    "ROM order has to be smaller than rank of input")

%% Preprocessing
Nsnapshots = size(DataMatrix,2);
T = (Nsnapshots-1)*dt;

X1=DataMatrix(:,1:end-1); % likely a very slow part of the code
X2=DataMatrix(:,2:end);


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

    meanL2 = abs(out.AvgB) .* sqrt( (exp(2*real(out.AvgOmega)*T)-1)./(2*real(out.AvgOmega)*T) );
    meanL2(isnan(meanL2)) = abs(out.AvgB(isnan(meanL2)));

    out.AvgMeanL2norm = meanL2;

end


% total-DMD Hemati debias
switch(options.rom_type)
    case 'tlsq'
        disp("Total least-squares truncation to order " + r)
        Z = [X1;X2];
        [~,~,Q] = svd(Z,'econ');
        Qr = Q(:,1:r);
        X1 = X1*Qr;
        X2 = X2*Qr;
    case 'lsq'
        disp("Standard least-squares truncation to order " + r)
end


%% Core DMD algorithm

% SVD
[U, S, V] = svd(X1, 'econ');
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:,1:r);

% Build Atilde
singular_values = diag(Sr);

Atilde = transpose(Ur)*X2*Vr*diag(1./singular_values);

assert( all(~isnan(Atilde),'all') && all(~isinf(Atilde),'all'), ...
    "Atilde shouldn't contain NaN or Inf - check invertibility of Sr")

assert( length(Atilde) == r, "Requested model rank not achieved")

% Compute DMD Modes
[W, D] = eig(Atilde);
Phi = X2*Vr*diag(1./singular_values)*W;

% IMPLEMENT FROM:
% Drmač, Zlatko, Igor Mezić, and Ryan Mohr. 2018.
% “Data Driven Modal Decompositions: Analysis and Enhancements.”
% SIAM Journal on Scientific Computing 40 (4): A2253–85. https://doi.org/10.1137/17M1144155.

switch (options.rrr)
  case 'compute':
    disp("PLACEHOLDER FOR COMPUTING RESIDUALS DMD_RRR")
  case 'optimize':
    disp("PLACEHOLDER FOR OPTIMIZING RESIDUALS DMD_RRR")
end

% normalize each column by its 2-norm (division by a constant)
[Phi,~] = normalizeModes( Phi ); % see appendix

%% Compute continuous-time eigenvalues
lambda = diag(D);
omega = log(lambda)/dt;

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
    LHS(:,:,k) = Phi*D^(options.step(k)-1);
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
if options.sortbyb
    [~,idx] = sort( abs(b), 'descend' );
else
    [~,idx] = sort( meanL2norm, 'descend' );
end

out.meanL2norm = meanL2norm(idx);
out.b = b(idx);
out.Phi = Phi(:,idx);
out.omega = omega(idx);
out.lambda = lambda(idx);
out.model_rank = length(Atilde);
end

function [NormMatrix, Coefficients ] = normalizeModes( ModeMatrix )

NormalizationPhase = median( angle(ModeMatrix), 1 );


Coefficients = vecnorm( ModeMatrix ) .* exp( 1j * NormalizationPhase );
NormMatrix = ModeMatrix ./ Coefficients;
end
