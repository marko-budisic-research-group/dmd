function [ROM,out] = reduce_order( Phi, omega, b, t, varargin )
% REDUCE_ORDER Assemble reduced-order data containing only
% values indexed by idx.
%
% Usage:
% ROM = reduce_order( Phi, omega, b, t, idx )
%
%     Phi (:,:) double {mustBeNumeric}    % modal matrix
%     omega (:,1) double {mustBeNumeric}  % list of continuous-time frequencies
%     b (:,1) double {mustBeNumeric}      % list of combination coefficients
%     t (1,:) double {mustBeReal}         % time vector on which modes are evaluated
%     idx (:,1) = true(size(b))           % logical index vector (or list of indices) of modes to be included (default: all)
%
%
%
% The output is a sum of k-indexed terms
% ROM(:,j) = SUM[ Phi(:,k) * exp(omega(k) * t(j) ) * b(k) ]
%
% Note that there are no assumptions that conjugate modes are found -- if
% omega contains complex conjugates, for accurate reconstruction of real
% valued data, both conjugates need to be listed in idx

% arguments
%     Phi (:,:) double {mustBeNumeric}    % modal matrix
%     omega (:,1) double {mustBeNumeric}  % list of continuous-time frequencies
%     b (:,1) double {mustBeNumeric}      % list of combination coefficients
%     t (1,:) double {mustBeReal}         % time vector on which modes are evaluated
%     idx (:,1) = true(size(b))           % logical index of modes to be used (default: all)
% end

p = inputParser;

p.addRequired('Phi',@(M)validateattributes(M, 'numeric',{'2d','finite'})  );

p.addRequired('omega',@(s)validateattributes(s, 'numeric', {'column','finite'}) );
p.addRequired('b',@(s)validateattributes(s, 'numeric', {'column','finite'}) );
p.addRequired('t',@(s)validateattributes(s, 'numeric', {'row','real','finite'}) );
p.addOptional('idx',true(size(b)),@(s)validateattributes(s, 'numeric', {} ));

parse(p,Phi, omega, b, t, varargin{:});
options = p.Results;
idx = options.idx;


%% select only requested modes
out.model_rank = size(Phi, 2);
Phi = Phi(:,idx);
omega = omega(idx);
b = b(idx);
N = numel(b);

%% assemble the modes
for mm = 1:N

    ROM_single = ( Phi(:, mm)  * exp( omega(mm) * t ) ) * b(mm); % note that Phi * exp generates a matrix
    if mm == 1
        ROM = ROM_single;
    else
        ROM = ROM + ROM_single;
    end

end

out.Phi = Phi;
out.omega = omega;
out.b = b;
out.data_rank = size(Phi,2);
