function [dx,M] = rhs_synthetic(t,x,p)
%%RHS_SYNTHETIC

dx = nan(size(x));

%% LINEAR COMPONENTS
% oscillating p.Component
C1 = oscillating_matrix( complex( log(2)/p.Comp1.DoublingTime, 2*pi/p.Comp1.Period) );

% real decaying p.Component
C2 = log(2)/p.Comp2.DoublingTime;

% oscillating p.Component with repeated roots
C3Block = oscillating_matrix( complex( log(2)/p.Comp3.DoublingTime, 2*pi/p.Comp3.Period) );
C3 = [ C3Block, [0,0;1,0]; zeros(2), C3Block ]; 

% oscillating p.Component (this will end up being truncated to create higher
% harmonics)
C4 = oscillating_matrix( complex( log(2)/p.Comp4.DoublingTime, 2*pi/p.Comp4.Period) );

% linear parts of the Hopf and Pitchfork components
CHopf = oscillating_matrix( complex(p.hopf, 2*pi/p.hopf_period ) );
CPfork = p.pitchfork;


M = blkdiag(C1, C2, C3, C4,CHopf, CPfork);

dx(1:12,:) = M * x(1:12,:);

%% NONLINEAR COMPONENTS

% Hopf nonlinearity 
xh = x(10,:);
yh = x(11,:);

dx(10,:) = dx(10,:) - xh.*( xh.^2 + yh.^2 );
dx(11,:) = dx(11,:) - yh.*( xh.^2 + yh.^2 );

% Pitchfork nonlinearity 
xp = x(12,:);
dx(12,:) = dx(12,:) - xp.^3;


end

% Auxiliary function returning a matrix with eigenvalue passed in the
% argument and its conjuage
function M = oscillating_matrix( z )


M = [ real(z), imag(z); -imag(z), real(z)  ] ;
%M(isinf(M)) = 0;

end