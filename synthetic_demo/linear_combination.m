% create spatiotemporal field by linear combination
clear
frequencies = [0,1.3j,-1.3j,2j,-2j,4.5j,-4.5j, -0.2+4.1j,-0.2-4.1j, -0.5+2.2j, -0.5-2.2j];
Nf = numel(frequencies);

N = 500;

x = linspace(0, Nf, N ).';


t = 0:0.25:50;
dt = t(2)-t(1);


Phi = [];
Lambda = [];

profile = @(z,z0)exp( -(z-z0).^2/0.5^2 );

for k = 1:numel(frequencies)
    
    Mag = 1/(1+abs(frequencies(k)));
    
    myprofile = triangularPulse( k-1, k, x);
    myprofile = profile(x, (k-1+k)/2 );
    
    Phi = [Phi, myprofile ];    
    Lambda = [Lambda; exp(frequencies(k)*t + 1j*(pi/2)*mod(k-1,2)) * Mag ];
    
end

Phi = Phi./vecnorm(Phi);

X0 = real( Phi * Lambda );


