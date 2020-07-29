function h = visualize_DMD_spectrum( omega, b, T, name)

halflife_cutoff = T/10;

sel = real(omega) > -log(2)/halflife_cutoff & imag(omega) >= 0;

meanL2norm = abs(b) .* sqrt( (exp(2*real(omega)*T)-1)./(2*real(omega)*T) );


h(1) = stem( imag(omega(sel)), abs(b(sel)), '-o',  'DisplayName', name + " |b|"); hold all;
h(2) = stem( imag(omega(sel)), meanL2norm(sel), 'x:', 'Color',h(1).Color,  'DisplayName', name + " mean L2"); hold off;

title({"Slowly decaying modes";"Halflife > " + halflife_cutoff} );
        xlabel("Im(\omega)")
