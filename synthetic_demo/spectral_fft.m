function [power, phase, F] = spectral_fft( x, dt, removemean )
%SPECTRAL_FFT Compute the one-sided power spectrum and phase spectrum of the
%scalar signal sampled at time step dt.
% function [power, phase, F] = spectral_fft( x, dt, removemean )
% x -- time series
% dt -- sampling step
% removemean - (optional; default:false) - remove mean before computing FFT

validateattributes(x, {'numeric'}, {'vector'});


if nargin==3 && removemean
    x = x - mean(x);
end

Nx = numel(x);

% compute fft and reassemble it in format that centers 0 component
X = fftshift( fft(x, 2*floor(Nx/2)+1) );
NX = numel(X);

F = ((-floor(NX/2)):floor(NX/2))/(dt*NX);

% compute the power (mean modulus squared) and phase of the FFT
power = (X .* conj(X))/NX;
phase = angle(X);

% keep values only positive frequencies that are smaller
% than the Nyquist frequency (1/2 dt)
sel = F >= 0 & F < 1/(2*dt);
F = F(sel);
power = power(sel);
phase = phase(sel);

%% plotting if no output is requested
if nargout == 0
  a1 = subplot(2,2,1);
  plot(F,power);
  xlabel('Frequency [cycles/time unit]');
  ylabel('Power');

  a2 = subplot(2,2,2);
  plot(F,phase/pi);
  xlabel('Frequency [cycles/time unit]');
  ylabel('Phase [fractions of  \pi ]');

  linkaxes([a1,a2],'x');

  a3 = subplot(2,2,3);
  plot(1./F,power);
  xlabel('Period [time unit]');
  ylabel('Power');

  a4 = subplot(2,2,4);
  plot(1./F,phase/pi);
  xlabel('Period [time]');
  ylabel('Phase [fractions of \pi ]');

  linkaxes([a3,a4],'x');

end
