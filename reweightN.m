% Reweight the magnitude spectrum of signal with no phase change
% Arguments:
%   x: signal to be reweighted
%   fs: signal sampling frequency
%   freqs: in Hz
%   levels: in dB
% Returns:
%   x2: reweighted signal
%   magFilt: magnitude of the applied filter, at f = (0:n-1)*fs/length(x)
function varargout = reweightN(x, fs, freqs, levels)
	n = length(x);
	f = (0:n-1)*fs/n;

	% Calculate FFT of the signal and decompose its magnitude and phase
	X = fft(x);
	Xmag = abs(X);
	Xang = angle(X);

	% Generate the inverse A-weighting filter response magnitude
	range = freqs(1)<=f & f<=freqs(end);		                           % Boolean array
	levels2 = interp1(freqs, levels, f(range), 'pchip');                   % Interpolate within the above range
	mags = 10.^(levels2/20);                                               % Convert dB levels to scale factors
	magFilt = ones(1,n);                                                   % Initialize the actual magnitude filter
	magFilt(range) = mags;                                                 % Set the first half of the filter magnitude
	magFilt(fliplr([range(2:end) false])) = fliplr(mags);                  % Set the second half of the filter magnitude

	% Filter the signal
	X2mag = Xmag .* magFilt;                                               % Filter the signal's magnitude
% 	meanX2mag = mean(X2mag);
% 	newX2 = X2mag - meanX2mag;
	X2 = X2mag .* exp(1i*Xang);                                            % Contruct the FFT of the filtered signal
	x2 = real(ifft(X2));                                                   % Calculate the actual filtered signal

	f1Range =find(f>50&f<12e3);
% 	ori = X2mag(f1Range);
% 	hello = newX2(f1Range);
% 	inverse = -1.*newX2(f1Range);
	
	%Plot the functions
    figure
     grid on
	plot(f(f1Range),20*log10(X2mag(f1Range)));
	hold on;
	plot(f(f1Range), 20*log10( movmean((X2mag(f1Range)), 100)), ...
		'k', 'linewidth',2);
	
	% Output according to the requested output arguments
	if nargout == 1
		varargout = {x2};                                                  % Only output the filtered signal
	elseif nargout == 2
		varargout = {x2, magFilt};                                         % Output both the filtered signal and the filter magnitude
	else
		error('Wrong number of output arguments');
	end
end
