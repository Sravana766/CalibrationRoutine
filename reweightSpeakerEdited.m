% Inverse filter speaker transfer function
% 
% Arguments:
% - x: signal to be reweighted, row 1 is left channel and row 2 is right
% - fs: sampling frequency (Hz)
% - spl: the level at which signal should be compensated for (dB SPL)
% - fl: the lowest frequency to apply reweighting at (Hz)
% - fh: the highest frequency to apply reweighting at (Hz)
%
% Returns:
% - x2: reweighted signal
% - magFilt: magnitude of the applied filter, at f = (0:n-1)*fs/length(x)
%
% Use case:
% - reweightSpeaker('load', 'emotiva1-back'): load calibration data
function varargout = reweightSpeakerEdited(x, fs, spl, fl, fh)
	persistent cal;
	
	if nargin == 1 && strcmpi(x, 'load')
		% load latest version of calibration data for the given device
		file = 'impz80.mat';
% 		if isempty(file)
% 			error(['[reweightSpeaker] no calibration data']);
% 		end
% 		names = {files};
% 		res = regexp(names,'^impz-','tokens');
% 		vers = cellfun(@(c)str2double(c{1}{1}), res);
% 		[~, inds] = sort(vers);
% 		file = names{inds(end)};
		load(file, 'cal');
		return;
	end
	
	if strcmpi(cal, 'none')
		varargout = {x};
		return;
	end
	
	if isempty(cal)
		error('[reweightSpeaker] no calibration data');
	end
	
	spls = [cal(:).spl];
	[~, i] = min(abs(spls-spl));    % find the closest spl
	freqs = cal(i).freqs;
	range = fl<=freqs & freqs<=fh;
	freqs = freqs(range);
	levels = cal(i).IMPZ0_LEVEL(:,range);
	% NOTE: in calibration code channel 1 is right and channel 2 is left
	% but matlab's audioplayer is the other way around
	% flip channels: from now on channel 1 is left and channel 2 is right
	levels = flip(levels,1);
	avg = mean(levels(:));
	levels = avg - levels;
	% zero pad at an octave above and below
	freqs = [min(freqs)/2, freqs, max(freqs)*2];
	zpad = zeros(size(levels,1),1);
	levels = [zpad, levels, zpad];
	
	x2 = [];
	magFilt = [];
	for i = 1:size(x,1)
		[x2(i,:), magFilt(i,:)] = reweightN(x(i,:), fs, freqs, levels(i,:));
	end
	
	if nargout == 1
		varargout = {x2};
	elseif nargout == 2
		varargout = {x2, magFilt};
	else
		error('Wrong number of output arguments');
	end	
end
