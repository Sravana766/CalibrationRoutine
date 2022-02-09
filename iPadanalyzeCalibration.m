clc;
clear;
%% Read calibration data measured by calibration.py script
%
% NOTE: throughout the entire script, channel 1 is the right ear and
% channel 2 is the left ear. matlab's audioplayer is the other way around

% clear;
% close all;

dataPath = 'data';
% plotPath = '';
plotPath = 'plots';
impzPath = 'impz';

files = {
%     'calibration-20210212-170022-10.0dbspl-100.mat'
%     'calibration-20210212-171945-15.0dbspl-100.mat'
%     'calibration-20210212-175300-20.0dbspl-100.mat'
%     'calibration-20210212-181117-25.0dbspl-100.mat'
%     'calibration-20210212-183819-30.0dbspl-100.mat'
%     'calibration-20210212-190054-35.0dbspl-100.mat'
%     'calibration-20210212-192656-40.0dbspl-100.mat'
%     'calibration-20210212-221304-45.0dbspl-2ndtry-100.mat'
%     'calibration-20210212-163557-50.0dbspl-40.mat'
%     'calibration-20210212-201648-55.0dbspl-50.mat'
%     'calibration-20210212-202938-60.0dbspl-50.mat'
%     'calibration-20210212-204137-65.0dbspl-50.mat'
%     'calibration-20210212-210021-70.0dbspl-100.mat'
%     'calibration-20210212-212322-75.0dbspl-100.mat'
    'calibration-20210212-215603-80.0dbspl-100.mat'
    
    
};

cal = struct();    % at the end store this struct in cal.mat

for fileID = 1:length(files)
	%% clear left over vars from previous iterations
	clearvars -except ...
		fileID files dataPath plotPath impzPath calName cal ...
		IMPZ0R_LEVEL_QUIET IMPZ0L_LEVEL_QUIET ...
		IMPZ0R_PHASE_QUIET IMPZ0L_PHASE_QUIET;
	
	%% load calibration measurement file
	if iscell(files{fileID})
		file = files{fileID}{1};
		level = files{fileID}{2};
	else
		file = files{fileID};
	end
	
	load(strjoin({dataPath,file}, '/'));

	if ~exist('level','var')    % 'level' is shown in plots
		if exist('stimSPL','var')
			if stimSPL <= 0; level = 'Quiet';
			else; level = [num2str(stimSPL) ' dB SPL'];
			end
		elseif exist('stimSL' ,'var')
			if stimSL == -100; level = 'Quiet';
			else; level = [num2str(stimSL ) ' dB SL' ];
			end
		else; level = '?';
		end
	end
	
	if ~exist('gapDur','var'); gapDur = 0; end    % compatibility
	
	fprintf('Analyzing calibration at %s\n', level);


	%% load MLS code
	load('mls18.mat');
	mls = mls18 - mean(mls18);
	mls = mls./max(abs(mls));
% 	mls = filtfilt(b,a,mls);
	% FFT of the MLS sequence:
	MLS = fft(-(2*(mls)-1));

	
	%% setup parameters
	f1 = [100, 25e3];  % min and max frequency in plots (Hz)
	f2 = [100, 12000];  % min and max frequency for dB calulations
	fs = inpFs;
	n  = 100;   % movmean
	[filt_b,filt_a] = butter(4, [20,20e3]/(fs/2),'bandpass');
	
	f = fs/length(mls)*(0:length(mls)-1);
	f1Range = f1(1)<=f&f<=f1(2);
	f2Range = f2(1)<=f&f<=f2(2);
	
	
	%% calculate impulse response
	if exist('inpAvgTrim', 'var')
		impz0 = zeros(2, length(mls));
		for j = 1:2
			sig = inpAvgTrim(j,1:length(mls));
			%sig = inpSumTrim(j,1:length(mls)) / double(countTrim);
			sig = sig - mean(sig);
			SIG = fft(sig);
			IMPZ0 = SIG.*conj(MLS);
			impz0(j,:) = real(ifft(IMPZ0));
		end
		impz0r = impz0(1,:);
		impz0l  = impz0(2,:);
		
	elseif exist('inpAll', 'var')
		% impulse response for each repeated measurement
		impz0 = zeros(count, 2, length(mls));
% 		fig2 = figure; hold on;
		for i = 1:count
			for j = 1:2
				sig = inpAll(i, j, :);
				sig = squeeze(sig)';
				sig = filtfilt(filt_b, filt_a, sig);
				sig = sig - mean(sig);
				
				% do not sync for quiet recording
				if strcmp(level, 'Quiet')
					indTrim = .5*fs;
				% sync relative to the first (right) channel
				elseif j==1
					look = round(syncDur + gapDur + 2*fs);
					xcor = conv(sig(1:look), fliplr(sync), 'same');
					[~, ind] = max(xcor);
% 					[~, ind] = max(left(1:look));
					indSync = ind - round(fs*syncDur/2);
					indTrim = ind + floor(fs*syncDur/2) + floor(fs*gapDur);
					if indSync < 1; continue; end
				% plot all recordings
% 				elseif j==2
% 					if indSync < 1; continue; end
% 					sigSync = sig(indSync:indSync+round(length(sig)/10));
% 					plot((0:length(sigSync)-1)/Fs, sigSync);
% 					sigTrim = sig(indTrim:indTrim+round(length(sig)/10));
% 					plot((0:length(sigTrim)-1)/Fs, sigTrim);
				elseif j==2
					if indSync < 1; continue; end
				end
				sig = sig(indTrim:end);
				sig = sig(1:length(mls));
				sig = sig - mean(sig);

				SIG = fft(sig);
				IMPZ0 = SIG.*conj(MLS);
				impz0(i,j,:) = real(ifft(IMPZ0));
			end
		end
		
		% average all impulse responses
		impz0_mean = zeros(2, length(mls));
		for j = 1:2    % per each right/left channel
			impz0_temp = squeeze(impz0(:,j,:));
			% drop outliers based on RMS z-score
			if ~strcmp(level, 'Quiet')
				z = zscore(rms(impz0_temp,2));
				msk = abs(z)<1.5;
				fprintf('Dropping %d recordings from channel %d\n', ...
					sum(~msk), j);
				impz0_temp = impz0_temp(msk,:);
			end
			impz0_mean(j,:) = mean(impz0_temp,1);
		end
		% separate right and left channel reponse
		impz0r = impz0_mean(1,:);
		impz0l  = impz0_mean(2,:);
% 		return;
	else
		error('Error');
	end
	
	
	%% setup figure	
	fig1 = figure;
	w = 1000; h = 900;
	set(gcf, 'position', [(1920-w)/2 (1080-h)/2 w h]);


	%% plot left
	IMPZ0L_LEVEL = 20*log10(abs(fft(impz0l)));
	IMPZ0L_PHASE = angle(fft(impz0l));
	
% 	IMPZ0L_LEVEL = IMPZ0L_LEVEL - mean(IMPZ0L_LEVEL(ind));
% 	IMPZ0L_LEVEL = IMPZ0L_LEVEL - IMPZ0L_LEVEL_QUIET;

	% Gentle boxcar smoothing over frequency domain
% 	b = [(0:N-1)/N ones(1,2*N) (N1-1:-1:0)/N];
% 	IMPZ0L_LEVEL = filter(b,1,IMPZ0L_LEVEL);
% 	IMPZ0L_LEVEL = IMPZ0L_LEVEL(2*N-1:end-2*N)/(4*N);

	%plot the final speaker transfer function
	ax1 = subplot(221);
    grid on
	plot(f(f1Range), IMPZ0L_LEVEL(f1Range));
	hold on;
	plot(f(f1Range), movmean(IMPZ0L_LEVEL(f1Range), n), ...
		'r', 'linewidth',2);
    grid on
	xlabel('Frequency [Hz]')
	ylabel('Magnitude [dB]')
	set(gca,'xscale','log','xtick',[500 1e3 2e3 4e3 8e3 16e3],...
		'xticklabel',{'500';'1k';'2k';'4k';'8k';'16k'})
	xlim([0.95*f1(1) 1.1*f1(2)])
	axis square tight;
	title({'Left', sprintf('Mean (%g-%g Hz): %.6g dB', ...
		f2, mean(IMPZ0L_LEVEL(f2Range)))});

	
	%% plot right
	IMPZ0R_LEVEL = 20*log10(abs(fft(impz0r)));
	IMPZ0R_PHASE = angle(fft(impz0r));
	
% 	IMPZ0R_LEVEL = IMPZ0R_LEVEL - mean(IMPZ0_RIGHT_LEVEL(ind));
% 	IMPZ0R_LEVEL = IMPZ0R_LEVEL - IMPZ0_RIGHT_LEVEL_QUIET;

	% Gentle boxcar smoothing over frequency domain
% 	b = [(0:n-1)/n ones(1,2*n) (n-1:-1:0)/n];
% 	IMPZ0R_LEVEL = filter(b,1,IMPZ0R_LEVEL);
% 	IMPZ0R_LEVEL = IMPZ0R_LEVEL(2*n-1:end-2*n)/(4*n);

	%plot the final speaker transfer function
	ax2 = subplot(222);
    
	plot(f(f1Range), IMPZ0R_LEVEL(f1Range));
	hold on;
	plot(f(f1Range), movmean(IMPZ0R_LEVEL(f1Range),n), ...
		'r', 'linewidth',2);
    grid on
	xlabel('Frequency [Hz]')
	ylabel('Magnitude [dB]')
	set(gca,'xscale','log','xtick',[500 1e3 2e3 4e3 8e3 16e3],...
		'xticklabel',{'500';'1k';'2k';'4k';'8k';'16k'})
	xlim([0.95*f1(1) 1.1*f1(2)])
	axis square tight;
	linkaxes([ax1 ax2], 'xy');
	title({'Right', sprintf('Mean (%g-%g Hz): %.6g dB', ...
		f2, mean(IMPZ0R_LEVEL(f2Range)))});

	
	%% plot ILD
	subplot(223); cla;
	% fig2 = figure;
	ILD = IMPZ0R_LEVEL - IMPZ0L_LEVEL;
	plot(f(f2Range), ILD(f2Range));
	hold on;
	plot(f(f2Range), movmean(ILD(f2Range),n), 'r', 'linewidth',2);
	xlabel('Frequency [Hz]')
	ylabel('Magnitude [dB]')
% 	set(gca,'xscale','lin','xtick',[500 1e3 2e3 4e3 8e3 16e3],...
% 		'xticklabel',{'500';'1k';'2k';'4k';'8k';'16k'})
	set(gca,'xscale','lin','xtick',[300, 600, 900, 1200]);
	axis square tight;
	xlim(f2);
% 	ylim([-1,+1]*std(ILD1(f2Range))+mean(ILD1(f2Range)));
	title({'ILD (R - L)', sprintf('Mean: %6g dB, SEM: %6g', ...
		mean(ILD(f2Range)), std(ILD(f2Range))/(sqrt(sum(f2Range))-1))});

	
	%% plot IPD
	figure(fig1);
	subplot(224); cla;
	IPD = IMPZ0L_PHASE - IMPZ0R_PHASE;
	IPD(IPD>pi) = IPD(IPD>pi) - 2*pi;
	IPD(IPD<-pi) = IPD(IPD<-pi) + 2*pi;
	plot(f(f2Range), IPD(f2Range));
	hold on;
	plot(f(f2Range), movmean(IPD(f2Range),n), 'r', 'linewidth',2);
	xlabel('Frequency [Hz]')
	ylabel('Phase [rad]')
% 	set(gca,'xscale','lin','xtick',[500 1e3 2e3 4e3 8e3 16e3],...
% 		'xticklabel',{'500';'1k';'2k';'4k';'8k';'16k'})
	set(gca,'xscale','lin','xtick',[300, 600, 900, 1200]);
	axis square tight;
	xlim(f2);
% 	ylim([-1,+1]*std(IPD(f2Range))+mean(IPD(f2Range)));
	title({'IPD (R - L)', ['Mean: ' num2str(mean(IPD(f2Range))) ' rad']});

	mtitle = ['Stimulus level: ' level ', Repeats: ' num2str(count)];
	axes('Position',[0 0 1 1],'Visible','off');
	text(.5,.98, mtitle, 'horizontalalignment', 'center', ...
		'fontweight', 'bold', 'fontsize',12, 'interpreter','none');
	uistack(gca, 'bottom');
	

	%% save
	if ~isempty(plotPath)
		fileName = files{fileID}(1:end-4);
		saveas(fig1, strjoin({plotPath,[fileName '.png']}, '/'));
		saveas(fig1, strjoin({plotPath,[fileName '.fig']}, '/'));
	end
	
	cal(fileID).spl = stimSPL;
	cal(fileID).fs = fs;
	cal(fileID).freqs = f;
	% channel 1 is right, channel 2 is left
	cal(fileID).impz0(1,:) = impz0r;
	cal(fileID).impz0(2,:) = impz0l;
	cal(fileID).IMPZ0_LEVEL(1,:) = movmean(IMPZ0R_LEVEL,n);
	cal(fileID).IMPZ0_LEVEL(2,:) = movmean(IMPZ0L_LEVEL,n);
	cal(fileID).IMPZ0_PHASE(1,:) = movmean(IMPZ0R_PHASE,n);
	cal(fileID).IMPZ0_PHASE(2,:) = movmean(IMPZ0L_PHASE,n);
	cal(fileID).ILD = ILD;
	cal(fileID).IPD = IPD;
end

save(strjoin({impzPath, ['impz-' '.mat']}, '/'), 'cal');

% IMPZ0R_LEVEL_QUIET = IMPZ0R_LEVEL;
% IMPZ0L_LEVEL_QUIET = IMPZ0L_LEVEL;
% IMPZ0R_PHASE_QUIET = IMPZ0R_PHASE;
% IMPZ0L_PHASE_QUIET = IMPZ0L_PHASE;