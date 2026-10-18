function [Baseline2, Final_specs] = BaselineCorrection(fids, F0, BW, LW, SNR, debugFLAG)
%% Load data
if nargin == 0
    debugFLAG = true;

    % Select file
    [file,path] = uigetfile({'*.mat','Matlab Files'},'File Selector','MultiSelect','on');

    % Obtain full name
    F = fullfile(path,file);
    T = load(F);
    if isfield(T,'Starting')
        F0 = T.Starting.F0/(10^6);  % Convert F0 to MHz
        BW = T.Starting.BW;

        % Ensure data is of proper shape, squeeze needs transposing (5120 data_pts x 1 acq to 1 acq x 5120 data_pts)
        if size(T.Starting.Data,2) == 1
            fids = squeeze(T.Starting.Data(1,:,:)).';
        else
            fids = squeeze(T.Starting.Data(1,:,:));
        end

    elseif isfield(T,'Processed')
        BW = T.Processed.BW_Hz;
        F0 = T.Processed.F0_MHz;
        LW = T.Processed.Linewidth_Hz;
        SNR = T.Processed.SNR;

        fids = [];
        fids(:,:) = T.Processed.Time_domain;
    end
end

% Combine averages (if necessary)
if size(fids,1) > 1
    fids = mean(fids,1);
end

% Convert to spectrum
f_comp = fftshift(fft(fids));

%% Perform baseline correction
% Find peaks
[TruePeaks, filtsz] = FindPeaks(f_comp, F0, LW, BW, SNR, debugFLAG);

% Calculate SNR
SNR = max(real(f_comp))/std(real(f_comp(end-length(f_comp)/8:end)));

% Create baseline without peaks
% Slightly filter data if SNR is too low
if SNR < 50
    padval = length(f_comp)/2;
    spec_corr = padarray(f_comp.', padval).';

    % Extract corrected fids
    fids_corr = ifft(ifftshift(spec_corr));

    % Temporarily filter (and flip) data
    spec_corr = fftshift(fft(exp(-(pi*filtsz/2).*((1:length(fids_corr))/(BW*2))).*fids_corr)).';

    % Remove zero padding
    Baseline0 = real(spec_corr(padval+1:end-padval));

    clear padval spec_corr fids_corr
else
    Baseline0 = real(f_comp).';
end

% Linear interpolation
Baseline1 = RemovePeak(TruePeaks, Baseline0);

% Create structure for model
if exist('TruePeaks','var'); clear TruePeaks; end
for i = 1:length(TruePeaks)
    TruePeaks(i).Start = TruePeaks(i).Start(1);
    TruePeaks(i).End = TruePeaks(i).Stop(1);
end

% Model baseline
Baseline2 = New_InterativeModel(Baseline1, 1000, 0, Baseline0, TruePeaks);

% Obtain corrected spectra
Final_specs = complex(real(f_comp).' - Baseline2, imag(f_comp).').';

if debugFLAG
    plot(Baseline2.'); 
    % figure; plot(real(f_comp)); axis tight; hold on; plot(Baseline1.');  plot(Baseline2.'); legend('Raw Data', 'Interpolated Baseline',  'Smoothed Baseline'); %#ok<UNRCH>
    % figure; plot(real(Final_specs)); axis tight; title('Corrected Spectra')
    % figure; plot(real(f_comp-Final_specs)); axis tight; title('Difference'); ylim([min(real(Final_specs)), max(real(Final_specs))])
end

%% Functions

    function Baseline = RemovePeak(Peaks, f)
        % Initialize baseline data
        Baseline = f;
        
        % Check for overlapping peaks
        Keeps = true(length(Peaks),1);
        for k = 2:length(Peaks)
            if Peaks(k-1).Stop(1) > Peaks(k).Start(1)
                % Combine peaks
                Peaks(k).Start(1) = Peaks(k-1).Start(1);
        
                % Flag old peak
                Keeps(k-1) = false;
            end
        end
        
        % Remove old peaks
        Peaks = Peaks(Keeps);
        
        for k = 1:length(Peaks)
            % Extract peak start & end from structure
            PKstart = Peaks(k).Start(1);
            PKend = Peaks(k).Stop(1);
        
            % Estimate slope across peak range
            PKslope = (f(PKend+1) - f(PKstart-1))/(PKend - PKstart + 2);
        
            % Check slope
            if PKslope == 0
                % Replace peak with flat line
                Baseline(PKstart:PKend) = f(PKend+1)*ones(PKend-PKstart+1, 1);
            else
                % Replace peak with linearly interpolated data
                kk = f(PKstart-1):PKslope:f(PKend+1);
                Baseline(PKstart:PKend) = kk(2:length(kk)-1);
            end
        end
    end
end