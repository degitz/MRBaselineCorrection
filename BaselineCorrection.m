function [Baseline2, Final_specs] = BaselineCorrection(fids, F0, BW, Nucleus, LW, pltFLAG, ID)
%% Load data
if nargin == 0
    pltFLAG = true;

    % Select file
    [file,path] = uigetfile({'*.mat','Matlab Files';'*.txt','Text Files'},'File Selector','MultiSelect','on');

    % Obtain full name
    F = fullfile(path,file);
    T = load(F);
    if isfield(T,'Starting')
        F0 = T.Starting.F0/(10^6);  % Convert F0 to MHz
        BW = T.Starting.BW;
        Nucleus = T.Starting.Nucleus;
        ID = char(T.Starting.Filename);

        % Ensure data is of proper shape, squeeze needs transposing (5120 data_pts x 1 acq to 1 acq x 5120 data_pts)
        if size(T.Starting.Data,2) == 1
            fids = squeeze(T.Starting.Data(1,:,:)).';
        else
            fids = squeeze(T.Starting.Data(1,:,:));
        end

    elseif isfield(T,'Processed')
        BW = T.Processed.BW_Hz;
        F0 = T.Processed.F0_MHz;
        Nucleus = T.Processed.Nucleus;
        LW = T.Processed.Linewidth_Hz;
        ID = char(T.Processed.Source_files);

        fids = [];
        fids(:,:) = T.Processed.Time_domain;
    end

    % Extract ID
    if strcmpi(ID(1:3),'MID')
        IDtmp = strsplit(ID, '_');
        ID = char(IDtmp{1});
        clear IDtmp
    else
        ID = '';
    end
end

% Calculate x-axis
Faxis = linspace(-BW/2, BW/2,size(fids,2))/F0;
if strcmp(Nucleus, '1H')
    Faxis = Faxis + 4.65;
end

% Determine x-axis step size in terms of number of points
Faxis_ss = (BW/length(Faxis))/F0; % ppm/pt

% Calculate number of scales for differentiation
scales = LW/(F0*Faxis_ss*8);
if scales < 3, scales = 3; end

% Combine averages (if necessary)
if size(fids,1) > 1
    fids = mean(fids,1);
end

% Convert to spectrum
f_comp = fftshift(fft(fids));

% Plot data
if pltFLAG
    pdata = [real(f_comp); imag(f_comp); abs(f_comp); angle(f_comp)]; 
    ptitle = ["Real", "Imaginary", "Magnitude", "Phase"];
    figure; tiledlayout(3,1,"TileSpacing","compact","Padding","tight");
    for i = [1, 3, 4]
        f = squeeze(pdata(i,:)); 
        [fSTD, df, dfAD, dfSTD, ~] = PeakStats(f, scales);

        nexttile;
        plot(f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on
        yline(fSTD, 'Color', [0 .447 .741], 'LineStyle', '--');
        plot(df, 'Color', [.635 .078 .184], 'LineWidth', 2);
        plot(dfAD, 'Color', [.635 .078 .184]);
        yline(dfSTD, 'Color', [.635 .078 .184], 'LineStyle', '--');
        ylabel(ptitle(i));
        if i == 1; title(ID); legend('f', 'fSTD', 'df', 'dfAD', 'dfSTD'); end; if i == 4; ylim([ -2*pi 2*pi]); end
        clear f fAD fSTD df dfAD dfSTD
    end
    f_og = real(f_comp); %#ok<NASGU>
    clear ptitle pdata i
end

% Filter data
[F, filtsz] = FilterSpec(f_comp, BW, 'filter', 'Complex');
f_mg = abs(F);
f = real(F);

% Calculate SNR
SNR = max(real(f_comp))/std(real(f_comp(end-length(f_comp)/8:end)));

% Determine linewidth of main peak (in Hz)
LW2 = Linewidth_calculator(f, F0, BW);

% Determine sample size in terms of number of points (we will use LW as sample size, can be changed later)
ss = ceil((LW2/F0)/Faxis_ss);

% Calculate statistical datasets
[fSTD, df, dfAD, dfSTD, dfsigns] = PeakStats(f, scales);

ss = ceil(ss);

% Obtain first standard deviation line for phase
% df_ph = diff(angle(f_comp));
df_ph = cwt_custom(angle(f_comp), scales);
dfSTD_ph = mean(df_ph) + std(df_ph);

%% Find easy magnitude peaks, min of 8, max of 30
% Obtain df & dfSTD via numerical differentiation
df_tmp = diff(f); df_nopeaks = df_tmp;
df_nopeaks(df_tmp > std(df_tmp)) = NaN; df_nopeaks(df_tmp < -std(df_tmp)) = NaN;
df_nopeaks(df_tmp > std(df_nopeaks)) = NaN; df_nopeaks(df_tmp < -std(df_nopeaks)) = NaN;
dfSTD_tmp = mean(df_nopeaks, 'omitnan') + std(df_nopeaks, 'omitnan');

% Find peaks using initial parameters
w0 = ss; p0 = dfSTD_tmp; % Original value
PK = findpeaks(f_mg, 'MinPeakWidth', w0, 'MinPeakProminence', p0, 'MinPeakDistance', ss);
PK = length(PK);

% Test if adequate number of peaks have been identified
if ~(PK >= 10 && PK <= 30)
    % Initialize variables
    iter = 0;
    delpkFLAG = true; % Flag indicating successful change in number of peaks
    dw = 1/2; dp = 1/2; % Change
    w = w0; p = p0; % Previous value

    % Iterate
    while true
        % Update changes
        if delpkFLAG
            % Revert to original values
            if PK < 10 % Increase
                dw = 1/2; dp = 1/2;
            else % Decrease
                dw = 2; dp = 2;
            end
        else
            % Change
            if PK < 10 % Increase
                dw = dw/2; dp = dp/2;
            else % Decrease
                dw = dw*2; dp = dp*2;
            end
        end

        % Change width
        w1 = w*dw; % New value
        tmp_PKw = findpeaks(f_mg, 'MinPeakWidth', w1, 'MinPeakProminence', p, 'MinPeakDistance', ss);
        tmp_PKw = length(tmp_PKw);

        % Change prominence
        p1 = p*dp; % New value
        tmp_PKp = findpeaks(f_mg, 'MinPeakWidth', w, 'MinPeakProminence', p1, 'MinPeakDistance', ss);
        tmp_PKp = length(tmp_PKp);

        % Check for a change with either
        if PK ~= tmp_PKw || PK ~= tmp_PKp
            % Set flag
            delpkFLAG = true;

            % Determine what to look for
            if PK < 10 % need to increase
                % Determine which changed the most
                if tmp_PKw > tmp_PKp % width caused greater change
                    % Set peak number
                    PK = tmp_PKw;

                    % Update width
                    w = w1;
                else % prominence caused greater change
                    % Set peak number
                    PK = tmp_PKp;

                    % Update prominence
                    p = p1;
                end
            else % need to decrease
                % Determine which changed the most
                if tmp_PKw < tmp_PKp % width caused greater change
                    % Set peak number
                    PK = tmp_PKw;

                    % Update width
                    w = w1;
                else % prominence caused greater change
                    % Set peak number
                    PK = tmp_PKp;

                    % Update prominence
                    p = p1;
                end
            end

            % Check number of peaks
            if PK >= 10 && PK <= 30 % no change needed
                break
            else
                iter = iter + 1;
            end
        else % No change
            % Set flag
            delpkFLAG = false;

            % Shift iterator
            iter = iter + 1;
        end

        % Check for max iterations
        if iter > 200
            fprintf('\nMaximum iterations reached. While loop will be stopped at %i peaks.', PK);
            break
        end
    end
    [pk, lc, wd, pr] = findpeaks(f_mg, 'MinPeakWidth', w, 'MinPeakProminence', p, 'MinPeakDistance', ss);
    clear w1 p1 w p w0 p0 dw dp iter delpkFLAG tmp_PKw tmp_PKp
else
    [pk, lc, wd, pr] = findpeaks(f_mg, 'MinPeakWidth', w0, 'MinPeakProminence', p0, 'MinPeakDistance', ss);
    clear w0 p0
end
clear PK df_tmp dfSTD_tmp

% Combine into structure
MagPeaks = struct('Height', pk(1), 'Location', lc(1), 'FWHM', wd(1), 'Prominence', pr(1));
for i = 2:length(pk)
    MagPeaks(i).Height = pk(i); MagPeaks(i).Location = lc(i); MagPeaks(i).FWHM = wd(i); MagPeaks(i).Prominence = pr(i);
end
clear pk lc wd pr

% Create upper bound phase check
ub = dfSTD_ph*4;
if ub > 2*pi
    ub = 3*pi;
    if ub > 2*pi
        ub = 2.5*pi;
        if ub > 2*pi
            ub = 2*pi;
        end
    end
end

% Remove obvious noise points
pFLAG = true(length(MagPeaks),1);
for i = 1:length(MagPeaks)
    % Check if too close to left hand side
    % - can't effectively evaluate, & unlikely peak would be here
    if MagPeaks(i).Location < ss
        pFLAG(i) = false;

    % Check if too close to right hand side
    % - can't effectively evaluate, & unlikely peak would be here
    elseif MagPeaks(i).Location > length(f)-ss
        pFLAG(i) = false;
        
    else % In acceptable area
        % Set region
        s = BuffIndex(MagPeaks(i).Location, ss, 2);

        % Find phase STD
        MagPeaks(i).PhaseSTD = std(df_ph(s));

        % Check STD of phase df & AD of real df
        % - STD greater than 1/4*global but less than 2*global is likely noise
        % - some peaks will be off a bit in magnitude data, so check if AD of at least half the points are less then dfSTD
        if MagPeaks(i).PhaseSTD > dfSTD_ph/3 && MagPeaks(i).PhaseSTD < ub && sum(dfAD(s)>dfSTD) < ss/2
            pFLAG(i) = false;
        end
    end
end
MagPeaks = MagPeaks(pFLAG);

% Extract locations and place into array
lc = zeros(length(MagPeaks),1);
for i = 1:length(MagPeaks)
    lc(i) = MagPeaks(i).Location;
end
lc = lc.';

clear ub

%% Now find locations in real data (might be off a bit due to smoothing & peak broadening)
% Plot data
if pltFLAG
    figure; plot(1:length(f), f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on; 
    yline(fSTD, 'Color', [0 .447 .741], 'LineStyle', '--');
    plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2);
    plot(1:length(f), dfAD.*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184]);
    yline(dfSTD*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184], 'LineStyle', '--');
    plot(1:length(f), dfsigns*(max(df)/2), 'Color', [.929 .694 .125]);
    L = legend('f', 'fSTD', 'df', 'dfAD', 'dfSTD', 'df Sign'); L.AutoUpdate = 'off'; title('All Peaks')
end

% Loop through located peaks to determine start & stop points
idx = 0;
Peaks0 = [];
for i = lc
    % Verify current peak is not within previous peak
    if idx > 0
        if i < Peaks0(idx).Stop(1)
            continue
        end
    end

    % Set indexer
    idx = idx + 1;
    
    % Find peak location in real data (will be slightly off from magnitude data due to smoothing/filtering
    % - First, grab location via identifying local max
    % - Avoids noise, but real peak may be outside sample
    ilc = find(f == max(f(BuffIndex(i,ss,4))));

    % Move point-by-point if needed
    sdir = 0;
    if f(ilc-1) > f(ilc)
        ilc = ilc-1;
        sdir = -1;
    elseif f(ilc+1) > f(ilc)
        ilc = ilc+1;
        sdir = 1;
    end
    if sdir ~= 0
        % Loop through points until peak is found
        while true
            % Check next point
            if f(ilc+sdir) > f(ilc)
                ilc = ilc + sdir;
            else
                break
            end
        end
    end

    % Add data into structure
    Peaks0(idx).Start = []; %#ok<*AGROW>
    Peaks0(idx).Peak(1) = ilc;
    Peaks0(idx).Stop = [];
    Peaks0(idx).Width = [];
    Peaks0(idx).StrikesSTD = [];
    Peaks0(idx).StrikesSign = [];
    Peaks0(idx).Reopened = [];

    % First, move left to find beginning
    Peaks0 = SearchLimit(Peaks0, idx, ss, f, fSTD, dfAD, dfSTD, dfsigns, false, pltFLAG);

    % Transfer strikes & reset in structure
    tmp_sSTD = Peaks0(idx).StrikesSTD;
    tmp_sSign = Peaks0(idx).StrikesSign;
    tmp_sOpen = Peaks0(idx).Reopened;
    Peaks0(idx).StrikesSTD = [];
    Peaks0(idx).StrikesSign = [];
    Peaks0(idx).Reopened = [];

    % Then, move right to find end
    Peaks0 = SearchLimit(Peaks0, idx, ss, f, fSTD, dfAD, dfSTD, dfsigns, true, pltFLAG);

    % Recombine strikes
    Peaks0(idx).StrikesSTD = [tmp_sSTD, Peaks0(idx).StrikesSTD];
    Peaks0(idx).StrikesSign = [tmp_sSign, Peaks0(idx).StrikesSign];
    Peaks0(idx).Reopened = [tmp_sOpen, Peaks0(idx).Reopened];

    % Ensure strikes are only for points contained within peak range
    Peaks0(idx).StrikesSTD = sort(Peaks0(idx).StrikesSTD(Peaks0(idx).StrikesSTD >= Peaks0(idx).Start(1)));
    Peaks0(idx).StrikesSign = sort(Peaks0(idx).StrikesSign(Peaks0(idx).StrikesSign >= Peaks0(idx).Start(1)));

    % Find true peak
    Peaks0(idx).Peak(1) = find(f == max(f(Peaks0(idx).Start(1):Peaks0(idx).Stop(1))));

    % Find value at indices
    Peaks0(idx).Start(2) = f(Peaks0(idx).Start(1));
    Peaks0(idx).Peak(2) = f(Peaks0(idx).Peak(1));
    Peaks0(idx).Stop(2) = f(Peaks0(idx).Stop(1));

    % Calculate peak info
    Peaks0 = CalculatePeakInfo(Peaks0, idx, f, df, dfsigns);

    % Check if bad peak
    if Peaks0(idx).TotalStrikes > Peaks0(idx).Width*0.75
        if pltFLAG
            % Remove from plot
            plotdata = findobj(gca, 'Type', 'Line'); 
            for j = 4:-1:1; delete(plotdata(j)); end
        end

        % Remove from structure
        if idx == 1
            Peaks0 = [];
        else
            Peaks0 = Peaks0(1:idx-1);
        end

        % Shift back indexer
        idx = idx - 1;

    else
        if pltFLAG
            % Plot peak
            PlotPoints(Peaks0(idx).Start(2), Peaks0(idx).Start(1), idx) 
            PlotPoints(Peaks0(idx).Stop(2), Peaks0(idx).Stop(1), idx)
        end
    end
end

% Verification to ensure peaks start & stop at proper locations
for i = 2:length(Peaks0)
    % Check if start & stop are flipped
    if Peaks0(i).Start(1) < Peaks0(i-1).Stop(1)
        % Flip points
        truestart = Peaks0(i-1).Stop;
        trueend = Peaks0(i).Start;
        Peaks0(i-1).Stop = trueend;
        Peaks0(i).Start = truestart;

        % Recalculate width
        Peaks0 = CalculatePeakInfo(Peaks0, i-1:i, f, df, dfsigns);
    end
end

% %% Find peaks (old)
% % Initialize stepping variable
% PrevPend = 1;
% 
% % Initialize flags
% OpenPeakFLAG = false;
% previouspeakcombination = false;
% maxSTD = 5;
% maxSign = 3;
% 
% % Plot data
% figure; plot(1:length(f), f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on; 
% yline(fSTD, 'Color', [0 .447 .741], 'LineStyle', '--'); 
% plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2); 
% plot(1:length(f), dfAD.*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184]); 
% yline(dfSTD*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184], 'LineStyle', '--'); 
% plot(1:length(f), dfsigns*(max(df)/2), 'Color', [.929 .694 .125]); 
% L = legend('f', 'fSTD', 'df', 'dfAD', 'dfSTD', 'df Sign'); L.AutoUpdate = 'off'; title('All Peaks')
% 
% % Initialize peak variables
% idx = 0;
% idx_old = 0;
% Peak = [];
% 
% % Debug counters
% IFS = nan(12,length(df));
% % Determine start & end points
% i = ss+1;
% while i <= length(df)-ss
%     % Set current sample
%     s = i:i+ss-1;
% 
%     % Check for no already open peak
%     if ~OpenPeakFLAG
%         % Check if current point is a previous peak
%         if ~finval(i)
%             % Shift peak indexers
%             idx = idx + 1;
%             idx_old = idx_old + 1;
% 
%             % Add in old peak
%             Peak(idx) = Peak_old(idx_old);
% 
%             % Shift index
%             i = Peak(idx).Stop + 1;
% 
%             % Reset sample
%             s = i:i+ss-1;
%         end
% 
%         % If 75% of next points are same sign                                              - or maybe sum(df(s)) >= 5*dfSTD
%         % - consistent slope could indicate start of peak
%         if sum(dfsigns(s)) >= ss*0.75 || sum(dfsigns(s)) <= ss*0.25
%             IFS(1,i) = 1;
%             % Check if 50% of slope & current point are positive                           - and next ss/10 points > 0
%             % - should be the case for start of a new peak
%             if sum(dfsigns(s)) >= ss*0.5 && dfsigns(i)
%                 IFS(2,i) = 1;
%                 % Check if i is within ss/2 of previous peak or current std is less than 2*previous std
%                 % - std should be high, relative to noise, near peaks
%                 % - if too close to other peaks, std check doesn't work
%                 if i-PrevPend <= ss/2 || (std(df(i:i+round(ss/2))) > 2*std(df(s-ss)) || mean(dfAD(s)) > 2*mean(dfAD(s-ss)))
%                     IFS(3,i) = 2;
%                     % Shift peak indexer
%                     idx = idx + 1;
% 
%                     % Update flags & strikes
%                     OpenPeakFLAG = true;
%                     Peak(idx).Start = []; %#ok<*SAGROW>
%                     Peak(idx).Peak = [];
%                     Peak(idx).Stop = [];
%                     Peak(idx).Width = [];
%                     Peak(idx).StrikesSTD = [];
%                     Peak(idx).StrikesSign = [];
%                     Peak(idx).Reopened = []; 
%                     maxSTD = 5;
%                     maxSign = 3;
% 
%                     % Set current sign
%                     pksign = true;
% 
%                     % Add peak to list
%                     Peak(idx).Start = i;
% 
%                     % Shift loop indexer
%                     i = i + ss;
% 
%                     % Plot
%                     if rem(idx,2) == 0
%                         plot(Peak(idx).Start, f(Peak(idx).Start), 'r.', 'MarkerSize', 20)
%                     else
%                         plot(Peak(idx).Start, f(Peak(idx).Start), 'b.', 'MarkerSize', 20)
%                     end
%                 else
%                     IFS(3,i) = 0;
%                     % Shift loop indexer
%                     i = i + 1;
%                 end
%             % Check if peak has already been found & average AD of surrounding pts is greater than STD
%             % - previous peak required for following tests
%             % - low deviation likely indicates just noise floor with coincidentally consistent slope
%             elseif idx > 0 && mean(dfAD(BuffIndex(i,ss,2))) > dfSTD
%                 IFS(2,i) = 0;
%                 IFS(4,i) = 1;
%                 % Case where slope is negative, check if previous peak was less than ss/2 away & current peak is less than 3*LW
%                 % - indicates peak was likely not completely captured or peaks are overlapping
%                 if i-PrevPend <= ss/2 && length(Peak(idx).Reopened) < 3 && i-Peak(idx).Start+1 < ss*3
%                     IFS(5,i) = 2;
%                     % Update flags & strikes
%                     OpenPeakFLAG = true;
%                     Peak(idx).Reopened = [Peak(idx).Reopened, i];
%                     if length(Peak(idx).StrikesSTD) >= maxSTD
%                         maxSTD = length(Peak(idx).StrikesSTD) + 1;
%                     end
%                     if length(Peak(idx).StrikesSign) >= maxSign
%                         maxSign = length(Peak(idx).StrikesSign) + 1;
%                     end
% 
%                     % Reset current sign
%                     pksign = false;
% 
%                     % Only do big jump if larger peak is before current point
%                     % - less likely to result in overjumping
%                     if max(f(i-round(ss/2):s(end))) > max(f(s(end)+1:s(end)+ss))
%                         % Shift loop indexer
%                         i = i + ss;
%                     else
%                         % Shift loop indexer
%                         i = i + round(ss/4);
%                     end
% 
%                     % Hide previous plot
%                     plotdata = get(gca, 'Children');
%                     set(plotdata(1), 'visible', 'off');
%                 else % Nearest peak is far away or reopened too many times or too large of LW
%                     IFS(5,i) = 0;
%                     % Shift loop indexer
%                     i = i + 1;
%                 end
%             else % No consistent increase & no previous peak, no peak
%                 IFS(2,i) = 0;
%                 IFS(4,i) = 0;
%                 % Shift loop indexer
%                 i = i + 1;
%             end
%         else % No long stretch of similar sign, no peak
%             IFS(1,i) = 0;
%             % Shift loop indexer
%             i = i + 1;
%         end
%     else % open peak
%         % Ensure current point is not on a previous peak
%         if finval(i)
%             % Check for actual peak, including a little past current index 
%             if all(finval(Peak(idx).Start:i+floor(ss/2)))
%                 Peak(idx).Peak = find(f == max(f(Peak(idx).Start:i+floor(ss/2))));
%             end
%             % Check if AD of any points within ss/4 are 2x greater than STD
%             % - means that there are very large changes, not noise
%             if any(dfAD(BuffIndex(i,ss,8)) >= 2*dfSTD)
%                 IFS(7,i) = 1;
%                 IFS(8,i) = 0;
% 
%                 % Check if next 75% are all positive, current peak width is 2*ss, and past peak peak
%                 % - means we have passed the peak & starting a new one
%                 if sum(dfsigns(s)) >= ss*0.75 && i-Peak(idx).Start+1 > ss*2 && i > Peak(idx).Peak
%                     IFS(8,i) = 1;
%                     % Add sign strikes to close peak
%                     itmp = i-1; % Temporary index
%                     while maxSign > length(Peak(idx).StrikesSign)
%                         itmp = itmp + 1;
%                         % Check if positive
%                         if dfsigns(itmp)
%                             % Add strike
%                             Peak(idx).StrikesSign = [Peak(idx).StrikesSign, itmp];
%                         end
%                     end
%                     clear itmp
%                 else
%                     % Update sign (to catch downward slope from other side of peak if not caught by other steps)
%                     pksign = dfsigns(i);
%                 end
% 
%                 % Shift loop indexer
%                 i = i + 1;
% 
%             % Check if AD of any points within ss/2 are greater than STD
%             % - means that there are large changes, possibly not noise 
%             elseif any(dfAD(BuffIndex(i,ss,4)) >= dfSTD)
%                 IFS(7,i) = 0;
%                 IFS(9,i) = 1;
%                 % Check for change in sign
%                 % - could indicate index is in noise region
%                 if dfsigns(i) ~= pksign
%                     IFS(10,i) = 1;
%                     % Check if next 75% are all same sign & before peak point
%                     % - would indicate getting close to peak point
%                     if (sum(dfsigns(s)) >= ss*0.75 || sum(dfsigns(s)) <= ss*0.25) && any(i-round(ss/4):i <= Peak(idx).Peak)
%                         IFS(11,i) = 1;
%                         % Set new sign
%                         pksign = ~pksign;
% 
%                         % Shift indexer
%                         i = i + round(ss/2);
%                     % Check if all signs are negative & pksign is positive
%                     % - could mean we missed the downward slope
%                     elseif sum(dfsigns(s)) == 0 && pksign
%                         IFS(11,i) = 0;
%                         IFS(12,i) = 1;
%                         % Set new sign
%                         pksign = false;
% 
%                         % Shift indexer
%                         i = i + ss;
%                     else
%                         IFS(11,i) = 0;
%                         IFS(12,i) = 0;
%                         % Add strike
%                         Peak(idx).StrikesSign = [Peak(idx).StrikesSign, i];
% 
%                         % Shift loop indexer
%                         i = i + 1;
%                     end
%                 else % AD of any nearby points is greater than STD but no change in sign, not noise
%                     IFS(10,i) = 0;
%                     % Shift loop indexer
%                     i = i + 1;
%                 end
%             else % AD of any nearby points are all less than STD, definitely noise
%                 IFS(7,i) = 0;
%                 IFS(9,i) = 0;
%                 % Add strike
%                 Peak(idx).StrikesSTD = [Peak(idx).StrikesSTD, i];
% 
%                 % Shift loop indexer
%                 i = i + 1;
%             end
%         else
%             % Add strikes to close peak
%             itmp = i-1; % Temporary index
%             while maxSign > length(Peak(idx).StrikesSign) || maxSTD > length(Peak(idx).StrikesSTD)
%                 itmp = itmp + 1;
% 
%                 % Check if positive
%                 if dfsigns(itmp)
%                     % Add strike
%                     Peak(idx).StrikesSign = [Peak(idx).StrikesSign, itmp];
%                 end
% 
%                 % Check STD
%                 if dfAD(itmp) < dfSTD
%                     % Add strike
%                     Peak(idx).StrikesSTD = [Peak(idx).StrikesSTD, itmp];
%                 end
%             end
%             clear itmp
%         end
% 
%         % Check for too many strikes
%         if length(Peak(idx).StrikesSign) >= maxSign || length(Peak(idx).StrikesSTD) >= maxSTD
%             % Close peak
%             [Peak, idx, PrevPend, OpenPeakFLAG] = ClosePeak(f, Peak, i, idx);
% 
%             % Shift back indexer (to account for neighboring peaks)
%             i = i - 1;
%         end
%     end
% end
% IFS = IFS + 1; IFS(isnan(IFS)) = 0; figure; imagesc(IFS); colormap hot;
% 
% % Verify last peak is closed
% if isempty(Peak(end).Stop)
%     Peak(end).Stop = length(df)-1;
% end
% 
% 

%% Combine doublets/triplets
if pltFLAG
    % Plot
    figure; plot(1:length(f), f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on; 
    plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2);
    plot(1:length(f), dfsigns*(max(df)/2), 'Color', [.929 .694 .125]);
    L = legend('f', 'df', 'df Sign'); L.AutoUpdate = 'off';
    title('Cleaned Peaks')
end

% Measure difference in height between left/peak, right/peak, & left/right
Ptdif = ones(length(Peaks0), 4);
for i = 1:length(Peaks0)
    Ptdif(i,1) = Peaks0(i).Peak(2) - Peaks0(i).Start(2);
    Ptdif(i,2) = Peaks0(i).Peak(2) - Peaks0(i).Stop(2);
    Ptdif(i,3) = abs(Peaks0(i).Stop(2) - Peaks0(i).Start(2));

    % Normalize
    Ptdif(i,:) = Ptdif(i,:)./max(Ptdif(i,:)); 

    % Find min
    Ptdif(i,4) = find(Ptdif(i,:) == min(Ptdif(i,:)));
end

% Check for doublets/triplets
if exist('Peaks1','var'); clear Peaks1; end
combFLAG = false(length(Peaks0),1);
idx = 0;
for i = 1:length(Peaks0)
    % Reset variables
    s = [];
    np = 0;

    % Check if peak was combined
    if ~combFLAG(i)
        % Check if LR difference is not the smallest for current peak
        if Ptdif(i,4) ~= 3
            % Check which side the min is
            if Ptdif(i,4) == 1
                dirFLAG = 'left';
            else
                dirFLAG = 'right';
            end

            % Loop based on min side
            if strcmp(dirFLAG, 'left') && i > 1
                for j = i-1:-1:1
                    % Check for doublet
                    % - peak immediately left should have a min on rightside
                    if j == i-1 && Ptdif(j,4) == 2
                        % Set doublet flag
                        np = -1;

                    % Check for triplet
                    % - peak two over left should have a min on rightside
                    elseif j == i-2 && Ptdif(j,4) == 2
                        % Set triplet flag
                        np = -2;
                    end
                end
            elseif strcmp(dirFLAG, 'right') && i < length(Peaks0)
                for j = i+1:length(Peaks0)
                    % Check for doublet
                    % - peak immediately right should have a min on leftside
                    if j == i+1 && Ptdif(j,4) == 1
                        % Set doublet flag
                        np = 1;

                    % Check for triplet
                    % - peak two over right should have a min on leftside
                    elseif j == i+2 && Ptdif(j,4) == 1
                        % Set triplet flag
                        np = 2;
                    end
                end
            else % too close to start or end
                % Shift iterator
                idx = idx + 1;

                % Append peak to structure
                Peaks1(idx) = Peaks0(i);

                if pltFLAG
                    % Plot
                    PlotPoints(Peaks1(idx).Start(2), Peaks1(idx).Start(1), idx) 
                    PlotPoints(Peaks1(idx).Stop(2), Peaks1(idx).Stop(1), idx)
                end

                % Skip remainder of loop to avoid index errors
                continue
            end

            % Set comparison points based on whether checking left or right
            if strcmp(dirFLAG, 'left')
                L = i + np;
                R = i;
                if np == 0
                    L = i - 1;
                end
            else
                L = i;
                R = i + np;
                if np == 0
                    R = i + 1;
                end
            end
            Lshorthalf = Peaks0(L).UpperHalf.Length;
            Rshorthalf = Peaks0(R).LowerHalf.Length;
            Llonghalf = Peaks0(L).LowerHalf.Length;
            Rlonghalf = Peaks0(R).UpperHalf.Length;
            Lgap = Peaks0(L).Stop(1);
            Rgap = Peaks0(R).Start(1);

            % Check for no doublets/triplets
            if np == 0
                % Check if poor shim
                % - the case if height difference between shared start & stop is small
                if abs(Peaks0(L).Start(2)-Peaks0(R).Start(2)) < (Peaks0(L).Peak(2)-Peaks0(L).Start(2))/2 && abs(Peaks0(L).Start(2)-Peaks0(R).Start(2)) < (Peaks0(R).Peak(2)-Peaks0(R).Stop(2))/2
                    % Set flag & sample
                    if strcmp(dirFLAG, 'left')
                        s = i-1:i;
                    else
                        s = i:i+1;
                    end
                    combFLAG(s) = true;
                end

            % Check if the shared halves of first & last peaks are 1/2 as long as unshared halves & distance between peak start & ends is less than ss
            elseif Rshorthalf < Rlonghalf/2 && Lshorthalf < Llonghalf/2 && Rgap-Lgap < ss
                % Set flag & sample
                if strcmp(dirFLAG, 'left')
                    s = i+np:i;
                else
                    s = i:i+np;
                end
                combFLAG(s) = true;
            end
        end

        % Shift iterator
        idx = idx + 1;

        % Check if need to combine peaks
        if ~isempty(s)
            % Check if need to remove already indexed peak
            if strcmp(dirFLAG, 'left')
                if pltFLAG
                    % Remove from plot
                    plotdata = findobj(gca, 'Type', 'Line'); 
                    for j = 2:-1:1; delete(plotdata(j)); end
                end

                % Determine number of points to remove
                numpts = length(s)-1;

                % Remove peaks
                Peaks1 = Peaks1(1:end-numpts);

                % Shift back iterator to account for removed peaks
                idx = idx - numpts;
            end

            % Combine peaks
            Peaks1(idx).Start = Peaks0(s(1)).Start;
            itmp = 0;
            for j = s
                itmp = itmp + 1;
                Peaks1(idx).Peak(:,itmp) = Peaks0(j).Peak;
            end
            clear itmp
            Peaks1(idx).Stop = Peaks0(s(end)).Stop;
            Peaks1(idx).Width = Peaks1(idx).Stop(1) - Peaks1(idx).Start(1) + 1;
            Peaks1(idx).StrikesSTD = []; Peaks1(idx).StrikesSign = []; Peaks1(idx).Reopened = []; Peaks1(idx).TotalStrikes = 0;
            for j = s
                Peaks1(idx).StrikesSTD = [Peaks1(idx).StrikesSTD, Peaks0(j).StrikesSTD];
                Peaks1(idx).StrikesSign = [Peaks1(idx).StrikesSign, Peaks0(j).StrikesSign];
                Peaks1(idx).Reopened = [Peaks1(idx).Reopened, Peaks0(j).Reopened];
                Peaks1(idx).TotalStrikes = Peaks1(idx).TotalStrikes + Peaks0(j).TotalStrikes;
            end
            Peaks1(idx).LowerHalf = Peaks0(i).LowerHalf;
            Peaks1(idx).UpperHalf = Peaks0(i+np).UpperHalf;

            % Ensure NaNs do not endure through combination
            if isnan(Peaks1(idx).LowerHalf.Length)
                Peaks1(idx).LowerHalf = struct('Length', Peaks1(idx).Peak(1,1) - Peaks1(idx).Start(1), ...
                                                  'Sum', sum(df(Peaks1(idx).Start(1):Peaks1(idx).Peak(1,1))), ...
                                                 'Sign', sum(dfsigns(Peaks1(idx).Start(1):Peaks1(idx).Peak(1,1))));
            end
            if isnan(Peaks1(idx).UpperHalf.Length)
                Peaks1(idx).UpperHalf = struct('Length', Peaks1(idx).Stop(1)+1 - Peaks1(idx).Peak(1,end), ...
                                                  'Sum', sum(df(Peaks1(idx).Peak(1,end):Peaks1(idx).Stop(1)+1)), ...
                                                 'Sign', sum(dfsigns(Peaks1(idx).Peak(1,end):Peaks1(idx).Stop(1)+1)));
            end

        else % not doublet/triplet peak
            % Append peak to structure
            Peaks1(idx) = Peaks0(i);
        end

        if pltFLAG
            % Plot
            PlotPoints(Peaks1(idx).Start(2), Peaks1(idx).Start(1), idx) 
            PlotPoints(Peaks1(idx).Stop(2), Peaks1(idx).Stop(1), idx)
        end
    end
end


%% Clean up peaks
if pltFLAG
    % Plot
    figure; plot(1:length(f), f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on; 
    plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2);
    plot(1:length(f), dfsigns*(max(df)/2), 'Color', [.929 .694 .125]);
    L = legend('f', 'df', 'df Sign'); L.AutoUpdate = 'off';
    title('Final Peaks')
end

% Clean up peaks & get more information for analysis
if exist('Peaks2','var'); clear Peaks2; end
idx = 0;
for i = 1:length(Peaks1)
    % Reset first if FLAG
    newpeak = false;

    % Check if one side is not >> than other side
    if Peaks1(i).UpperHalf.Length/Peaks1(i).LowerHalf.Length < 5 && Peaks1(i).LowerHalf.Length/Peaks1(i).UpperHalf.Length < 5
        % Reset flags
        leftFLAG = false;
        rightFLAG = false;

        % Check if peak is a true peak (lower half should be positive slope)
        if Peaks1(i).LowerHalf.Sign >= Peaks1(i).LowerHalf.Length*0.75
            leftFLAG = true;
        % Some allowence is made for doublet/triplet peaks & poor shim
        elseif Peaks1(i).LowerHalf.Sign >= Peaks1(i).LowerHalf.Length*0.675
            % Create normalized peak value
            normPk = Peaks1(i).Peak(2) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Find max value to left of peak
            dbPk = max(f(round(Peaks1(i).Start(1)+ss/4):round(Peaks1(i).Peak(1)-ss/4))) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Check for possible doublet (if normalized point is 90% of main peak)
            if dbPk/normPk >= 0.9
                leftFLAG = true;
            end
        end

        % Check if peak is a true peak (upper half should be negative slope)
        if Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25
            rightFLAG = true;
        % Some allowence is made for doublet/triplet peaks & poor shim
        elseif Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.325
            % Create normalized peak value
            normPk = Peaks1(i).Peak(2) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Find max value to right of peak
            dbPk = max(f(round(Peaks1(i).Peak(1)+ss/4):round(Peaks1(i).Stop(1)-ss/4))) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Check for possible doublet (if normalized point is 90% of main peak)
            if dbPk/normPk >= 0.9
                rightFLAG = true;
            end
        end

        if leftFLAG && rightFLAG
            % Shift indexer & append peak data to final array
            idx = idx + 1;
            Peaks2(idx) = Peaks1(i); %#ok<*SAGROW>

            % Set flag
            newpeak = true;

            % Plot peak
            if debugFLAG
                PlotPoints(Peaks2(idx).Start(2), Peaks2(idx).Start(1), idx)
                PlotPoints(Peaks2(idx).Stop(2), Peaks2(idx).Stop(1), idx)
            end
        end
    end

    % Check if upper half is negative slope & slope between two peak bases is mostly positive
    if idx > 0 && ~newpeak && Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25 && sum(dfsigns(Peaks2(idx).Stop(1)+1:Peaks1(i).Start(1))) >= (Peaks1(i).Start(1) - (Peaks2(idx).Stop(1)+1))*0.75 && ~previouspeakcombination
        % Combine two peaks (because peaks are likely overlapping due to poor shim and/or poor spectral resolution)
        Peaks2(idx).Stop = Peaks1(i).Stop;
        Peaks2(idx).StrikesSTD = unique(sort([Peaks1(i).StrikesSTD, Peaks2(idx).StrikesSTD]));
        Peaks2(idx).StrikesSign = unique(sort([Peaks1(i).StrikesSign, Peaks2(idx).StrikesSign]));
        Peaks2 = CalculatePeakInfo(Peaks2, idx, f, df, dfsigns);

        % Set flag
        previouspeakcombination = true;

        % Plot points
        if debugFLAG
            PlotPoints(Peaks2(idx).Stop(2), Peaks2(idx).Stop(1), idx, 1)
        end
    else
        % Reset flag
        previouspeakcombination = false;
    end

    % % See if peak is a true peak (lower half should be positive slope, upper half should be negative slope, & one side is not >> than other side)
    % % Some allowence is made for doublet/triplet peaks & poor shim
    % if Peaks1(i).LowerHalf.Sign >= Peaks1(i).LowerHalf.Length*0.75 && Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25 && abs(Peaks1(i).UpperHalf.Length/Peaks1(i).LowerHalf.Length) < 5 && abs(Peaks1(i).LowerHalf.Length/Peaks1(i).UpperHalf.Length) < 5
    %     % Shift indexer & append peak data to final array
    %     idx = idx + 1;
    %     Peaks2(idx) = Peaks1(i); %#ok<*SAGROW>
    % 
    %     % Reset flag
    %     previouspeakcombination = false;
    % 
    %     if pltFLAG
    %         % Plot peak
    %         PlotPoints(Peaks2(idx).Start(2), Peaks2(idx).Start(1), idx) 
    %         PlotPoints(Peaks2(idx).Stop(2), Peaks2(idx).Stop(1), idx)
    %     end
    % 
    % % Check if upper half is negative slope & slope between two peak bases is mostly positive
    % elseif idx > 0 && Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25 && sum(dfsigns(Peaks2(idx).Stop(1)+1:Peaks1(i).Start(1))) >= (Peaks1(i).Start(1) - (Peaks2(idx).Stop(1)+1))*0.75 && ~previouspeakcombination
    %     % Combine two peaks (because peaks are likely overlapping due to poor shim and/or poor spectral resolution)
    %     Peaks2(idx).Stop = Peaks1(i).Stop;
    %     Peaks2(idx).StrikesSTD = unique(sort([Peaks1(i).StrikesSTD, Peaks2(idx).StrikesSTD]));
    %     Peaks2(idx).StrikesSign = unique(sort([Peaks1(i).StrikesSign, Peaks2(idx).StrikesSign]));
    %     Peaks2 = CalculatePeakInfo(Peaks2, idx, f, df, dfsigns);
    % 
    %     % Change flag
    %     previouspeakcombination = true;
    % 
    %     if pltFLAG
    %         % Plot points
    %         PlotPoints(Peaks2(idx).Stop(2), Peaks2(idx).Stop(1), idx, 1)
    %     end
    % else
    %     % Reset flag
    %     previouspeakcombination = false;
    % end
end

%% Shift starts & stops to *very* local mins
for idx = 1:length(Peaks2)
    f_tmp = NaN(size(f_comp));
    f_tmp(BuffIndex(Peaks2(idx).Start(1),ss,8)) = dfAD(BuffIndex(Peaks2(idx).Start(1),ss,8));
    [~, Peaks2(idx).Start(1)] = min(f_tmp, [], 'all', 'omitnan');
    f_tmp = NaN(size(f_comp));
    f_tmp(BuffIndex(Peaks2(idx).Stop(1),ss,8)) = dfAD(BuffIndex(Peaks2(idx).Stop(1),ss,8));
    [~, Peaks2(idx).Stop(1)] = min(f_tmp, [], 'all', 'omitnan');
end

% Do final verification to ensure peaks start & stop at proper locations
for i = 2:length(Peaks2)
    % Check if start & stop are flipped
    if Peaks2(i).Start(1) < Peaks2(i-1).Stop(1)
        % Flip points
        truestart = Peaks2(i-1).Stop(1);
        trueend = Peaks2(i).Start(1);
        Peaks2(i-1).Stop(1) = trueend;
        Peaks2(i).Start(1) = truestart;

        % Recalculate width
        Peaks2 = CalculatePeakInfo(Peaks2, i, f, df, dfsigns);
        Peaks2 = CalculatePeakInfo(Peaks2, i-1, f, df, dfsigns);
    end
end

% Check for overlapping peaks
keeps = true(length(Peaks2),1);
for i = 2:length(Peaks2)
    if Peaks2(i-1).Stop(1) > Peaks2(i).Start(1)
        % Combine peaks
        Peaks2(i).Start(1) = Peaks2(i-1).Start(1);

        % Flag old peak
        keeps(i-1) = false;
    end
end

% Remove old peaks
Peaks2 = Peaks2(keeps);

%% Perform baseline correction
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
Baseline1 = RemovePeak(Peaks2, Baseline0);

% Create structure for model
if exist('TruePeaks','var'); clear TruePeaks; end
for i = 1:length(Peaks2)
    TruePeaks(i).Start = Peaks2(i).Start(1);
    TruePeaks(i).End = Peaks2(i).Stop(1);
end

% Model baseline
Baseline2 = New_InterativeModel(Baseline1, 1000, 0, Baseline0, TruePeaks);

% Obtain corrected spectra
Final_specs = complex(real(f_comp).' - Baseline2, imag(f_comp).').';

if pltFLAG
    plot(Baseline2.'); 
    % figure; plot(real(f_comp)); axis tight; hold on; plot(Baseline1.');  plot(Baseline2.'); legend('Raw Data', 'Interpolated Baseline',  'Smoothed Baseline'); %#ok<UNRCH>
    % figure; plot(real(Final_specs)); axis tight; title('Corrected Spectra')
    % figure; plot(real(f_comp-Final_specs)); axis tight; title('Difference'); ylim([min(real(Final_specs)), max(real(Final_specs))])
end
%% Functions

    function coefs = cwt_custom(ySIG, scales)
        % Number of iterations
        iters = 10;
        
        % haar reconstruction LPF & HPF
        Lo_R = [1/sqrt(2), 1/sqrt(2)];
        Hi_R = [1/sqrt(2), -1/sqrt(2)];
        
        
        % Compute N-step reconstructed haar wavelet coefficients
        psi = 1;
        lf = length(Lo_R);
        lpsi = 2+lf-2;
        psi = wconv1(dyadup(psi,0),Hi_R); psi = wkeep1(psi,lpsi,'c',0); % Upsample and convolute
        for k = 2:iters
            lpsi = 2*length(psi)+lf-2;
            psi = wconv1(dyadup(psi,0),Lo_R); psi = wkeep1(psi,lpsi,'c',0); % Upsample and convolute
        end
        
        % Approximate wavelet function
        psi = (sqrt(2)^iters)*psi;
        nb = 2^10;
        nbpts = nb+2;
        psi = [0 wkeep1(psi,nb) 0];
        
        % Generate 2^ITER points grid
        xWAV = linspace(0,(nbpts-1)/(2^iters),nbpts);
        
        % Integrate wavelet function
        stepWAV = xWAV(2)-xWAV(1);
        val_WAV = cumsum(psi)*stepWAV;
        
        % Calculate cwt of ySIG
        xWAV = xWAV-xWAV(1);
        xMaxWAV = xWAV(end);
        ySIG   = ySIG(:)';
        k = 1+floor((0:scales*xMaxWAV)/(scales*stepWAV));
        if isscalar(k), k = [1 1]; end
        coefs = sqrt(scales)*wkeep1(diff(wconv1(ySIG,fliplr(val_WAV(k)))),length(ySIG));
    end
    
    function [fSTD, df, dfAD, dfSTD, dfsigns] = PeakStats(f, scales)
        % Take derivative of filtered data
        % df = [diff(f), 0];
        df = cwt_custom(f, scales);
        
        % Determine derivative sign (true = pos, false = neg/0)
        dfsigns = sign(df) == 1;
        
        % Remove definite peaks
        df_np = df;
        df_np(df > std(df)) = NaN;
        df_np(df < -std(df)) = NaN;
        df_np(df > std(df_np)) = NaN;
        df_np(df < -std(df_np)) = NaN;
        
        % Obtain absolute deviation & standard deviation for df
        dfAD = abs(df - mean(df_np, 'omitnan'));
        dfSTD = mean(df_np, 'omitnan') + std(df_np, 'omitnan');
        
        % Obtain first standard deviation line for spectra
        fSTD = mean(f) + std(f);
    end
        
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
    
    function Peak = SearchLimit(Peak, idx, ss, f, fSTD, dfAD, dfSTD, dfsigns, rdir, pltFLAG)
        % Set shift amount
        if rdir
            sval = 1;
        else
            sval = -1;
        end
        
        % Flip dfsigns data if going left
        if ~rdir
            dfsigns = ~dfsigns;
        end
        
        % Set bounds with nearby peaks
        ppFLAG = false;
        if idx > 1 && ~rdir
            ppFLAG = true;
            iMIN = Peak(idx-1).Stop(1);
        end
        
        % Check for noisy peak (or possible doublet/triplet)
        if rdir
            pkrng = Peak(idx).Peak(1):Peak(idx).Peak(1)+round(ss/5);
        else
            pkrng = Peak(idx).Peak(1)-round(ss/5):Peak(idx).Peak(1);
        end
        pksigns = sum(dfsigns(pkrng))/length(pkrng);
        if pksigns > 0.4 && pksigns < 0.6
            % Perform initial shift to get off noisy peak
            k = Peak(idx).Peak(1) + sval*length(pkrng);
        else
            % Perform initial shift to get off peak point
            k = Peak(idx).Peak(1) + sval;
        end
        
        % Loop through points
        while true
            % Check if index is encroaching on a previous peak
            if ppFLAG
                if k <= iMIN
                    % Force shut the peak
                    Peak(idx).Reopened = [k, k, k];
                end
            end
        
            % Only check for strikes if spectra point is below first standard deviation
            % - otherwise, it might flag points around the peak
            if f(k) < fSTD && length(Peak(idx).Reopened) < 3
                % Check if point is within ss/8 of peak
                % - slope right at peak is very small, so df checks are not very accurate
                if abs(k - Peak(idx).Peak(1)) < ss/8
                    % Check if dfAD is greater than dfSTD/2
                    % - there should still be some variation at peak
                    if dfAD(k) > dfSTD/2
                        % Check if point is positive
                        % - positive slope would mean noise or doublet/triplet (but haven't coded for that yet)
                        if dfsigns(k)
                            % Add strike
                            Peak(idx).StrikesSign = [Peak(idx).StrikesSign, k];
                        end
                    else
                        % Add strike
                        Peak(idx).StrikesSTD = [Peak(idx).StrikesSTD, k];
                    end
        
                % Check if AD of any points within ss/4 are 2x greater than STD
                % - means that there are very large changes, not noise
                elseif any(dfAD(BuffIndex(k,ss,8)) >= 2*dfSTD)
                    % Set current sample
                    if rdir
                        s1 = k:k+ss-1;
                    else
                        s1 = k-ss+1:k;
                    end
        
                    % Check if next 75% are all positive and current peak width is ss
                    % - means we have passed the peak & starting a new one
                    if sum(dfsigns(s1)) >= ss*0.75 && abs(k-Peak(idx).Peak(1)+1) > ss
                        % Add strikes to close peak
                        while length(Peak(idx).StrikesSign) < 3 && length(Peak(idx).StrikesSTD) < 5
                            % Shift loop indexer
                            k = k + sval;
        
                            % Check if positive
                            if dfsigns(k)
                                % Add strike
                                Peak(idx).StrikesSign = [Peak(idx).StrikesSign, k];
                            end
        
                            % Check STD
                            if dfAD(k) < dfSTD
                                % Add strike
                                Peak(idx).StrikesSTD = [Peak(idx).StrikesSTD, k];
                            end
                        end
        
                    % Check if next 50% are all positive
                    % - means we are getting too close to noise
                    elseif sum(dfsigns(s1)) >= ss*0.5
                        % Check if current AD is less than STD or if sign is positive
                        if dfAD(k) < dfSTD || dfsigns(k)
                            % Add strike
                            Peak(idx).StrikesSign = [Peak(idx).StrikesSign, k];
                        end
                    end
            
                % Check if AD of any points within ss/2 are greater than STD
                % - means that there are large changes, possibly not noise 
                elseif any(dfAD(BuffIndex(k,ss,4)) >= dfSTD)
                    % Check for positive slope
                    % - could indicate index is in noise region
                    if dfsigns(k)
                        % Add strike
                        Peak(idx).StrikesSign = [Peak(idx).StrikesSign, k];
                    end
        
                else % AD of any nearby points are all less than STD, definitely noise
                    % Add strike
                    Peak(idx).StrikesSTD = [Peak(idx).StrikesSTD, k];
                end
            end
        
            % Check for too many strikes
            closeFLAG = false;
            if length(Peak(idx).StrikesSign) >= 3 || length(Peak(idx).StrikesSTD) >= 5 || length(Peak(idx).Reopened) >= 3
                % Check if current slope is positive
                if dfsigns(k-sval)
                    closeFLAG = true;
    
                % Check if last 10 points had STD strike
                elseif length(Peak(idx).StrikesSTD) >= 10
                    if all(Peak(idx).StrikesSTD(end-9:end) == k-sval*10:sval:k-sval)
                        closeFLAG = true;
                    end
                end
            end
    
            % Close peak
            if closeFLAG
                % Chech if already reopened too many times
                if length(Peak(idx).Reopened) < 3 && abs(k-Peak(idx).Peak(1)) < ss
                    % Set evaluation bounds for reopening
                    if rdir
                        lbd = k;
                        ubd = Peak(idx).Peak(1) + ss;
                    else
                        lbd = Peak(idx).Peak(1) - ss;
                        ubd = k;
            
                        % Shift bound back if index is encroaching on a previous peak
                        if idx > 1 && ppFLAG; lbd = ReBound(lbd, ubd, iMIN); end
                    end
            
                    % Check if there is a possible peak near location
                    % - should be local max value & have very low dfAD
                    posPeak = find(f == max(f(lbd:ubd)));
                    if posPeak ~= k && dfAD(posPeak) < dfSTD/2
                        % Set evaluation bounds for next check
                        if rdir
                            lbd = posPeak;
                            ubd = posPeak + round(ss/3);
                        else
                            lbd = posPeak - round(ss/3);
                            ubd = posPeak;
            
                            % Shift bound back if index is encroaching on a previous peak
                            if idx > 1 && ppFLAG; lbd = ReBound(lbd, ubd, iMIN); end
                        end
            
                        % Check if a value slightly after dfAD has high AD & is negative
                        posSlope = find(dfAD == max(dfAD(lbd:ubd)));
                        if dfAD(posSlope) > 2*dfSTD && ~dfsigns(posSlope)
                            % Check if closely following points are also negative
                            if rdir
                                lbd = posSlope;
                                ubd = posSlope + round(ss/4);
                            else
                                lbd = posSlope - round(ss/4);
                                ubd = posSlope;
                            end
                            if all(~dfsigns(lbd:ubd))
                                % Reopen peak
                                Peak(idx).Reopened = posSlope;
        
                                % Shift down strikes to be just under max
                                if length(Peak(idx).StrikesSign) >= 3 
                                    badstrike = Peak(idx).StrikesSign(dfAD(Peak(idx).StrikesSign) == min(dfAD(Peak(idx).StrikesSign)));
                                    Peak(idx).StrikesSign = Peak(idx).StrikesSign(Peak(idx).StrikesSign ~= badstrike);
                                elseif length(Peak(idx).StrikesSTD) >= 5
                                    badstrike = Peak(idx).StrikesSTD(dfAD(Peak(idx).StrikesSTD) == max(dfAD(Peak(idx).StrikesSTD)));
                                    Peak(idx).StrikesSTD = Peak(idx).StrikesSTD(Peak(idx).StrikesSTD ~= badstrike);
                                end
        
                                % Shift loop indexer
                                k = posSlope + sval;
        
                                % Skip remainder of loop
                                continue
                            end
                        end
                    end
                end
        
                % Close peak
                Peak = CloseFlipPeak(Peak, f, k, idx, rdir, pltFLAG);
        
                % End loop
                break
            else
                % Shift loop indexer
                k = k + sval;
            end
        end
    
        function lbd = ReBound(lbd, ubd, iMIN)
            % Check if index is encroaching on a previous peak
            if lbd <= iMIN
                % Shift index back
                lbd = iMIN + 1;
            end
    
            % Check if index needs to be padded back
            if lbd >= ubd
                lbd = ubd;
            end
        end
        % function i = ReBound(i, sval, rdir, iMAX)
        %     % Check if index is encroaching on a previous peak
        %     if ~rdir && ~fFLAG) || (rdir && fFLAG)
        %         % If moving left with not flipped data or moving right with flipped data
        %         if i <= iMAX && ~fFLAG) || (i >= iMAX && fFLAG)
        %             % Shift index back
        %             i = iMAX - sval*2;
        %         end
        %     end
        % end
    
        function Peak = CloseFlipPeak(Peak, f, i, idx, rdir, pltFLAG)
            % % Correct flipped data
            % if fFLAG
            %     f = flip(f);
            %     Peak(idx).Peak(1) = length(f) - Peak(idx).Peak(1) + 1;
            %     i = length(f) - i + 1;
            %     if ~isempty(Peak(idx).StrikesSTD); Peak(idx).StrikesSTD = length(f) - Peak(idx).StrikesSTD + 1; end
            %     if ~isempty(Peak(idx).StrikesSign); Peak(idx).StrikesSign = length(f) - Peak(idx).StrikesSign + 1; end
            %     if ~isempty(Peak(idx).Reopened); Peak(idx).Reopened = length(f) - Peak(idx).Reopened + 1; end
            % end
    
            % Add current point to list
            if ~rdir % && ~fFLAG) || (rdir && fFLAG) % moving right w/ flipped data or moving left w/ unflipped data
                Peak(idx).Start(1) = i;
            else % moving left w/ flipped data or moving right w/ unflipped data
                Peak(idx).Stop(1) = i;
            end
    
            if pltFLAG
                % Plot closed peak
                PlotPoints(f(i), i, idx)
            end
        end
    end
    
    function ib = BuffIndex(i, ss, b)
        if round(ss/b) == 0
            ib = (i-1):(i+1);
        else
            ib = (i-round(ss/b)):(i+round(ss/b));
        end
    end
    
    function Pk = CalculatePeakInfo(Pk, I, f, df, dfsigns)
        for k = I
            % Calculate full width @ base
            Pk(k).Width = Pk(k).Stop(1) - Pk(k).Start(1) + 1;
        
            % Find true peak
            Pk(k).Peak(1) = find(f == max(f(Pk(k).Start(1):Pk(k).Stop(1))));
        
            % Measure amplitude
            base = (Pk(ii).Stop(2) + Pk(ii).Start(2))/2;
            Pk(ii).Amplitude = Pk(ii).Peak(2) - base;

            % Measure total amount of strikes
            Pk(k).TotalStrikes = length(Pk(k).StrikesSTD) + length(Pk(k).StrikesSign);
        
            % Calculate halves
            if Pk(k).Start(1) ~= Pk(k).Peak(1)
                Pk(k).LowerHalf = struct('Length', Pk(k).Peak(1) - Pk(k).Start(1), ...
                                            'Sum', sum(df(Pk(k).Start(1):Pk(k).Peak(1))), ...
                                           'Sign', sum(dfsigns(Pk(k).Start(1):Pk(k).Peak(1))));
            else
                Pk(k).LowerHalf = struct('Length', NaN, ...
                                           'Sign', NaN);
            end
            if Pk(k).Peak(1) ~= Pk(k).Stop(1)+1 % Plus one because peak point is included as part of second half
                Pk(k).UpperHalf = struct('Length', Pk(k).Stop(1)+1 - Pk(k).Peak(1), ...
                                            'Sum', sum(df(Pk(k).Peak(1):Pk(k).Stop(1)+1)), ...
                                           'Sign', sum(dfsigns(Pk(k).Peak(1):Pk(k).Stop(1)+1)));
            else
                Pk(k).UpperHalf = struct('Length', NaN, ...
                                           'Sign', NaN);
            end
        end
    end
    
    function PlotPoints(f, i, idx, n)
        % Check if need to remove points
        if nargin == 4
            pltdata = findobj(gca, 'Type', 'Line');
            for pt = n:-1:1
                delete(pltdata(pt));
            end
        end
        
        % Plot points
        if rem(idx,2) == 0
            plot(i, f, 'r.', 'MarkerSize', 20)
        else
            plot(i, f, 'b.', 'MarkerSize', 20)
        end
    end
    
    function [spec, filtsz] = FilterSpec(spec0, BW, fmethod, dtype, F0)
        % Inputs:
        %   spec0: complex data
        %      BW: real, imag, mag, or phase data (depends on what you want to smooth)
        %  method: filter or smooth
        %   dtype: what data you want outputted, (i.e. real, imaginary, magnitude, phase, or complex). Smooth cannot output complex data
        %      F0: center frequency
        
        % Isolate real data
        spec_re = real(spec0);
        
        switch fmethod
            case 'filter'
                % Calculate SNR
                S = max(spec_re);
                N = std(spec_re(end-length(spec_re)/8:end));
                snr = S/N;
        
                if exist('filtsz', 'var')
                    clear filtsz
                end
        
                % Select filter size based on SNR
                if snr > 500
                    filtsz = 2;
                elseif snr > 250
                    filt_ub = 500;
                    filt_lb = 230;
                    filt_minsz = 2;
                    filt_maxsz = 3;
                elseif snr < 75
                    filtsz = 10;
                else
                    if snr > 100
                        filt_ub = 250;
                        filt_lb = 100;
                        filt_minsz = 2;
                        filt_maxsz = 5;
                    elseif snr > 75
                        filt_ub = 100;
                        filt_lb = 75;
                        filt_minsz = 5;
                        filt_maxsz = 10;
                    end
                end
                if ~exist('filtsz', 'var')
                    filtsz = round(((filt_ub - snr)/(filt_ub - filt_lb))*(filt_maxsz - filt_minsz) + filt_minsz);
                end
        
                % Zero pad spectra (to account for effect of filter on sides)
                pval = length(spec0)/2;
                spec_corr_pad = padarray(spec0.', pval).';
        
                % Extract corrected fids
                fids_corrected = ifft(ifftshift(spec_corr_pad));
        
                % Temporarily filter (and flip) data
                spec_corrected = fftshift(fft(exp(-(pi*filtsz).*((1:length(fids_corrected))/(BW*2))).*fids_corrected)).';
        
                % Isolate final data
                switch dtype
                    case 'Real'
                        spec = real(spec_corrected);
                    case 'Imaginary'
                        spec = imag(spec_corrected);
                    case 'Magnitude'
                        spec = abs(spec_corrected);
                    case 'Phase'
                        spec = angle(spec_corrected);
                    case 'Complex'
                        spec = spec_corrected;
                end
        
            case 'smooth'
                % Calculate linewidth (# pts)
                LW_pts = round(Linewidth_calculator(spec_re,F0,BW)*(length(spec_re)/BW));
        
                % Isolate data type
                switch dtype
                    case 'Real'
                        spec = real(spec0);
                    case 'Imaginary'
                        spec = imag(spec0);
                    case 'Magnitude'
                        spec = abs(spec0);
                    case 'Phase'
                        spec = angle(spec0);
                end
        
                % Pad data to remove inconsistencies at edges
                pval = 2;
                spec = padarray(spec, [0 pval], mean(spec(1:round(length(spec0)/40))), "pre");
                spec = padarray(spec, [0 pval], mean(spec(end-round(length(spec0)/40):end)), "post");
        
                % Smooth data
                filtsz = LW_pts;
                spec = smooth(spec, filtsz,'moving');
        end
        
        % Trim off zero padding
        spec = spec(pval+1:end-pval).';
    end
    
    function FWHM = Linewidth_calculator(data, F0, BW)
        % Calculate x-axis in Hz
        faxis = (linspace(-BW/2, BW/2, length(data)) + F0*1e6)/1e6;
        
        % Find peak
        [peak_max, peak_pos] = max(data);
        
        % Calculate half max
        HM = peak_max/2;
        
        % Initialize data & flags
        flag_FWHM = false;
        flag_RHM = false;
        flag_LHM = false;
        RHM = 0;
        LHM = 0;
        steps = 1;
        
        % run continuously until flag_FWHM is true
        while ~flag_FWHM 
            if peak_pos-steps > 1 || peak_pos+steps < length(data) && HM > 0 % If locations are within 1 and less than the length of the provided vector
                if ~flag_LHM && peak_pos-steps > 1
                    % Identifies peak height of point(s) descending from peak value on the left side for half-height
                    test_LHM = data(peak_pos-steps);
                else % z == 0 && peak_pos-t < 1
                    test_LHM = HM;
                end
        
                if ~flag_RHM && peak_pos+steps < length(data)
                    % Identifies peak height of point descending from peak value on the right side for half-height
                    test_RHM = data(peak_pos+steps);
                else % y == 0 && peak_pos+t > length(X)
                    test_RHM = HM;
                end
        
                % If test_left point is less than the half-maximum value, then interpolate for left half-max position
                if ~flag_LHM && test_LHM < HM
                    x_points = [faxis(peak_pos-steps) faxis(peak_pos-steps+1)];
                    y_points = [test_LHM data(peak_pos-steps+1)];
                    a = polyfit(x_points,y_points,1);
                    LHM = (HM-a(2))/a(1);
                    flag_LHM = true;
                end
        
                % If test_right point is less than the half-maximum value, then interpolate for right half-max position
                if ~flag_RHM && test_RHM < HM
                    x_points = [faxis(peak_pos+steps-1) faxis(peak_pos+steps)];
                    y_points = [data(peak_pos+steps-1) test_RHM];
                    a1 = polyfit(x_points,y_points,1);
                    RHM = (HM-a1(2))/a1(1);
                    flag_RHM = true;
                end
        
                % If positions are within bounds, assign values and stop while loop.
                if LHM ~= 0 && RHM ~= 0
                    FWHM = (RHM - LHM)*10^6;
                    flag_FWHM = true;
                end
        
                steps = steps+1;  % Increase value of steps or position away from peak with each iteration of while loop
        
            else  % No peak is found within the region
                flag_FWHM = true;
                FWHM = 0;
            end
        end
    end
end
