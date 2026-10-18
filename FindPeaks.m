function [TruePeaks, filtsz, ss] = FindPeaks(f_comp, F0, BW, opts)
%% Initialize
[f_comp, LW, SNR, filtsz, debugFLAG, A0, A1, PKmin, PKmax] = parseInputs(f_comp, opts, F0, BW);

% Determine x-axis step size in terms of number of points
Faxis_ss = BW/length(f_comp); % Hz/pt

% Calculate number of scales for differentiation
scales = LW/(Faxis_ss*8);
if scales < 3, scales = 3; end

% Plot data
if debugFLAG
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
        if i == 1; legend('f', 'fSTD', 'df', 'dfAD', 'dfSTD'); end; if i == 4; ylim([-2*pi 2*pi]); end
        clear f fAD fSTD df dfAD dfSTD
    end
    f_og = real(f_comp); %#ok<NASGU>
    clear ptitle pdata i
end

% Filter data
[F, filtsz] = FilterSpec(f_comp, BW, 'filter', struct('SNR', SNR, 'filtsz', filtsz));

% % Apply phasing
% Axis_for_phasing = (-BW/2):BW/length(f_comp):(BW/2 - 1/length(f_comp));
% F = F.*exp(1j*((A1*2*pi*Axis_for_phasing)+(A0*pi/180)));

f_mg = abs(F);
f = real(F);

% Determine linewidth of main peak (in Hz)
LW2 = Linewidth_calculator(f, F0, BW);

% Determine sample size in terms of number of points (we will use LW as sample size, can be changed later)
ss = ceil(LW2/Faxis_ss);

% Calculate statistical datasets
[fSTD, df, dfAD, dfSTD, dfsigns] = PeakStats(f, scales);

% Obtain first standard deviation line for phase
df_ph = cwt_custom(angle(f_comp), scales);
dfSTD_ph = mean(df_ph) + std(df_ph);

%% Find easy peaks in magnitude & real data
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
if ~(PK >= PKmin && PK <= PKmax)
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
            if PK <= PKmin % Increase
                dw = 1/2; dp = 1/2;
            else % Decrease
                dw = 2; dp = 2;
            end
        else
            % Change
            if PK <= PKmin % Increase
                dw = dw/2; dp = dp/2;
            else % Decrease
                dw = dw*2; dp = dp*2;
            end
        end

        % Change width
        w1 = w*dw; % New value
        tmp_PKw = FindMgRePeaks(f_mg, f, w1, p, ss);
        tmp_PKw = length(tmp_PKw);

        % Change prominence
        p1 = p*dp; % New value
        tmp_PKp = FindMgRePeaks(f_mg, f, w, p1, ss);
        tmp_PKp = length(tmp_PKp);

        % Check for a change with either
        if PK ~= tmp_PKw || PK ~= tmp_PKp
            % Set flag
            delpkFLAG = true;

            % Determine what to look for
            if PK <= PKmin % need to increase
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
            if PK >= PKmin && PK <= PKmax % no change needed
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
    MPlocation = FindMgRePeaks(f_mg, f, w, p, ss);
    clear w1 p1 w p w0 p0 dw dp iter delpkFLAG tmp_PKw tmp_PKp
else
    MPlocation = FindMgRePeaks(f_mg, f, w0, p0, ss);
    clear w0 p0
end
clear df_tmp dfSTD_tmp

% Remove points to close to edges
MPlocation = MPlocation(MPlocation>ss*2 & MPlocation<(length(f)-ss*2));
PK = length(MPlocation);

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
MPphasestd = zeros(PK,1);
pFLAG = true(PK,1);
for i = 1:PK
    % Set region
    s = BuffIndex(MPlocation(i), ss/2);

    % Find phase STD
    MPphasestd(i) = std(df_ph(s));

    % Check STD of phase df & AD of real df
    % - STD greater than 1/3*global but less than 2*global is likely noise
    % - some peaks will be off a bit in magnitude data, so check if AD of at least half the points are less then dfSTD
    if MPphasestd(i) < ub
        if MPphasestd(i) > dfSTD_ph/3 && sum(dfAD(s)>dfSTD) < ss/2
            % Measure phase STD of neighboring regions
            pSTDl = std(df_ph(BuffIndex(MPlocation(i)-ss*2, ss)));
            pSTDu = std(df_ph(BuffIndex(MPlocation(i)+ss*2, ss)));
    
            % Check STD of phase df of real df (but bounds*1.5)
            % - smaller peaks in noisy data often aren't caught by first check
            % - AD is especially problematic with noisy data
            if MPphasestd(i) < dfSTD_ph/2
                % Check STD of phase df against that of neighboring region
                % - ensure region is not just a noise region
                if (MPphasestd(i) > pSTDl/3 && pSTDl < ub) && ...
                   (MPphasestd(i) > pSTDu/3 && pSTDu < ub)
                    pFLAG(i) = false;
                end
            else
                pFLAG(i) = false;
            end
        end
    end
end
lc = MPlocation(pFLAG);

%% Now find locations in real data (might be off a bit due to smoothing & peak broadening)
% Plot data
if debugFLAG, PlotPeaks('All Peaks', f, df, dfsigns, fSTD, dfSTD, dfAD); end

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
    ilc = find(f == max(f(BuffIndex(i,ss/4))));

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
    Peaks0 = SearchLimit(Peaks0, idx, ss, f, fSTD, dfAD, dfSTD, dfsigns, false, debugFLAG);

    % Transfer strikes & reset in structure
    tmp_sSTD = Peaks0(idx).StrikesSTD;
    tmp_sSign = Peaks0(idx).StrikesSign;
    tmp_sOpen = Peaks0(idx).Reopened;
    Peaks0(idx).StrikesSTD = [];
    Peaks0(idx).StrikesSign = [];
    Peaks0(idx).Reopened = [];

    % Then, move right to find end
    Peaks0 = SearchLimit(Peaks0, idx, ss, f, fSTD, dfAD, dfSTD, dfsigns, true, debugFLAG);

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
        % Remove from plot
        if debugFLAG
            plotdata = findobj(gca, 'Type', 'Line');
            for j = 2:-1:1; delete(plotdata(j)); end
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
        % Plot peak
        if debugFLAG
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

%% Combine doublets/triplets
% Plot
if debugFLAG, PlotPeaks('Combined Peaks', f, df, dfsigns); end

% Measure difference in height between left/peak, right/peak, & left/right
Ptdif = ones(length(Peaks0), 4);
for i = 1:length(Peaks0)
    Ptdif(i,1) = Peaks0(i).Peak(2) - Peaks0(i).Start(2); % left/lower half
    Ptdif(i,2) = Peaks0(i).Peak(2) - Peaks0(i).Stop(2); % right/upper half
    Ptdif(i,3) = abs(Peaks0(i).Stop(2) - Peaks0(i).Start(2)); % width

    % Normalize
    Ptdif(i,:) = Ptdif(i,:)./max(Ptdif(i,:));

    % Find min
    Ptdif(i,4) = find(Ptdif(i,:) == min(Ptdif(i,:)));
end

% Check for doublets/triplets
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
                J = i-1:-1:1; if length(J) > 3, J = J(1:3); end
                for j = J
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
                J = i+1:length(Peaks0); if length(J) > 3, J = J(1:3); end
                for j = J
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

                % Plot
                if debugFLAG
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
            RLgap = Peaks0(R).Start(1) - Peaks0(L).Stop(1);

            % Check for no doublets/triplets
            if np == 0
                % Check if overlapping peaks
                % - the case if height difference between shared start & stop is small
                if abs(Peaks0(L).Start(2)-Peaks0(R).Start(2)) < (Peaks0(L).Peak(2)-Peaks0(L).Start(2))/2 && ...
                   abs(Peaks0(L).Start(2)-Peaks0(R).Start(2)) < (Peaks0(R).Peak(2)-Peaks0(R).Stop(2))/2
                    % Set flag & sample
                    if strcmp(dirFLAG, 'left')
                        s = i-1:i;
                    else
                        s = i:i+1;
                    end
                    combFLAG(s) = true;
                end

            % Check if the shared halves of first & last peaks are 1/2 as long as unshared halves & distance between peak start & ends is less than ss
            elseif Peaks0(R).LowerHalf.Length < Peaks0(R).UpperHalf.Length/2 && ...
                   Peaks0(L).UpperHalf.Length < Peaks0(L).LowerHalf.Length/2 && ...
                   RLgap < ss
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
                % Remove from plot
                if debugFLAG
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
            Peaks1(idx).LowerHalf = Peaks0(s(1)).LowerHalf;
            Peaks1(idx).UpperHalf = Peaks0(s(end)).UpperHalf;

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

        % Plot
        if debugFLAG
            PlotPoints(Peaks1(idx).Start(2), Peaks1(idx).Start(1), idx)
            PlotPoints(Peaks1(idx).Stop(2), Peaks1(idx).Stop(1), idx)
        end
    end
end

%% Clean up peaks
% Plot
if debugFLAG, PlotPeaks('Cleaned Peaks', f, df, dfsigns); end

% Clean up peaks & get more information for analysis
idx = 0;
for i = 1:length(Peaks1)
    % Reset first if FLAG
    newpeak = false;

    % Check if one side is not >> than other side
    if Peaks1(i).UpperHalf.Length/Peaks1(i).LowerHalf.Length < 5 || Peaks1(i).LowerHalf.Length/Peaks1(i).UpperHalf.Length < 5
        % Reset flags
        leftFLAG = false;
        rightFLAG = false;

        % Check if peak is a true peak (lower half should be positive slope)
        if Peaks1(i).LowerHalf.Sign >= Peaks1(i).LowerHalf.Length*0.75
            leftFLAG = true;
        % Some allowence is made for doublet/triplet peaks & poor shim
        elseif Peaks1(i).LowerHalf.Sign >= Peaks1(i).LowerHalf.Length*0.675
            % Create normalized peak value
            Pkmain = Peaks1(i).Peak(2) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Find max value to left of peak
            Pkrange = max(f(round(Peaks1(i).Start(1)+ss/4):round(Peaks1(i).Peak(1)-ss/4))) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Check for possible doublet (if normalized point is 90% of main peak)
            if Pkrange/Pkmain >= 0.9, leftFLAG = true; end
        end

        % Check if peak is a true peak (upper half should be negative slope)
        if Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25
            rightFLAG = true;

        % Some allowence is made for doublet/triplet peaks & poor shim
        elseif Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.325
            % Create normalized peak value
            Pkmain = Peaks1(i).Peak(2) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Find max value to right of peak
            Pkrange = max(f(round(Peaks1(i).Peak(1)+ss/4):round(Peaks1(i).Stop(1)-ss/4))) - (Peaks1(i).Start(2) + Peaks1(i).Stop(2))/2;

            % Check for possible doublet (if normalized point is 90% of main peak)
            if Pkrange/Pkmain >= 0.9, rightFLAG = true; end
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
    if idx > 0 && ~newpeak && ~previouspeakcombination &&...
       Peaks1(i).UpperHalf.Sign <= Peaks1(i).UpperHalf.Length*0.25 && ...
       sum(dfsigns(Peaks2(idx).Stop(1)+1:Peaks1(i).Start(1))) >= (Peaks1(i).Start(1)-(Peaks2(idx).Stop(1)+1))*0.75
        
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
end

%% Combine overlapping peaks
% Clean up peaks & get more information for analysis
if length(Peaks2) > 1
    % Plot
    if debugFLAG
        PlotPeaks('Final Peaks', f, df, dfsigns);
    end

    % Initialize indexers
    i = 0;
    idx = 0;
    idx1 = 1;

    % Loop through peaks
    while i < length(Peaks2)
        % Shift indexer
        i = i + 1;

        % Reset flag
        combFLAG = false;

        % Check if there is a peak to the left
        if idx > 0
            % Check if peak ends close to start of current peak
            if Peaks2(i).Start(1) - Peaks3(idx).Stop(1) < ss/3
                % Check if left is significantly shorter & closer to peak than right
                if Peaks2(i).UpperHalf.Length/Peaks2(i).LowerHalf.Length >= 3 && ...
                   Peaks2(i).Peak(2)-Peaks2(i).Start(2) < 3*(Peaks2(i).Peak(2)-Peaks2(i).Stop(2))
                    % Update flag
                    combFLAG = true;

                    % Plot point
                    if debugFLAG
                        PlotPoints(Peaks2(i).Stop(2), Peaks2(i).Stop(1), idx, 1)
                    end
                end
            end

        % Check if there is a peak to the right
        elseif i < length(Peaks2)
            % Check if peak ends close to start of current peak
            if Peaks2(i+1).Start(1) - Peaks2(i).Stop(1) < ss/3
                % Check if right is significantly shorter & closer to peak than left
                if Peaks2(i).LowerHalf.Length/Peaks2(i).UpperHalf.Length >= 2.5 && ...
                  (Peaks2(i).Peak(2)-Peaks2(i).Stop(2)) < 3*(Peaks2(i).Peak(2)-Peaks2(i).Start(2))
                    % Update indexers & flag
                    i = i + 1;
                    idx = idx + 1;
                    combFLAG = true;

                    % Index peak
                    Peaks3(idx) = Peaks2(i-1);

                    % Plot point
                    if debugFLAG
                        PlotPoints(Peaks3(idx).Start(2), Peaks3(idx).Start(1), idx)
                        PlotPoints(Peaks2(i).Stop(2), Peaks2(i).Stop(1), idx)
                    end
                end
            end
        end

        % Check if combined
        if combFLAG
            % Update indexer
            idx1 = idx1 + 1;

            % Combine data
            Peaks3(idx).Stop = Peaks2(i).Stop;
            Peaks3(idx).Width = Peaks3(idx).Stop(1) - Peaks3(idx).Start(1);
            Peaks3(idx).Peak(idx1:(idx1+size(Peaks2(i).Peak,1)-1),:) = Peaks2(i).Peak;
            Peaks3(idx).Amplitude(idx1:(idx1+length(Peaks2(i).Amplitude)-1)) = Peaks2(i).Amplitude;
            Peaks3(idx).UpperHalf = Peaks2(i).UpperHalf;
            Peaks3(idx).StrikesSTD = [Peaks3(idx).StrikesSTD, Peaks2(i).StrikesSTD];
            Peaks3(idx).StrikesSign = [Peaks3(idx).StrikesSign, Peaks2(i).StrikesSign];

        else
            % Reset subindexer
            idx1 = 1;

            % Append to array
            idx = idx + 1;
            Peaks3(idx) = Peaks2(i);

            % Plot peak
            if debugFLAG
                PlotPoints(Peaks3(idx).Start(2), Peaks3(idx).Start(1), idx)
                PlotPoints(Peaks3(idx).Stop(2), Peaks3(idx).Stop(1), idx)
            end
        end
    end
else
    Peaks3 = Peaks2;
end

%% Final Verifications
% % Shift starts & stops to *very* local mins
% for idx = 1:length(Peaks3)
%     f_tmp = NaN(size(f_comp));
%     f_tmp(BuffIndex(Peaks3(idx).Start(1),ss,8)) = dfAD(BuffIndex(Peaks3(idx).Start(1),ss,8));
%     [~, Peaks3(idx).Start(1)] = min(f_tmp, [], 'all', 'omitnan');
%     f_tmp = NaN(size(f_comp));
%     f_tmp(BuffIndex(Peaks3(idx).Stop(1),ss,8)) = dfAD(BuffIndex(Peaks3(idx).Stop(1),ss,8));
%     [~, Peaks3(idx).Stop(1)] = min(f_tmp, [], 'all', 'omitnan');
% end
% 
% % Check if start & stop of neighboring peaks are flipped
% for i = 2:length(Peaks3)
%     if Peaks3(i).Start(1) < Peaks3(i-1).Stop(1)
%         % Flip points
%         truestart = Peaks3(i-1).Stop(1);
%         trueend = Peaks3(i).Start(1);
%         Peaks3(i-1).Stop(1) = trueend;
%         Peaks3(i).Start(1) = truestart;
% 
%         % Recalculate width
%         Peaks3 = CalculatePeakInfo(Peaks3, i, f, df, dfsigns);
%         Peaks3 = CalculatePeakInfo(Peaks3, i-1, f, df, dfsigns);
%     end
% end
% 
% % Check for overlapping peaks
% keeps = true(length(Peaks3),1);
% for i = 2:length(Peaks3)
%     if  Peaks3(i).Start(1) < Peaks3(i-1).Stop(1)
%         % Combine peaks
%         Peaks3(i).Start(1) = Peaks3(i-1).Start(1);
% 
%         % Flag old peak
%         keeps(i-1) = false;
%     end
% end
% 
% % Remove old peaks
% TruePeaks = Peaks3(keeps);
TruePeaks = Peaks3;

%% Functions

    function [fcomp, LW, SNR, filtsz, debugFLAG, a0, a1, PKmin, PKmax] = parseInputs(fcomp, opts, F0, BW)
        % Ensure input is structure
        if ~isstruct(opts), error('Options variable must be a structure'); end

        % Check input
        if length(size(fcomp)) > 2
            error('Input must be a 1D array')
        elseif size(fcomp,1) > 1
            if size(fcomp,2) > 1
                error('Input must be a 1D array')
            else
                fcomp = fcomp.';
            end
        end
        freal = real(fcomp);

        % Set default parameters
        filtsz = 100;
        debugFLAG = false;
        a1 = 0;
        a0 = 0;
        PKmin = 10;
        PKmax = 30;

        % Parse through input structure
        if isfield(opts, 'LW')
            LW = opts.LW;
        else
            LW = Linewidth_calculator(freal, F0, BW);
        end
        if isfield(opts, 'SNR')
            SNR = opts.SNR;
        else
            SNR = max(freal)/std(freal(end-length(freal)/8:end));
        end
        if isfield(opts, 'filtsz'), filtsz = opts.filtsz; end
        if isfield(opts, 'debugFLAG'), debugFLAG = opts.debugFLAG; end
        if isfield(opts, 'a1'), a1 = opts.a1; end
        if isfield(opts, 'a0'), a0 = opts.a0; end
        if isfield(opts, 'PKmin'), PKmin = opts.PKmin; end
        if isfield(opts, 'PKmax'), PKmax = opts.PKmax; end
    end

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
           
    function LC = FindMgRePeaks(f_mg, f_re, w, p, ss)
        % Find peaks in magnitude & real data
        [~, LC_mg] = findpeaks(f_mg, 'MinPeakWidth', w, 'MinPeakProminence', p, 'MinPeakDistance', ss);
        [~, LC_re] = findpeaks(f_re, 'MinPeakWidth', w, 'MinPeakProminence', p, 'MinPeakDistance', ss);

        % Combine peaks
        LC = sort([LC_mg, LC_re]);

        % Combine peaks that are too close to each other
        ix = 0;
        while ix < length(LC)-1
            ix = ix + 1;
            if LC(ix+1)-LC(ix) < ss
                LC(ix) = round((LC(ix+1)+LC(ix))/2);
                LC(ix+1) = 0;
                ix = ix + 1;
            end
        end
        LC = LC(LC~=0);
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
                elseif any(dfAD(BuffIndex(k,ss/8)) >= 2*dfSTD)
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
                elseif any(dfAD(BuffIndex(k,ss/4)) >= dfSTD)
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
                if dfsigns(k)
                    closeFLAG = true;
    
                % Check if last 10 points had STD strike
                elseif length(Peak(idx).StrikesSTD) >= 10
                    if all(Peak(idx).StrikesSTD(end-9:end) == k-sval*9:sval:k)
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
    
        function Peak = CloseFlipPeak(Peak, f, i, idx, rdir, pltFLAG)
            % Add current point to list
            if ~rdir
                Peak(idx).Start(1) = i;
            else
                Peak(idx).Stop(1) = i;
            end
    
            if pltFLAG
                % Plot closed peak
                PlotPoints(f(i), i, idx)
            end
        end
    end
    
    function ib = BuffIndex(i, ss)
        ss = round(ss);
        if ss == 0
            ib = (i-1):(i+1);
        else
            ib = (i-ss):(i+ss);
        end
    end
    
    function Pk = CalculatePeakInfo(Pk, I, f, df, dfsigns)
        for k = I
            % Calculate full width @ base
            Pk(k).Width = Pk(k).Stop(1) - Pk(k).Start(1) + 1;
        
            % Find true peak
            Pk(k).Peak(1) = find(f == max(f(Pk(k).Start(1):Pk(k).Stop(1))));
        
            % Measure amplitude
            base = (Pk(k).Stop(2) + Pk(k).Start(2))/2;
            Pk(k).Amplitude = Pk(k).Peak(2) - base;

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
    
    function PlotPeaks(ptitle, f, df, dfsigns, fSTD, dfSTD, dfAD)
        figure; plot(1:length(f), f, 'Color', [0 .447 .741], 'LineWidth', 2); axis tight; hold on;
        if nargin == 6
            yline(fSTD, 'Color', [0 .447 .741], 'LineStyle', '--');
            plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2);
            plot(1:length(f), dfAD.*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184]);
            yline(dfSTD*((max(f)/2)/max(dfAD)), 'Color', [.635 .078 .184], 'LineStyle', '--');
            Ls = {'f', 'fSTD', 'df', 'dfAD', 'dfSTD', 'df Sign'};
        else
            plot(1:length(f), df, 'Color', [.635 .078 .184], 'LineWidth', 2);
            Ls = {'f', 'df', 'df Sign'};
        end
        plot(1:length(f), dfsigns*(max(df)/2), 'Color', [.929 .694 .125]);
        LL = legend(Ls); LL.AutoUpdate = 'off';
        title(ptitle)
    end

    function [spec, filtsz] = FilterSpec(spec0, BW, fmethod, opts)
        % Inputs:
        %   spec0: complex data
        %      BW: real, imag, mag, or phase data (depends on what you want to smooth)
        %  method: filter or smooth
        %    opts: options structure containing
        %           - F0, only used if method is smooth
        %           - SNR: (optional), only used if method is filter 
        %           - filtsz: (optional), only used if method is filter

        % Isolate real data
        spec_re = real(spec0);

        % Extract options
        switch fmethod
            case 'filter'
                if isfield(opts, 'SNR')
                    snr = opts.SNR;
                else
                    S = max(spec_re);
                    N = std(spec_re(end-length(spec_re)/8:end));
                    snr = S/N;
                end
                if isfield(opts, 'filtsz')
                    filtsz_old = opts.filtsz;
                else
                    filtsz_old = 100;
                end
            case 'smooth'
                if isfield(opts, 'F0')
                    f0 = opts.F0;
                else
                    error('F0 is a required input for smoothing')
                end
        end

        switch fmethod
            case 'filter'
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

                % Overwrite if new filter size is larger than previous
                if filtsz > filtsz_old, filtsz = filtsz_old; end
        
                % Zero pad spectra (to account for effect of filter on sides)
                pval = length(spec0)/2;
                spec_pad = padarray(spec0.', pval).';
        
                % Extract padded fids
                fids_pad = ifft(ifftshift(spec_pad));
        
                % Filter (and flip) data
                spec = fftshift(fft(exp(-(pi*filtsz).*((1:length(fids_pad))/(BW*2))).*fids_pad)).';
        
            case 'smooth'
                % Calculate linewidth (# pts)
                LW_pts = round(Linewidth_calculator(spec_re,f0,BW)*(length(spec_re)/BW));
        
                % Isolate data type
                spec_re = real(spec0);
                spec_im = imag(spec0);
        
                % Pad data to remove inconsistencies at edges
                pval = 2;
                spec_re = padarray(spec_re, [0 pval], mean(spec_re(1:round(length(spec0)/40))), "pre");
                spec_re = padarray(spec_re, [0 pval], mean(spec_re(end-round(length(spec0)/40):end)), "post");
                spec_im = padarray(spec_im, [0 pval], mean(spec_im(1:round(length(spec0)/40))), "pre");
                spec_im = padarray(spec_im, [0 pval], mean(spec_im(end-round(length(spec0)/40):end)), "post");
        
                % Smooth data
                filtsz = LW_pts;
                spec_re = smooth(spec_re, filtsz, 'moving');
                spec_im = smooth(spec_im, filtsz, 'moving');

                spec = spec_re + 1i.*spec_im;
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
        
            else % No peak is found within the region
                flag_FWHM = true;
                FWHM = 0;
            end
        end
    end
end