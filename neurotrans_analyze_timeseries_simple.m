% mdsts is the denoised timeseries which has been normalized with global f0
% and global fmax (ionomycin signal)
% spave gives the mean values between the breakpoints 
% udsts gives the time series which has not been normalized
% this is for computing df/f0 = (f - f0)/f0 based on a local f0
% msts is the original time series normalized by ionomycin signal
% msts is the original time series nomalized by f0
function [mdsts,spave,udsts,msts,mstsf0,uts] = neurotrans_analyze_timeseries_simple(mts,background,idx,bgidx,sptimepoints,completetimepoints,ioind);

plotting = 0;
macro_thr = [10 0.3];
micro_thr = [2.5 0.2];
% remove rows containing NaN's
oldlen = length(idx);
okidx = ones(length(idx),1,'logical');
for i = 1:length(idx) 
    ts = mts(:,idx(i));
    if any(isnan(ts))
        okidx(i) = 0;
    end
end
idx = idx(okidx);
newlen = length(idx);
if oldlen ~= newlen 
    disp('Warning: columns with NaNs: removing them')
end

mdsts = zeros(size(mts,1),length(idx));
msts = zeros(size(mts,1),length(idx));
mstsfo = zeros(size(mts,1),length(idx));
udsts = zeros(size(mts,1),length(idx));

tslen = length(sptimepoints);
if ~exist('ioind','var')  % to keep compatibility 
    ioind = 1690:1730;
end
bgts = mean(background(:,bgidx),2);

dsptp = sptimepoints(2:end) - sptimepoints(1:(end - 1));
breakpointsB = sptimepoints(find(dsptp > 1));
breakpointsE = sptimepoints(find(dsptp > 1) + 1) - 1;
spave = zeros(length(idx),length(breakpointsB));

if length(bgts) < size(mts,1)
    disp('Warning: background time course too short. Correcting for it')
    df = size(mts,1) - length(bgts);
    bgts = [bgts; bgts((end - df + 1):end)];
else
   if length(bgts) > size(mts,1)
       disp('Warning: background time course too long. Correcting for it')
       bgts(1:size(mts,1));
   end
end
% remove rows containing NaN's
okidx = ones(length(idx),1);
for i = 1:length(idx) 
    ts = mts(:,idx(i));
    if any(isnan(ts))
        okidx(i) = 0;
    end
end

        

for i = 1:length(idx) 
    ts = mts(:,idx(i));
    uts(:,i) = ts;
    sts = ts - bgts; % substract the background

    dsts = wden(sts,'modwtsqtwolog','s','mln',3,'haar'); % wavelet denoising
    dsts = dsts';
    udsts(:,i) = dsts;
    [mii min_idx] = min(dsts(sptimepoints));
    min_idx = sptimepoints(min_idx);
    m1 = max(min_idx - 2,1);
    m2 = min(min_idx + 2,max(sptimepoints));
    mi = median(dsts(m1:m2));   
    dsts = dsts - mi; % remove baseline; this has changed 240518; was before min(dsts) 
    ma = max(dsts(ioind));
    dsts = dsts/max(dsts(ioind)); % divide by ion signal
    sts = (sts - mi);
    mstsf0(:,i) = sts/mi;
    msts(:,i) = sts/ma;
    mdsts(:,i) = dsts(completetimepoints);
    for j = 1:length(breakpointsB)
        spave(i,j) = mean(dsts(breakpointsB(j):breakpointsE(j)));
    end    
    
end


