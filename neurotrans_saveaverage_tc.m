function [ucode,N] = neurotrans_saveaverage_tc(xlsfn,ts,code,index,xtime,treated)

if ~exist('treated','var')
    treated{1} = [1:3 7:9]; % untreated cells
    treated{2} = [4:6 10:12]; % treated cells
end



ucode = unique(code);
ttt = 1;
f{ttt} = 'Time(s)'
ttt = 2;
for i = 1:length(ucode)
    f{ttt} = ['SZ',num2str(ucode(i)),'untreated'];
    f{ttt + 1} = ['SZ',num2str(ucode(i)),'treated']; 
    ttt = ttt + 2;
end
xlswrite(xlsfn,f,1,'A1');

mats = [xtime];
for i = 1:length(ucode)
    for k = 1:2
        idx = find(code == ucode(i) & ismember(index,treated{k}));
        ats = ts{idx(1)};
        for j = 2:length(idx)
            ats = [ats ts{idx(j)}];
        end
        mats = [mats mean(ats,2)];
        N(i,k) = size(ats,2);
    end
end
xlswrite(xlsfn,mats,1,'A2');