function [ucode,N] = neurotrans_saveaverage_tc(xlsfn,ts,code,index,xtime,treated)

if ~exist('treated','var')
    treated{1} = [1:3 7:9]; % untreated cells, for compatibility 
    treated{2} = [4:6 10:12]; % treated cells, for compatibility
end

Ntc = 10;
trstr{1} = 'untreated';
trstr{2} = 'treated';

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
clear f
mats = [xtime];
for i = 1:length(ucode)
    ats10 = [xtime];
    f{1} = 'Time(s)'; 
    for k = 1:2
        
        idx = find(code == ucode(i) & ismember(index,treated{k}));
        ats = ts{idx(1)};
        for j = 2:length(idx)
            ats = [ats ts{idx(j)}];
        end
        ats10 = [ats10 mean(ats,2)];
        ridx = randperm(size(ats,2),Ntc);
        mats = [mats mean(ats,2)];
        N(i,k) = size(ats,2);
        f{2 + (k -1)*(Ntc +1)} = ['Mean ',trstr{k}];
        
        for j = 1:Ntc
            f{j + 2 + (k - 1)*(Ntc +1)} = ['SZ',num2str(ucode(i)),trstr{k},num2str(ridx(j))];
        end
        ats10 = [ats10 ats(:,ridx)];
    end
    xlsfn10 = [xlsfn(1:(end - 4)),'SZ',num2str(ucode(i)),'.xls'];
    xlswrite(xlsfn10,f,1,'A1');
    xlswrite(xlsfn10,ats10,1,'A2');
end
xlswrite(xlsfn,mats,1,'A2');