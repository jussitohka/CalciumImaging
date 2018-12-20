% finding the group code and file index based on filename
% for schizophrenia dataset 
% code: group code
% index: day-wise index
% date: measurement date
% indexDate: an unique index to each dish
function [code,index,dateCode,indexDate] = neurotrans_fn2code_sz(fn)
   nfn = length(fn);
   code = zeros(nfn,1);
   index = zeros(nfn,1);
   dateCode.month = zeros(nfn,1);
   dateCode.day = zeros(nfn,1); 
   indexDate = zeros(nfn,1);
   for i = 1:nfn
       idxus = strfind(fn{i},'_');
       idxus = idxus(1);
      
       idxsz = strfind(fn{i},'SZ');
       idxsz = idxsz(1);
       code(i) = str2num(fn{i}((idxsz +2):(idxus - 1)));
       index(i) = str2num(fn{i}(1:(idxsz -1)));
       if ~isempty(strfind(fn{i},'Dec')) | ~isempty(strfind(fn{i},'dec'))
           dateCode.month(i) = 12; % December
           idxmo = strfind(fn{i},'ec');
       else
           dateCode.month(i) = 11; % November
           idxmo = strfind(fn{i},'ov');
       end
       
       dateCode.day(i) = str2num(fn{i}((idxus + 1):(idxmo - 2)));
   end
   dc = dateCode.month*100 +  dateCode.day;
   udc = unique(dc);
   sdc = zeros(nfn,1);
   for i = 1:nfn
       sdc(i) = find(udc == dc(i)); % codes 1 to 6 according to the measurement date
   end
   indexDate = sdc*100 + index;
   
   
end

