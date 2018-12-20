% finding the group code and file index based on filename
% for schizophrenia dataset 
% code: group code
% index: day-wise index
% date: measurement date
% indexDate: an unique index to each dish
function [code,index,dateCode,indexDate] = neurotrans_fn2code_sz_males(fn)
   nfn = length(fn);
   code = zeros(nfn,1);
   index = zeros(nfn,1);
   dateCode.month = zeros(nfn,1);
   dateCode.day = zeros(nfn,1);
   dateCode.year = zeros(nfn,1);
   indexDate = zeros(nfn,1);
   for i = 1:nfn
       idxus = strfind(fn{i},'_');
       idxus = idxus(1);
      
       idxsz = strfind(fn{i},'SZ');
       idxsz = idxsz(1);
       code(i) = str2num(fn{i}((idxsz +2):(idxus - 1)));
       index(i) = str2num(fn{i}(1:(idxsz -1)));
       dateCode.year(i) = str2num(fn{i}((idxus +1):(idxus + 4)));
       tmpstr = fn{i}((idxus +6):(idxus + 7));
       if strcmp(tmpstr(2),'.')
          dateCode.month(i) = str2num(tmpstr(1:1));
          incr = 1;
       else
          dateCode.month(i) = str2num(tmpstr(1:2));
          incr = 2; 
       end
       tmpstr = fn{i}((idxus +6 + incr +1):(idxus + 6 + incr + 2));
       if strcmp(tmpstr(2),'.')
           dateCode.day(i) = str2num(tmpstr(1:1));
       else
           dateCode.day(i) = str2num(tmpstr(1:2));
       end
       
   end
   dc = dateCode.month*100 +  dateCode.day;
%    udc = unique(dc);
%    sdc = zeros(nfn,1);
%    for i = 1:nfn
%        sdc(i) = find(udc == dc(i)); % codes 1 to 6 according to the measurement date
%    end
   indexDate = dc*100 + index;
   
   
end

