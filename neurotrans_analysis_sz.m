% Data analysis main script, for female cells
% (C) Jussi Tohka 2018 University of Eastern Finland
% see LICENSE for license

close all
clear
if isunix
   addpath /research/users/justoh/matlab/export_fig-master
   load('/research/groups/tohkagroup/neurotrans/sz/SZ_cell_stimulus_amplitudes240518_g_allnorm.mat')
   % load('C:\Users\justoh\Results\rashid\schizo120618\SZ_cell_stimulus_amplitudes240518_g_allnorm.mat')
   figdir = '/research/groups/tohkagroup/neurotrans/schizo120618/figures130618'
   
else
    addpath C:\Users\justoh\matlab\export_fig-master % uses export_fig: https://github.com/altmany/export_fig
    load('C:\Users\justoh\Results\neurotrans\sz\SZ_cell_stimulus_amplitudes240518_g_allnorm.mat')
    figdir = 'C:\Users\justoh\Results\neurotrans\sz\figures'
end
validate = 0; % if this is 1, validates the nested statistical test used (that the alpha level is correct if there is no effect). 
              % If this is 0, then performs the statistical test.  
gtype{1} = 'Global';
gtype{2} = 'Local';
gtype{3} = 'Sp';


gtype2{1} = '';
gtype2{2} = '';
gtype2{3} = 'Sp';

g = g_lnorm;
ggg = 2;
j = 1;

plotfigs = 0;
plotfigs2 = 0;
do_nested = 1;  

remove_outliers = 1;

%between_factors = [2 1 1 2 2 1 3 3 3]';% this is the original ordering
ucode = [1 2 9 10 11 12 14 15 16];
% reorder the data to get groups in the order SZ, TW, CN
g = g(:,:,[2 3 6 1 4 5 9 7 8]); % note the change from the previous version.
celllines = [1 2 3 1 2 3 1 2 3];
between_factors = [1 1 1 2 2 2 3 3 3]'; % This for reordered g
%bf{1} = 'SZ';
%bf{2} = 'twin';
%bf{3} = 'Ctrl';

bf{1} = 'S';
bf{2} = 'T';
bf{3} = 'C';


% find the dish identifiers
for i = 1:9
   ufile{i} = unique(g(:,1,i));
   iii = find(ufile{i} > 0);
   ufile{i} = ufile{i}(iii);
end

 
treated{1} = [1:3 7:9]; % untreated cells
treated{2} = [4:6 10:12]; % treated cells
% SZ1: healthy twin (Non-SZ)   T0
% SZ2: schizophrenic twin (SZ) S0
% SZ15: healthy control (Ctrl) C1 not used in the paper
% 
% SZ9: schizophrenic twin (SZ)  S1
% SZ10: healthy twin (Non-SZ)   T1 
% SZ16: healthy control (Ctrl)  C2 not used in the paper
% 
% SZ11: healthy twin (Non-SZ)    T2
% SZ12: schizophrenic twin (SZ)  S2
% SZ14: healthy control (Ctrl)   C0  not used in the paper



bfadd = zeros(3,1);
% between_facros_str_ind is an unique string for each subj
for i = 1:9
    between_factors_str{i} = bf{between_factors(i)};
    between_factors_str_ind{i} = strcat(bf{between_factors(i)},num2str(bfadd(between_factors(i))));
    bfadd(between_factors(i)) = bfadd(between_factors(i)) + 1;
end


lab = {'FileID'
       'DateCode'
       'DishNo';
       'Glu+Gly (no Mg)'
       'Glu+Gly'
       'GABA'
       'KCl'
       'ionomycin'};
stattype = {'KCL normalized'};
stattype{j} = [gtype{ggg},' ',stattype{j}];





% MAKE DATA MATRICES

for m = 1:1
   for j = 1:3
        datamat{m,j}= zeros(9,2); 
        maddatamat{m,j}= zeros(9,2); 
        for i = 1:9
            kt = 1;
             for k = 1:2
            
                sss = [];
                ttt = [];
                gg = ismember(g(:,3,i),treated{k});            
               
                idxx = find(g(:,1,i) > 0 & gg); 
                sss = [sss; g(idxx,j + 3,i)./g(idxx,7,i)];
    
               
               ttt = [ttt; g(idxx,1,i)];
               if remove_outliers
                   if length(sss) > 5
                       omed = median(sss);
                       omad = 1.4826*mad(sss);
                       inlieridx = find(abs(sss - omed) < 3*omad);
                       sss = sss(inlieridx);
                       ttt = ttt(inlieridx);
                   end
               end
               
                
                datamat{m,j}(i,k) =  median(sss);
                maddatamat{m,j}(i,k) = mad(sss);
                bigdatamat{m,j,i,k} = sss;
                bigdatamat_dishcodes{m,j,i,k} = ttt;
                kt = kt + 1;
            end
        end
        
        
   end
end

if plotfigs2
    nsg = [1 1 2 2 3 3 1 2 3]; % the number of subject in group
    col = {'r.','r.','b.','b.','g.','g.','r.','b.','g.'};
    for m = 1:2
        for j = 1:3
            figure
            hold
           for i = 1:9
               for k = 1:2
                   plot(between_factors(i) + 0.5*(k - 1) + 0.1*nsg(i),bigdatamat{m,j,i,k},col{i})
               end
           end
        end
    end
end



if do_nested
   
    for m = 1:1
        for j = 1:3
            y = [];
            subj = [];
            treat = [];
            dx = [];
            dish = [];
            clear tsubj
            for i = 1:9
                for k = 1:2
                    y = [y; bigdatamat{m,j,i,k}];
                    subj = [subj; i*ones(size(bigdatamat{m,j,i,k}))];
                    treat = [treat; k*ones(size(bigdatamat{m,j,i,k}))];
                    dx  = [dx; between_factors(i)*ones(size(bigdatamat{m,j,i,k}))];
                    dish = [dish;bigdatamat_dishcodes{m,j,i,k}];
                end
            end
            for kkk = 1:length(subj)
                tsubj{kkk} = between_factors_str_ind{subj(kkk)};
            end
			% do also the t-test
			idx = find(subj < 7); % SZs and twins
			if validate 
                zp{m,j} = 0;
                zzz{m,j} = 0;
                iter = 10000;
                for aaa = 1:iter
			         yyy = 3*randn(size(y)) + 1; 
			         [zzzz,zpp] = nested_ttest_twins(yyy(idx),subj(idx),treat(idx),{[1 4],[2 5],[3 6]});
                     zp{m,j} = zp{m,j} + double(zpp < 0.05)/iter;
                     zzz{m,j} = zzz{m,j} + double(zpp < 0.01)/iter;
                end
		    else
			    [zzz{m,j},zp{m,j},tstats{m,j}] = nested_ttest_twins_nochange(y(idx),subj(idx),treat(idx),{[1 4],[2 5],[3 6]},fullfile(figdir,strcat('nested_t_',stattype{m},lab{j + 3})),[lab{j + 3} gtype2{ggg}]);
			end
            
             
        end
    end
end
% ******* END OF MAKE DATA MATRICES FOR ANOVA AND NESTED ANOVA1 ********
if ~validate
    save(fullfile(figdir,'nested_t_stats'),'zzz','tstats')
end





