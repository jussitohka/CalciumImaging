% Data analysis main script, for male cells
% (C) Jussi Tohka 2018 University of Eastern Finland
% see LICENSE for license
close all
clear
if isunix
    addpath /research/users/justoh/matlab/export_fig-master
    load('/research/groups/tohkagroup/neurotrans/sz_males/SZ_cell_stimulus_amplitudes240518_g_allnorm.mat')
    
    figdir = '/research/groups/tohkagroup/neurotrans/sz_males/figures'
   
else
    addpath C:\Users\justoh\matlab\export_fig-master
    load('C:\Users\justoh\Results\neurotrans\sz_males\SZ_cell_stimulus_amplitudes_males200818_g_allnorm.mat')
    figdir = 'C:\Users\justoh\Results\neurotrans\sz_males\figures'
end
validate = 0; % if this is 1, validates the nested statistical test used (that the alpha level is correct if there is no effect). 
              % If this is 0, then performs the statistical test.  
              % validation correct
gtype{1} = 'Global';
gtype{2} = 'Local';
gtype{3} = 'Sp';


gtype2{1} = '';
gtype2{2} = '';
gtype2{3} = 'Sp';
ggg = 2;
g = g_lnorm;


plotfigs = 0;
plotfigs2 = 0;
do_nested = 1;  

remove_outliers = 1;

%between_factors = [2 1 1 2 2 1 3 3 3]';% this is the original ordering
% these are codes for subjects
% line 1 : SZ7 (sick twin)
%          SZ8 (healthy twin)
%          SZ13 (control)
%          SZ20 (control)
% line 2: SZ17 (sic k twin)
%         SZ18 (healthy twin)
%         SZ19 (control)


ucode = [7  8  13  17  18  19  20];
% ucode = [1 2 9 10 11 12 14 15 16];
% reorder the data to get groups in the order SZ, TW, CN
g = g(:,:,[1 4 2 5 3 7 6]); % note the change from the previous version.
celllines = [1 2 1 2 1 1 2];
between_factors = [1 1 2 2 3 3 3]'; % This for reordered g
%bf{1} = 'SZ';
%bf{2} = 'twin';
%bf{3} = 'Ctrl';

bf{1} = 'S';
bf{2} = 'T';
bf{3} = 'C';


% find the dish identifiers
for i = 1:size(g,3)
   ufile{i} = unique(g(:,1,i));
   iii = find(ufile{i} > 0);
   ufile{i} = ufile{i}(iii);
end

 
treated{1} = [1 2]; % untreated cells
treated{2} = [3 4]; % treated cells

bfadd = zeros(3,1);
% between_facros_str_ind is an unique string for each subj
for i = 1:size(g,3)
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
stattype = {'KCL normalized',
            'KCL normalized (log)'
            'Ionom normalized (log)'
            'Ionom normalized'};
        
for j = 1:4
    stattype{j} = [gtype{ggg},' ',stattype{j}];
end
% study correlations  between other and KCL
gsz = size(g,3);
if 0
for i = 1:gsz
   idxx = find(g(:,1,i) > 0 & g(:,4,i) < 2 & g(:,5,i) < 2 & g(:,6,i) < 2);
   
   for j = 1:3
      figure
      cc(i,j) = corr(g(idxx,j + 3,i),g(idxx,7,i));
      scatter(log(g(idxx,j + 3,i)),log(g(idxx,7,i)))
      xlabel(lab{j + 3});      
      ylabel(lab{7});
      legend(['Code ' num2str(ucode(i))])
   end
end
end



for m = 1:1
   for j = 1:3
        datamat{m,j}= zeros(gsz,2); 
        maddatamat{m,j}= zeros(gsz,2); 
        for i = 1:gsz
            kt = 1;
             for k = 1:2
            % for k = 1:length(ufile{i})    
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
               % ttt = [ttt; kt*ones(length(idxx),1)];
                
                datamat{m,j}(i,k) =  median(sss);
                maddatamat{m,j}(i,k) = mad(sss);
                bigdatamat{m,j,i,k} = sss;   % m is measurement type, j is neurotransmitter, i is subject
                bigdatamat_dishcodes{m,j,i,k} = ttt;
                kt = kt + 1;
            end
        end
        if plotfigs 
            figure
         %   boxplot(sss,ttt,'labels',ucode);
            boxplot(sss,ttt,'labels',lll);
            title([lab{j + 3} ' KCL normalized'])
               
            
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
           for i = 1:gsz
               for k = 1:2
                   plot(between_factors(i) + 0.5*(k - 1) + 0.1*nsg(i),bigdatamat{m,j,i,k},col{i})
               end
           end
        end
    end
end



if do_nested
    % make y for anovan  
    for m = 1:1
        for j = 1:3
            y = [];
            subj = [];
            treat = [];
            dx = [];
            dish = [];
            clear tsubj
            for i = 1:gsz
                for k = 1:2
                    y = [y; bigdatamat{m,j,i,k}];  % m is measurement type, j is neurotransmitter, i is subject
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
			idx = find(subj < 5); % SZs and twins
            %***************************
            % % NESTED T-TESTS HERE
            % ***********************
			if validate 
                zp{m,j} = 0;
                zzz{m,j} = 0;
                iter = 10000;
                for aaa = 1:iter
			         yyy = 3*randn(size(y)) + 1; 
			         [zzzz,zpp] = nested_ttest_twins(yyy(idx),subj(idx),treat(idx),{[1 3],[2 4]});
                     zp{m,j} = zp{m,j} + double(zpp < 0.05)/iter;
                     zzz{m,j} = zzz{m,j} + double(zpp < 0.01)/iter;
                end
		    else
			    [zzz{m,j},zp{m,j},tstats{m,j}] = nested_ttest_twins_nochange(y(idx),subj(idx),treat(idx),{[1 3],[2 4]},fullfile(figdir,strcat('nested_t_',stattype{m},lab{j + 3})),[lab{j + 3} gtype2{ggg}]);
			end
           % [ppp{m,j}, ntbl{m,j}, stats{m,j}] = anovan(y,{treat,subj,dx},'nested',[0 0 0 ;0 0 1 ;0 0 0],'random',[2],'model','interaction','varnames',{'Treat','Subj','Dx'},'display','off');
           % cap = getString(message('stats:anovan:ConstrainedTypeIIISumsOfSquares'));
            digits = [(-1) (-1) (-1) (-1) (-1) 2 4 ];
            % figh = statdisptable_jt(ntbl{m,j}, getString(message('stats:anovan:NWayANOVA')), ['AnalysisOfVariance ' strcat(stattype{m},lab{j + 3})],cap,digits);
            %  text(0,0,strcat(stattype{m},lab{j + 3})) 
            % export_fig(fullfile(figdir,strcat('cellANOVA',stattype{m},lab{j + 3},'.png')) ,'-r300','-painters'); 
			% print('-dpng','-r300','-painters','-noui',fullfile(figdir,strcat('cellANOVA',stattype{m},lab{j + 3},'.png'))); 
            % f = figure('visible','off');
            % boxplot(y,{tsubj,treat})
            % title(strcat(stattype{m},lab{j + 3})) 
           % export_fig(fullfile(figdir,strcat('cellBox',stattype{m},lab{j + 3},'.png')) ,'-r300','-painters'); 
			% print('-dpng','-r300','-painters','-noui',fullfile(figdir,strcat('cellBox',stattype{m},lab{j + 3},'.png')));
            %   [ppp{m,j}, ntbl{m,j}, stats{m,j}] = anovan(y,{treat,subj,dx,dish},'nested',[0 0 0 0;0 0 1 0;0 0 0 0;1 1 0 0],'random',[2 4],'model','interaction','varnames',{'Treat','Subj','Dx','Dish'});
             
             
        end
    end
end
if ~validate
    save(fullfile(figdir,'nested_t_stats'),'zzz','tstats')
end

