% ***************************************************************
% QUANTIFICATION MAIN SCRIPT, MALES
% ***************************************************************
% Expression of functional receptors to the main neurotransmitters
% glutamate and GABA was tested in neuronal cultures by measuring
% changes in the level of Intracellular calcium during fast 
% application of respective agonists. To measure the intracellular
% Ca2+ transients cells were loaded with the cell-permeable 
% indicator Fluo-4AM (F10471, Life Technologies, Eugene, OR, USA).
% This fluorescent dye was applied for 30 min at 37蚓, followed by
% a 10 min washout with a basic solution containing (in mM): 
% 152 NaCl, 10 glucose, 10 HEPES, 2.5 KCl, 2 CaCl2, 1 MgCl2; 
% pH 7.4. The cells were then transferred to the recording chamber 
% continuously rinsed with the same basic solution. Cell fluorescent 
% signals were acquired by a microscopy (Olympus, IX70) on the Till
% Photonics imaging setup (FEI GmbH, Munich, Germany) equipped with 
% a 12-bit CCD Camera (SensiCam, Kelheim, Germany) with a light 
% excitation wavelength of 494nm and adequate filters. Glutamate 
% (100 然 with the co-agonist glycine 10 然) or GABA (100 然) were 
% applied for 2s by a fast perfusion system (RSC-200, BioLogic 
% Science Instruments, Grenoble, France). KCl (30mM) was then 
% applied for neuronal cell type identification. At the end cells 
% were challenged with ionomycin (10 然, 2s) for normalization of 
% tested signals.

% All samples (autumn and April ) with numbers 1 or 2 are non treated. 
% The samples with numbers 3 and 4 are treated with just one exception: 3 SZ19 from2018.11.3 is not treated. 


clear
savefiles = 1;
saveaveragetc = 0;
% addpath C:\Users\justoh\matlab\FluoroSNNAP15.04.18\FluoroSNNAP15.04.18\FluoroSNNAP_code\oopsi-master
if isunix
    datadir{1} = '/research/groups/tohkagroup/rashid/data/schizo_males';
	resdir = '\research/groups/tohkagroup/neurotrans/schizo_males\';
else
	datadirbase = 'C:\Users\justoh\data\Rashid\schizo_males';
    datadir{1} = fullfile(datadirbase,'2017.10.30'); % protocol 1
    datadir{2} = fullfile(datadirbase,'2017.11.3');  % protocol 1
    datadir{3} = fullfile(datadirbase,'2018.4.10');  % protocol 2
    datadir{4} = fullfile(datadirbase,'2018.4.11');  % protocol 2
    datadir{5} = fullfile(datadirbase,'2017.4.17');  % protocol 2
    datadir{6} = fullfile(datadirbase,'2017.4.18');  % protocol 2
    
	resdir = 'C:\Users\justoh\results\neurotrans\sz_males\';
end
	


resfn_sti = 'SZ_cell_stimulus_amplitudes_males200818';
resfn_summaries = 'SZ_stimulus_amplitudes_summaries_males200818';
resfn_tc = 'SZ_average_tc';

xlsfn_sti = fullfile(resdir,resfn_sti);
xlsfn_sti_sum = fullfile(resdir,resfn_summaries);
xlsfn_tc = fullfile(resdir,resfn_tc);

iter = 0;
for j = 1:length(datadir)
    d = dir(datadir{j});
    for i = 1:length(d)
      if  length(d(i).name) > 5
          if isempty(strfind(d(i).name,'background')) & strcmp(d(i).name((end - 3):end),'xlsx')
             iter = iter + 1;
             if j < 3
                 protocol(iter) = 1;
             else
                 protocol(iter) = 2;
             end
             fn{iter} = d(i).name;
             [num{iter},txt{iter},raw{iter}] = xlsread(fullfile(datadir{j},fn{iter}),'Sheet1'); % chanelles num(:,2:2:)
             fnbg{iter} = [fn{iter}(1:(end - 5)) '_background' fn{iter}((end -4):end)];
             [numbg{iter},txtbg{iter},rawbg{iter}] = xlsread(fullfile(datadir{j},fnbg{iter}),'Sheet1');

          end          
      end
    end
end
% The protocol for the earlier (2017) days: 
% BS no Mg 60 sec  60 
% Glutamate no Mg 100uM 2 sec 62
% BS no Mg 60 sec 122
% BS 120 sec 242
% Glutamate 100uM 2 sec 244
% BS 180 sec 424
% GABA 100uM 2 sec 426
% BS 180 sec 606 
% KCL 30mM 2 sec 608 
% BS 120 sec 728 
% ionomycin 2 sec 730 
% BS 60 sec 790

% The protocol for the latter (2018) days 
%  5 min BS        300
% 2 sec KCL 30mM   302
% 3 min BS         482
% 2 sec GABA 100uM 484 
% 2 min BS         604
% 1min BS no Mg    664
% 2 sec glutamate 100uM & glycine 10uM 666
% 1 min BS no Mg 726
% 2 min BS 846
% 2 sec ionomycin 848
% 1 min BS 908

% THIS IS THE OLD PROTOCOL: 
% *********************************************************
% 1 min BS  60 121
% 2 sec glutamate 100uM & glycine 10uM no Mg 62 125
% 1 min BS no Mg 122 245
% 2 min BS 242
% 2 sec glutamate 100uM & glycine 10uM 244
% 3 min BS 424
% 2 sec GABA 100uM 426
% 3 min BS 606
% 2 sec KCL 608
% 2 min BS 728
% 2 sec ionomycin 730
% 1 min BS 790
% ********************************************************
times{1} = [0 60 62 122 242 244 424 426 606 608 728 730 790];
times{2} = [0 300 302 482 484 604 664 666 726 846 848 908];
a_sptimesB{1} = [0 62   244 426 608 730]; % sp = baseline
a_sptimesE{1} = [60 242 424 606 728 790];
a_sptimesB{2} = [0 302 484 666 848];
a_sptimesE{2} = [300 482 664 846 908];

delay = 40;
onset = 10;

for j = 1:iter
    tp = num{j}(:,1);
    sptimesB  = a_sptimesB{protocol(j)};
    sptimesE  = a_sptimesE{protocol(j)};
    tpspB = zeros(size(sptimesB));
    tpspE = zeros(size(sptimesE));
    for i = 1:length(sptimesB)  
        [~,tpspB(i)] = min(abs(tp - sptimesB(i)));
        [~,tpspE(i)] = min(abs(tp - sptimesE(i)));        
    end
   sptimepoints = [tpspB(1):(tpspE(1) - onset)];
   for i = 1:length(sptimesB)
       sptimepoints = [sptimepoints (tpspB(i) + delay):(tpspE(i) - onset)];
   end           
   sz = size(num{j});
   ioind = (tpspE(end -1) - onset):(tpspB(end) + delay);
   [dts{j},spaves{j},udts{j},msts{j},mstsf0{j},uts{j}] = neurotrans_analyze_timeseries_simple(num{j},numbg{j},[2:2:sz(2)],[2:2:6],sptimepoints,1:sz(1),ioind);
% s3{j} = s3{j}/sbase;
end

% neurotrans_write_exel(xlsfn,stats,fn);
% neurotrans_plot_heatmaps(fn,resdir,1:30,dts);
% tois in order Glu no mg, Glu, GABA, KCl ionomycin
% Note that the only Glu  measurement in the protocol 2 is "no mg" measurement; 
% with Mg missing from protocol 2 
% Dummy Glu no Mg added instead 
toi{1} = [60 242 424 606 728];
toi{2} = [664 664 482 300 846];
for j = 1:iter
   [ amplitudes_gnorm{j}, amplitudes_spnorm{j}, amplitudes_lnorm{j}] = neurotrans_stimulus_amplitude(dts{j},num{j}(:,1),toi{protocol(j)},udts{j});
end
[code,index,dateCode,indexDate] = neurotrans_fn2code_sz_males(fn);
treated{1} = [1 2]; % untreated cells
treated{2} = [3 4]; % treated cells
if saveaveragetc
    for j = 1:2
        idx = find(protocol == j);
        [ucode1,N1] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_denoised_normalized_ion_wmean',num2str(j),'.xls'),dts(idx),code(idx),index(idx),num{idx(1)}(:,1),treated);
        [ucode2,N2] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_normalized_ion_wmean',num2str(j),'.xls'),msts(idx),code(idx),index(idx),num{idx(1)}(:,1),treated);
        [ucode3,N3] =neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_normalized_f0_wmean',num2str(j),'.xls'),mstsf0(idx),code(idx),index(idx),num{idx(1)}(:,1),treated);
        [ucode4,N4] =neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_wmean',num2str(j),'.xls'),uts(idx),code(idx),index(idx),num{idx(1)}(:,1),treated);
    end
end
ucode = unique(code); % this gives the ordering of subjects in the g -matrix
if savefiles
    g_gnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_gnorm.xls'),amplitudes_gnorm,fn,code,indexDate,index,dateCode);
    g_lnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_lnorm.xls'),amplitudes_lnorm,fn,code,indexDate,index,dateCode);
    g_spnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_spnorm.xls'),amplitudes_spnorm,fn,code,indexDate,index,dateCode);
    save(strcat(xlsfn_sti,'_g_allnorm'),'g_gnorm','dts','spaves','udts','g_lnorm','g_spnorm','ucode') 
end
