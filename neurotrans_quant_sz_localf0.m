% QUANTIFICATION MAIN SCRIPT, FEMALES
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
clear
savefiles = 1; % save files?
saveaveragetc = 0; % save average time-courses

if isunix
    datadir{1} = '/research/groups/tohkagroup/rashiddata/schizo';
	resdir = '\research/groups/tohkagroup/neurotrans/sz';
else
	datadir{1} = 'C:\Users\justoh\data\rashid\schizo';
	resdir = 'C:\Users\justoh\results\neurotrans\sz';
end
	

% resfn = 'results280218.xls';
resfn_sti = 'SZ_cell_stimulus_amplitudes240518';
resfn_summaries = 'SZ_stimulus_amplitudes_summaries240518';
resfn_tc = 'SZ_average_tc';
% xlsfn = fullfile(resdir,resfn);
xlsfn_sti = fullfile(resdir,resfn_sti);
xlsfn_sti_sum = fullfile(resdir,resfn_summaries);
xlsfn_tc = fullfile(resdir,resfn_tc);
%fn = 'spont activity network 17-07-2017.xlsx';
iter = 0;
for j = 1:length(datadir)
    d = dir(datadir{j});
    for i = 1:length(d)
      if  length(d(i).name) > 5
          if ~isempty(strfind(d(i).name,'2017')) & strcmp(d(i).name((end - 3):end),'xlsx')
             iter = iter + 1;
             fn{iter} = d(i).name;
             [num{iter},txt{iter},raw{iter}] = xlsread(fullfile(datadir{j},fn{iter}),'Sheet1'); % chanelles num(:,2:2:)
             fnbg{iter} = [fn{iter}(1:(end - 9)) 'background' fn{iter}((end -4):end)];
             [numbg{iter},txtbg{iter},rawbg{iter}] = xlsread(fullfile(datadir{j},fnbg{iter}),'Sheet1');

          end          
      end
    end
end
% The protocol is the same for all dates and dishes: 
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
times = [0 60 62 122 242 244 424 426 606 608 728 730 790];
sptimesB = [0 62   244 426 608 730];
sptimesE = [60 242 424 606 728 790];
tpspB = zeros(size(sptimesB));
tpspE = zeros(size(sptimesE));
delay = 40;
onset = 10;

for j = 1:iter
    tp = num{j}(:,1);
   
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

end


toi = [60 242 424 606 728];
for j = 1:iter
   [ amplitudes_gnorm{j}, amplitudes_spnorm{j}, amplitudes_lnorm{j}] = neurotrans_stimulus_amplitude(dts{j},num{j}(:,1),toi,udts{j});
end
[code,index,dateCode,indexDate] = neurotrans_fn2code_sz(fn);
if saveaveragetc
    j = 1;
    idx = 1:length(dts);
    [ucode1,N1] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_denoised_normalized_ion_wmean',num2str(j),'.xls'),dts(idx),code(idx),index(idx),num{idx(1)}(:,1));
    [ucode2,N2] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_normalized_ion_wmean',num2str(j),'.xls'),msts(idx),code(idx),index(idx),num{idx(1)}(:,1));
    [ucode3,N3] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_normalized_f0_wmean',num2str(j),'.xls'),mstsf0(idx),code(idx),index(idx),num{idx(1)}(:,1));
    [ucode4,N4] = neurotrans_saveaverage_tc10(strcat(xlsfn_tc,'_tc_wmean',num2str(j),'.xls'),uts(idx),code(idx),index(idx),num{idx(1)}(:,1));  
end

if savefiles 
    g_gnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_gnorm.xls'),amplitudes_gnorm,fn,code,indexDate,index,dateCode);
    g_lnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_lnorm.xls'),amplitudes_lnorm,fn,code,indexDate,index,dateCode);
    g_spnorm = neurotrans_write_excel_stimulus_cell_sz(strcat(xlsfn_sti,'_spnorm.xls'),amplitudes_spnorm,fn,code,indexDate,index,dateCode);
    save(strcat(xlsfn_sti,'_g_allnorm'),'g_gnorm','dts','spaves','udts','g_lnorm','g_spnorm') 
end

