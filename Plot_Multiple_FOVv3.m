% Plot Multiple FISH FOV
% Copyright (C) 2017 David McSwiggen

%%%%%%%%%%%%%%%%%%%%%%%% GNU LICENSE OVERVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, please see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; close all;



% All_Replicates{1} = {'Quant_FISH_rep_1_pos_1','Quant_FISH_rep_1_pos_2','Quant_FISH_rep_1_pos_3','Quant_FISH_rep_1_pos_4',...
%     'Quant_FISH_rep_1_pos_5','Quant_FISH_rep_1_pos_6','Quant_FISH_rep_1_pos_7'};%,'Quant_FISH_rep_1_pos_8','Quant_FISH_rep_1_pos_9',...
% % 'Quant_FISH_rep_1_pos_10','Quant_FISH_rep_1_pos_11','Quant_FISH_rep_1_pos_12','Quant_FISH_rep_1_pos_13','Quant_FISH_rep_1_pos_14',...
% %     'Quant_FISH_rep_1_pos_15'};%
All_Replicates{1} = {'Quant_FISH_rep_2_pos_1','Quant_FISH_rep_2_pos_2','Quant_FISH_rep_2_pos_3',...
    'Quant_FISH_rep_2_pos_5','Quant_FISH_rep_2_pos_6','Quant_FISH_rep_2_pos_7','Quant_FISH_rep_2_pos_8',...
    'Quant_FISH_rep_2_pos_9','Quant_FISH_rep_2_pos_10','Quant_FISH_rep_2_pos_13','Quant_FISH_rep_2_pos_14',...
    'Quant_FISH_rep_2_pos_15','Quant_FISH_rep_2_pos_16'};
All_Replicates{2} = {'Quant_FISH_rep_3_pos_1','Quant_FISH_rep_3_pos_2','Quant_FISH_rep_3_pos_3','Quant_FISH_rep_3_pos_4',...
     'Quant_FISH_rep_3_pos_5','Quant_FISH_rep_3_pos_6','Quant_FISH_rep_3_pos_7',...
     'Quant_FISH_rep_3_pos_9','Quant_FISH_rep_3_pos_10','Quant_FISH_rep_3_pos_11','Quant_FISH_rep_3_pos_12','Quant_FISH_rep_3_pos_13',...
     'Quant_FISH_rep_3_pos_14','Quant_FISH_rep_3_pos_15','Quant_FISH_rep_3_pos_16','Quant_FISH_rep_3_pos_17'};


Directory_input_path = 'PUT INPUT FILEPATH HERE';
Directory_output_path = 'PUT OUTPUT FILEPATH HERE';
Data_set_name = 'CHOOSE A NAME FOR THE DATA SET HERE';

Save_data = 1;

n_bootstrap = 500;
FNR_TRITC = 0.122;

PlotTitle = {'mRNA FISH, alternating probes'};
All_paired_spots_total = [];
All_TRITC_solo_total = [];
All_Cy5_solo_total = [];
All_counts_total = [];
FISH_data_condition_all = struct;

for replicate = 1:size(All_Replicates,2)
    disp(['Working on replicate number ', num2str(replicate)]);
    Replicate_paired_spots_total{replicate} = [];
    Replicate_TRITC_solo_total{replicate} = [];
    Replicate_Cy5_solo_total{replicate} = [];
    Replicate_counts_total{replicate} = [];
    
    %Load each cell of array with the aggregated data from one replicate
    for file=1:size(All_Replicates{replicate},2)
        load([Directory_input_path, cell2mat(All_Replicates{replicate}(file))]);
        current_struct = FISH_data_individual_cells;
        current_paired_spot_matrix = All_paired_spot_matrix;
        current_solo_TRITC_spot_matrix = All_solo_TRITC_spot_matrix;
        current_solo_Cy5_spot_matrix = All_solo_Cy5_spot_matrix;
        current_counting_stats_matrix = All_counting_stats;
        if file == 1
            FISH_data_condition_replicate = current_struct;
        else
            FISH_data_condition_replicate = horzcat(FISH_data_condition_replicate, current_struct);
        end
        Replicate_paired_spots_total{replicate} = vertcat(Replicate_paired_spots_total{replicate}, current_paired_spot_matrix);
        Replicate_TRITC_solo_total{replicate} = vertcat(Replicate_TRITC_solo_total{replicate}, current_solo_TRITC_spot_matrix);
        Replicate_Cy5_solo_total{replicate} = vertcat(Replicate_Cy5_solo_total{replicate}, current_solo_Cy5_spot_matrix);
        Replicate_counts_total{replicate} = vertcat(Replicate_counts_total{replicate}, current_counting_stats_matrix);
    end
    
    Replicate_counts_total{replicate}(:,9) = Replicate_counts_total{replicate}(:,3) ./ Replicate_counts_total{replicate}(:,5);
    Replicate_counts_total{replicate}(:,8) = 1-(Replicate_counts_total{replicate}(:,7) + Replicate_counts_total{replicate}(:,9));
    
    
     
    All_paired_spots_total = vertcat(All_paired_spots_total, Replicate_paired_spots_total{replicate});
    All_TRITC_solo_total = vertcat(All_TRITC_solo_total,Replicate_TRITC_solo_total{replicate});
    All_Cy5_solo_total = vertcat(All_Cy5_solo_total,Replicate_Cy5_solo_total{replicate});
    All_counts_total = vertcat(All_counts_total,Replicate_counts_total{replicate});
    
    if replicate == 1
        FISH_data_condition_all = FISH_data_condition_replicate;
    else
        FISH_data_condition_all = horzcat(FISH_data_condition_all,FISH_data_condition_replicate);
    end
    
    %TEST UNPAIRED SPOTS FOR POISSON-NESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_unpaired_Cy5 = All_counts_total(:,3);
    n_unpaired_TRITC = All_counts_total(:,2);
    n_unpaired_total = n_unpaired_Cy5 + n_unpaired_TRITC;
    % Generate frequency lists for each spot type, all with the same
    % binning
    edges = [0:max(n_unpaired_total)+1];
    [hist_unpaired_total, edges] = histcounts(n_unpaired_total,edges);
    poiss_total = fitdist(edges(1:end-1)','Poisson','Frequency',hist_unpaired_total');
    expected_total = size(n_unpaired_total,1)*pdf(poiss_total,edges(1:end-1));
    [hyp,pval,statistics] = chi2gof(edges(1:end-1),'Ctrs',edges(1:end-1),'Frequency',hist_unpaired_total,...
        'Expected',expected_total,'NParams',1);
    if pval < 0.05
        disp('Good News! The total number of unpaired spots in this condition have a Poisson distribution!');
    else
        disp('Bad News! The total number of unpaired spots do not have a Poisson distribution.');
    end
    
    hist_unpaired_TRITC = histcounts(n_unpaired_TRITC,edges);
    poiss_TRITC = fitdist(edges(1:end-1)','Poisson','Frequency',hist_unpaired_TRITC');
    expected_TRITC = size(n_unpaired_TRITC,1)*pdf(poiss_TRITC,edges(1:end-1));
    [hyp,pval,statistics] = chi2gof(edges(1:end-1),'Ctrs',edges(1:end-1),'Frequency',hist_unpaired_TRITC,...
        'Expected',expected_TRITC,'NParams',1)
    if pval < 0.05
        disp('Good News! The unpaired TRITC spots in this condition have a Poisson distribution!');
    else
        disp('Bad News! The unpaired TRITC spots do not have a Poisson distribution.');
    end
    
    
    [hist_unpaired_Cy5] = histcounts(n_unpaired_Cy5,edges);
    poiss_Cy5 = fitdist(edges(1:end-1)','Poisson','Frequency',hist_unpaired_Cy5');
    expected_Cy5 = size(n_unpaired_Cy5,1)*pdf(poiss_Cy5,edges(1:end-1));
    [hyp,pval,statistics] = chi2gof(edges(1:end-1),'Ctrs',edges(1:end-1),'Frequency',hist_unpaired_Cy5,...
        'Expected',expected_Cy5,'NParams',1);
    if pval < 0.05
        disp('Good News! The unpaired Cy5 spots in this condition have a Poisson distribution!');
    else
        disp('Bad News! The unpaired Cy5 spots do not have a Poisson distribution.');
    end


    %BOOTSTRAP RESAMPLE EACH REPLICATE TO GENERATE MEANS AND STDEV%%%%%%%%%
    data_set_size = size(Replicate_counts_total{replicate},1);
    tic;
    bootstrapped_counts = NaN(floor(data_set_size/2),18);
    
    for number_to_choose = 1:floor(data_set_size/2)
        for resample = 2:n_bootstrap
            random_resampled_data = datasample(Replicate_counts_total{replicate},number_to_choose,1,'Replace',false);
            %random_data = datasample(All_counting_stats,resample,1,'Replace',false);
            sampling(resample,:) = mean(random_resampled_data,1,'omitnan');
        end
        
        bootstrapped_counts(number_to_choose,1:2:17) = mean(sampling,'omitnan');
        bootstrapped_counts(number_to_choose,2:2:18) = 2*(std(sampling,'omitnan'));
    end
    
    bootstrapped_replicates(replicate,:) = bootstrapped_counts(end,:);
    toc;
    
    
    if Save_data == 1
        save([Directory_output_path,Data_set_name,'_rep_',num2str(replicate)],'Replicate_paired_spots_total','Replicate_TRITC_solo_total',...
            'Replicate_Cy5_solo_total','Replicate_counts_total','FISH_data_condition_replicate','bootstrapped_counts');
    end
    
end

% Generate a matrix to use for boxpolts comparing different replicates

Boxplot_replicates_matrix = NaN(size(All_counts_total,1),(size(All_Replicates,2))*4);
Boxplot_labels = cell(1,size(All_Replicates,2)+1);
for rep = 1:size(All_Replicates,2)
   Boxplot_replicates_matrix(1:size(Replicate_counts_total{rep},1),1+((rep-1)*4)) = Replicate_counts_total{rep}(:,1);
   Boxplot_replicates_matrix(1:size(Replicate_counts_total{rep},1),2+((rep-1)*4)) = Replicate_counts_total{rep}(:,2);
   Boxplot_replicates_matrix(1:size(Replicate_counts_total{rep},1),3+((rep-1)*4)) = Replicate_counts_total{rep}(:,3);
   Boxplot_replicates_matrix(1:size(Replicate_counts_total{rep},1),4+((rep-1)*4)) = Replicate_counts_total{rep}(:,6);
   Boxplot_labels{rep} = ['Rep ',num2str(rep)];
end

Boxplot_labels{end} = ['Pool'];


bootstrapped_all_mean = mean(bootstrapped_replicates,1, 'omitnan');
bootstrapped_all_stdevs = (sum(bootstrapped_replicates.^2,1)/size(All_Replicates,2)).^(1/2);
% %Column 1: mean paired spots
% %Column 2: Stdev paired spots
% %Column 3: mean TRITC only spots
% %Column 4: Stdev TRITC only spots
% %Column 5: mean Cy5 only spots (ORF)
% %Column 6: Stdev Cy5 only spots (ORF)
% %Column 7: mean all TRITC spots
% %Column 8: Stdev all TRITC spots
% %Column 9: mean Cy5 total
% %Column 10: Stdev Cy5 total
% %Column 11: mean total
% %Column 12: Stdev total
% %Column 13: mean fraction paired (LUTI)
% %Column 14: Stdev fraction paired (LUTI)
%bootstrap_all_stats(:,1) = ['LUTI mRNA','TRITC only','ORF mRNA','all TRITC','all Cy5','all spots','Fraction LUTI']
bootstrap_all_stats = zeros(9,2);
bootstrap_all_stats(:,1) = bootstrapped_all_mean(1:2:17)';
bootstrap_all_stats(:,2) = bootstrapped_all_stdevs(2:2:18)';



All_paired_spot_matrix = All_paired_spots_total;
All_solo_TRITC_spot_matrix = All_TRITC_solo_total;
All_solo_Cy5_spot_matrix = All_Cy5_solo_total;
All_counting_stats_raw = All_counts_total;
% All_counting_stats_raw(:,8) = All_counting_stats_raw(:,2)./All_counting_stats_raw(:,4);
% All_counting_stats_raw(:,9) = All_counting_stats_raw(:,3)./All_counting_stats_raw(:,5);
%Column 1: Number paired (LUTI)
%Column 2: Number TRITC alone (False neg)
%Column 3: Number Cy5 alone (ORF)
%Column 4: Total TRITC
%Column 5: Total Cy5
%Column 6: Total detections
%Column 7: Fraction LUTI
%Column 8: Fraction of unpaired TRITC
%Column 9: Fraction of unpaired Cy5

n_cells = size(All_counts_total,1);
Boxplot_all_matrix = horzcat(Boxplot_replicates_matrix,[All_counting_stats_raw(:,1:3),All_counting_stats_raw(:,6)]);

if Save_data == 1
        save([Directory_output_path,Data_set_name,'_all_reps'],'All_paired_spot_matrix','All_solo_TRITC_spot_matrix','All_solo_Cy5_spot_matrix',...
            'All_counting_stats_raw','FISH_data_condition_replicate','Replicate_paired_spots_total','Replicate_TRITC_solo_total',...
            'Replicate_Cy5_solo_total','Replicate_counts_total','FISH_data_condition_replicate','bootstrap_all_stats');
    end


%% Determine Stats as a function of expression level

max_total_mRNA = max(All_counting_stats_raw(:,6));
stats_per_expression_level = NaN(max_total_mRNA, 18);
for expression_level = 1:max_total_mRNA
    cell_index = find(All_counting_stats_raw(:,6) == expression_level);
    if size(cell_index,1) > 1
        stats_per_expression_level(expression_level,1:2:17) = mean(All_counting_stats_raw(cell_index,:),'omitnan');
        stats_per_expression_level(expression_level,2:2:18) = std(All_counting_stats_raw(cell_index,:),1,'omitnan');
    elseif size(cell_index,1) == 1
        stats_per_expression_level(expression_level,1:2:17) = All_counting_stats_raw(cell_index,:);
        stats_per_expression_level(expression_level,2:2:18) = 0;
    end
    
end
Percent_per_expression = [stats_per_expression_level(:,15),stats_per_expression_level(:,17),stats_per_expression_level(:,13)];

% %% Bootstrap Resampling of histogram bins
% bin_edges = [1:2:max(All_counting_stats_raw(:,6))];
% Histogram_bins_confidence.Cy5_only_counts = NaN(n_bootstrap, size(bin_edges,2));
% Histogram_bins_confidence.TRITC_only_counts = NaN(n_bootstrap, size(bin_edges,2));
% Histogram_bins_confidence.paired_counts = NaN(n_bootstrap, size(bin_edges,2));
% 
% for sample = 1:n_bootstrap
%     resampled_pooled_data = datasample(All_counting_stats_raw,floor(size(All_counting_stats_raw,1)/2),1,'Replace',false);
%     Histogram_bins_confidence.Cy5_only_counts(sample,1:(end-1)) = histcounts(resampled_pooled_data(:,3),bin_edges,'Normalization','probability');
%     Histogram_bins_confidence.TRITC_only_counts(sample,1:(end-1)) = histcounts(resampled_pooled_data(:,2),bin_edges,'Normalization','probability');
%     Histogram_bins_confidence.paired_counts(sample,1:(end-1)) = histcounts(resampled_pooled_data(:,1),bin_edges,'Normalization','probability');
% 
% end
% 
% Histogram_bins_confidence.Cy5_bootstrapped_stats =...
%     [mean(Histogram_bins_confidence.Cy5_only_counts,'omitnan');2*std(Histogram_bins_confidence.Cy5_only_counts,1,'omitnan')];
% Histogram_bins_confidence.TRITC_bootstrapped_stats =...
%     [mean(Histogram_bins_confidence.TRITC_only_counts,'omitnan');2*std(Histogram_bins_confidence.TRITC_only_counts,1,'omitnan')];
% Histogram_bins_confidence.paired_bootstrapped_stats =...
%     [mean(Histogram_bins_confidence.paired_counts,'omitnan');2*std(Histogram_bins_confidence.paired_counts,1,'omitnan')];
% 
% errorbar(bin_edges(:) + 0.5, Histogram_bins_confidence.Cy5_bootstrapped_stats(1,:), Histogram_bins_confidence.Cy5_bootstrapped_stats(2,:),'r');
% hold on
% errorbar(bin_edges(:) + 0.5, Histogram_bins_confidence.paired_bootstrapped_stats(1,:), Histogram_bins_confidence.paired_bootstrapped_stats(2,:),'k');

%% Single cell heatmap analysis

matrix_size = max([max(All_counting_stats_raw(:,1)),max(All_counting_stats_raw(:,3))]);
single_cell_heatmap = zeros(matrix_size,matrix_size);
for row = 1:matrix_size
    for col = 1:matrix_size
        pixel_val = size(find(and(All_counting_stats_raw(:,3) == row,All_counting_stats_raw(:,1) == col)),1);
        %% 
        single_cell_heatmap(row,col) = pixel_val;
    end
end

imagesc(single_cell_heatmap);
axis image
        

%% BAYESIAN TESTING FOR CREDIBLE INTERVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This module uses Bayesian inference to predict, based on a false negative
% rate (FNR), the true number of Cy5 spots that should have been paired up
% with a TRITC spot. Presumably the rest of the spots are truly Cy5 only,
% and represent only the ORF form of the transcript. The FNR will be
% determined by the Odd/Even data set.
All_counting_stats_corrected = zeros(size(Replicate_counts_total,1),10);

for cell = 1:n_cells
    n_red_only_obs = All_counting_stats_raw(cell,3);
    n_grn_only_obs = All_counting_stats_raw(cell,2);
    n_pairs_obs = All_counting_stats_raw(cell,1);
    n_red_all_obs = All_counting_stats_raw(cell,5);
    n_grn_all_obs = All_counting_stats_raw(cell,4);
    n_cell_total = All_counting_stats_raw(cell,6);
    
    pair_red_prob_array = zeros(n_red_only_obs);
    if n_grn_all_obs > n_pairs_obs
        if n_red_only_obs > 0
            for j = 0:(n_red_only_obs)
                numerator = ((FNR_TRITC^j)*(1-FNR_TRITC)^(n_red_only_obs-j)*nchoosek(n_red_only_obs,j)); %probability (assume binomial) of missing j TRITC spots out of n red-only spots
                for i = 0:n_red_only_obs
                    denom(i+1) = ((FNR_TRITC^i)*(1-FNR_TRITC)^(n_red_only_obs-i)*nchoosek(n_red_only_obs,i)); %Sum of all probabilities for j missed TRITC spots
                end
                sum_denom = sum(denom(:));
                prob(n_pairs_obs + j +1) = numerator/sum_denom; %conditional probability for j paired spots, given n_pairs_obs observed pairs
            end
            
            probability_dist = prob(n_pairs_obs +1 :n_red_all_obs +1 ).*[1:n_red_only_obs+1];
            [lambdahat,CI] = poissfit(probability_dist,0.05);
            
            
            
            %Rounding to the nearest whole mRNA
            n_red_most_probable = n_red_only_obs - lambdahat; %most likely number of red spots
            if n_red_only_obs - CI(2) >= 0
                n_red_lower = n_red_only_obs - CI(2); %2.5th percentile for number of red-only spots
            else
                n_red_lower = 0;
            end
            if n_red_only_obs - CI(1) >= 0
                n_red_upper = n_red_only_obs - CI(1);%97.5th percentile for number of red-only spots
            else
                n_red_upper = 0;
            end
            n_pair_most_probable = n_red_all_obs - n_red_most_probable;
            n_pair_lower = n_red_all_obs - n_red_upper;
            n_pair_upper = n_red_all_obs - n_red_lower;
            
        else
            n_red_most_probable = n_red_only_obs;
            n_red_lower = n_red_only_obs;
            n_red_upper = n_red_only_obs;
            n_pair_most_probable = n_pairs_obs;
            n_pair_lower = n_pairs_obs;
            n_pair_upper = n_cell_total;
            
        end
    else
        n_red_most_probable = n_red_only_obs;
        n_red_lower = n_red_only_obs;
        n_red_upper = n_red_only_obs;
        n_pair_most_probable = n_pairs_obs;
        n_pair_lower = n_pairs_obs;
        n_pair_upper = n_cell_total;
        
    end
    
    
    All_counting_stats_corrected(cell,:) = [n_grn_only_obs, n_red_only_obs, n_red_lower, n_red_most_probable, n_red_upper,...
        n_pairs_obs, n_pair_lower, n_pair_most_probable, n_pair_upper, (n_pair_most_probable/n_cell_total)];
    %Column 1: Number of TRITC observed (false positive)
    %Column 2: Number of Cy5 only observations (ORF only)
    %Column 3: 2.5th percentile number of Cy5 alone
    %Column 4: Most probable number of ORF mRNA
    %Column 5: 97.5th percentile number of Cy5 alone
    %Column 6: Observed pairs (LUTI)
    %Column 7: 2.5th percentile number of LUTI RNA
    %Column 8: Most probable number of LUTI RNA
    %Column 9: 97.5th percentile number of LUTI RNA
    %Column 10: Most probable fraction of LUTI RNA
end

%% SAVE ALL THE DATA
if Save_data == 1
    save([Directory_output_path,Data_set_name,'_all_reps'],'All_paired_spot_matrix','All_solo_TRITC_spot_matrix','All_solo_Cy5_spot_matrix',...
        'All_counting_stats_raw','bootstrap_all_stats','All_counting_stats_corrected');
    
    dlmwrite([Directory_output_path,Data_set_name,'_raw_stats.txt'],All_counting_stats_raw,'delimiter','\t');
    dlmwrite([Directory_output_path,Data_set_name,'_corrected_stats.txt'],All_counting_stats_corrected,'delimiter','\t');
    dlmwrite([Directory_output_path,Data_set_name,'_bootstrap_averages.txt'],bootstrap_all_stats  ,'delimiter','\t');
    
end

%% PLOT "FINAL" SPOT PAIRS
% Compound_figure = figure('position',[20 100 2200 1900]);
% s1 = subplot(8,8,1:3);
% % Plot the pairwise distances of "true pairs" and other spots
% spot_distances_pairs = histogram(All_paired_spot_matrix(:,7));
% spot_distances_pairs.Normalization = 'pdf';
% spot_distances_pairs.BinWidth = 0.5;
% spot_distances_pairs.FaceColor = [0.8185,0.7327,0.3498];
% hold on;
% spot_distances_TRITC = histogram(All_solo_TRITC_spot_matrix(:,7));
% spot_distances_TRITC.Normalization = 'pdf';
% spot_distances_TRITC.BinWidth = 0.5;
% spot_distances_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
% spot_distances_Cy5 = histogram(All_solo_Cy5_spot_matrix(:,7));
% spot_distances_Cy5.Normalization = 'pdf';
% spot_distances_Cy5.BinWidth = 0.5;
% spot_distances_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% axis([0 20 0 1]);
% title ('Distance from Cy5 spot to closest neighbor');
% title ('Distance to nearest neighbor','FontSize',10, 'FontName', 'Helvetica');
% xlabel('Distance (px)', 'FontSize',8, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',8, 'FontName', 'Helvetica');
% legend('Paired detections','Unpaired TRITC detections','Unpaired Cy5 detections','Box','off','FontSize',14,'Location','North');

Compound_figure = figure('position',[20 100 2200 1900]);
s1 = subplot(8,8,1:3);
% Plot the pairwise distances of "true pairs" and other spots
spot_distances_pairs_CDF = cdfplot(All_paired_spot_matrix(:,7));
spot_distances_pairs_CDF.Color = [0.8185,0.7327,0.3498];
spot_distances_pairs_CDF.LineWidth = 2;
hold on;
spot_distances_TRITC_CDF = cdfplot(All_solo_TRITC_spot_matrix(:,7));
spot_distances_TRITC_CDF.Color = [0.1801, 0.7177, 0.6424];
spot_distances_TRITC_CDF.LineWidth = 2;
spot_distances_Cy5_CDF = cdfplot(All_solo_Cy5_spot_matrix(:,7));
spot_distances_Cy5_CDF.Color = [0.208, 0.1663, 0.5292];
spot_distances_Cy5_CDF.LineWidth = 2;
hold off;
axis([0 20 0 1]);
title ('Distance from Cy5 spot to closest neighbor');
title ('Distance to nearest neighbor','FontSize',10, 'FontName', 'Helvetica');
xlabel('Distance (px)', 'FontSize',8, 'FontName', 'Helvetica');
ylabel('Probability', 'FontSize',8, 'FontName', 'Helvetica');
legend('Paired detections','Unpaired TRITC detections','Unpaired Cy5 detections','Box','off','FontSize',14,'Location','North');

s2 = subplot(8,8,[5:8,13:16]);
% Plot the spot intensities paired and un-paired spots with TRITC
spot_intensity_pairs = histogram(All_paired_spot_matrix(:,10),50);
normalized_spot_intensity_pairs_tritc = spot_intensity_pairs.Values ./sum(spot_intensity_pairs.Values);
hold on;
spot_intensity_pairs.Normalization = 'pdf';
spot_intensity_pairs.BinWidth = 250;
spot_intensity_pairs.FaceColor = [0.8185,0.7327,0.3498];
spot_intensity_TRITC = histogram(All_solo_TRITC_spot_matrix(:,3),50);
normalized_spot_intensity_solo_tritc = spot_intensity_TRITC.Values ./sum(spot_intensity_TRITC.Values);
spot_intensity_TRITC.Normalization = 'pdf';
spot_intensity_TRITC.BinWidth = 250;
spot_intensity_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
hold off;
title ('Integrated spot intensity for paired and unpaired','FontSize',10, 'FontName', 'Helvetica');
ylabel('Probability', 'FontSize',8, 'FontName', 'Helvetica');
legend('Paired detections','Unpaired detections','Box','off','Location','North');

s3 = subplot(8,8,[21:24,29:32]);
% Plot the spot intensities paired and un-paired spots with Cy5
spot_intensity_pairs = histogram(All_paired_spot_matrix(:,3),50);
normalized_spot_intensity_pairs_cy5 = spot_intensity_pairs.Values ./sum(spot_intensity_pairs.Values);
hold on;
spot_intensity_pairs.Normalization = 'pdf';
spot_intensity_pairs.BinWidth = 250;
spot_intensity_pairs.FaceColor = [0.8185,0.7327,0.3498];
spot_intensity_Cy5 = histogram(All_solo_Cy5_spot_matrix(:,3),50);
normalized_spot_intensity_solo_cy5 = spot_intensity_Cy5.Values ./sum(spot_intensity_Cy5.Values);
spot_intensity_Cy5.Normalization = 'pdf';
spot_intensity_Cy5.BinWidth = 250;
spot_intensity_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
hold off;
xlabel('Integrated Spot Intensity (A.U.)', 'FontSize',10, 'FontName', 'Helvetica');
ylabel('Probability', 'FontSize',8, 'FontName', 'Helvetica');
legend('Paired detections','Unpaired detections','Box','off','Location','North');

linkaxes([s2, s3], 'x');
set(s2, 'xticklabel',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the four boxplot sets
s4 = subplot(8,8,[41,49,57]);
Box_plot = boxplot(Boxplot_all_matrix(:,1:4:end),'notch','on','labels',Boxplot_labels,'colors','k','outliersize',4,'symbol','ko');
title('Paired','FontSize',10, 'FontName', 'Helvetica');
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.8185,0.7327,0.3498],'FaceAlpha',.5);
 end

s5 = subplot(8,8,[42,50,58]);
Box_plot = boxplot(Boxplot_all_matrix(:,2:4:end),'notch','on','labels',Boxplot_labels,'colors','k','outliersize',4,'symbol','ko');
title('TRITC only','FontSize',10, 'FontName', 'Helvetica');
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.1801, 0.7177, 0.6424],'FaceAlpha',.5);
 end
 
 s6 = subplot(8,8,[43,51,59]);
Box_plot = boxplot(Boxplot_all_matrix(:,3:4:end),'notch','on','labels',Boxplot_labels,'colors','k','outliersize',4,'symbol','ko');
title('Cy5 only','FontSize',10, 'FontName', 'Helvetica');
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.208, 0.1663, 0.5292],'FaceAlpha',.5);
 end
 
  s7 = subplot(8,8,[44,52,60]);
Box_plot = boxplot(Boxplot_all_matrix(:,4:4:end),'notch','on','labels',Boxplot_labels,'colors','k','outliersize',4,'symbol','ko');
title('All spots','FontSize',10, 'FontName', 'Helvetica');
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.6, 0.6, 0.6],'FaceAlpha',.5);
 end
 
%Plot the bootstrapping results
s8 = subplot(8,8,[46:48,54:56,62:64]);
x_vector = 1:size(bootstrapped_counts,1);
hold on;
lineProps.col{1} = [0.8185,0.7327,0.3498];
lineProps.width = 3;
mseb(x_vector,bootstrapped_counts(:,13)',bootstrapped_counts(:,14)',lineProps,1);
lineProps.col{1} = [0.1801, 0.7177, 0.6424];
lineProps.width = 3;
mseb(x_vector,bootstrapped_counts(:,15)',bootstrapped_counts(:,16)',lineProps,1);
lineProps.col{1} = [0.208, 0.1663, 0.5292];
lineProps.width = 3;
mseb(x_vector,bootstrapped_counts(:,17)',bootstrapped_counts(:,18)',lineProps,1);
title(['Bootstrapping results for ',num2str(n_bootstrap),' bootstrap iterations']);
xlabel('number of cells randomly sampled','FontSize',8, 'FontName', 'Helvetica');
ylabel('mean percent of mRNA');
axis([0 ((size(bootstrapped_counts,1))*1.05) 0 1]);
hold off

hold on;
s9 = subplot(8,8,[17:19,25:27]);
stacked_bar = bar(Percent_per_expression(:,1:3),'stacked','LineWidth',0.10);
%axis([0 1 0 min(All_counting_stats_raw(:,6)) max(All_counting_stats_raw(:,6))]);
title('Fraction of paired and unpaired spots, by total mRNA amount');
stacked_bar(1).EdgeColor = [0.1001, 0.6377, 0.5624];
stacked_bar(1).FaceColor = [0.1801, 0.7177, 0.6424];
stacked_bar(2).EdgeColor = [0.128, 0.1263, 0.4492];
stacked_bar(2).FaceColor = [0.208, 0.1663, 0.5292];
stacked_bar(3).EdgeColor = [0.855,0.87,0.6127];
stacked_bar(3).FaceColor = [0.935,0.95,0.6927];
legend('TRITC only','Cy5 only','Paired','Location','Northwest')
hold off;

img = getframe(gcf);
imwrite(img.cdata, [Directory_output_path,Data_set_name,'_comp_fig.png']);

% Plot number of mRNA detections per cell for
% subplot(1,2,1)
% spot_counts_total_LUTI_raw = histogram(All_counting_stats_raw(:,1));
% hold on;
% spot_counts_total_LUTI_raw.Normalization = 'pdf';
% spot_counts_total_LUTI_raw.BinWidth = 1;
% spot_counts_total_LUTI_raw.FaceColor = [0.1801, 0.7177, 0.6424];
% spot_counts_total_ORF_raw = histogram(All_counting_stats_raw(:,3));
% spot_counts_total_ORF_raw.Normalization = 'pdf';
% spot_counts_total_ORF_raw.BinWidth = 1;
% spot_counts_total_ORF_raw.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% title ('Number of mRNA per cell','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Total number of mRNA per cell', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('LUTI mRNA','ORF mRNA', 'Box','off','FontSize',14,'Location','North');


%% Other figs
% figure('position',[20 100 1600 1200])
% % Plot number of mRNA detections per cell for
% subplot(1,2,1)
% spot_counts_total_LUTI_raw = histogram(All_counting_stats_raw(:,7));
% spot_counts_total_LUTI_raw.Normalization = 'pdf';
% spot_counts_total_LUTI_raw.BinWidth = 0.05;
% spot_counts_total_LUTI_raw.FaceColor = [0.1801, 0.7177, 0.6424];
% 
% 
% subplot(1,2,2)
% spot_counts_total_LUTI_corrected = histogram(All_counting_stats_corrected(:,10));
% spot_counts_total_LUTI_corrected.Normalization = 'pdf';
% spot_counts_total_LUTI_corrected.BinWidth = 0.05;
% spot_counts_total_LUTI_corrected.FaceColor = [0.1801, 0.7177, 0.6424];
% title ('Number of mRNA per cell','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Total number of mRNA per cell after bayesian correction', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('LUTI mRNA','ORF mRNA', 'Box','off','FontSize',14,'Location','North');



% figure;
% [a,Vals] = rose(All_paired_spot_matrix(:,14), 20);
% Vals = Vals./sum(Vals);
% MaxLim = max([0.025 1.01*max(Vals)]);
% h = polar(0,MaxLim);
% delete(h);
% set(gca, 'Nextplot','add')
% %# draw patches instead of lines: polar(t,r)
% [x,y] = pol2cart(a,Vals);
% h = patch(reshape(x,4,[]), reshape(y,4,[]), [0.1801, 0.7177, 0.6424]);
% title(['Angular displacement between TRITC and Cy5 Channels'], 'FontSize', 10, 'FontName', 'Helvetica');
% xlabel('angle in degrees', 'FontSize',10, 'FontName', 'Helvetica');
% 
% figure
% for replicate = 1:size(All_Replicates,2)
%     violin(All_counting_stats_raw(:,1),'facecolor',[0.1801, 0.7177, 0.6424],...
%         'edgecolor','k','facealpha',0.3,'mc','k','medc',[]);
%     title(['Replicate number ',replicate,' ORF transcripts']);
% end
