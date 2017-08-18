%Random resample list comparison
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
clear all;clc


%%%% Bootstrapping %%%%
%%%randomly sample 1/2 of the Cy3 and Cy5 data sets. Compare mean, median,
%%%and overlap
bootstrap_sample_size = 450;
n_bootstrap = 1000;

Output_filepath = 'CHOOSE AN OUTPUT FILEPATH';
Case_name = '8144log_vs_8144_6h';
%Column 1: Number paired (LUTI)
%Column 2: Number TRITC alone (False neg)
%Column 3: Number Cy5 alone (ORF)
%Column 4: Total TRITC
%Column 5: Total Cy5
%Column 6: Total detections
%Column 7: Fraction LUTI
%Column 8: Fraction of unpaired TRITC
%Column 9: Fraction of unpaired Cy5

Filename_1 = 'SET 1'; %Reference data set
load('CHOOSE AN INTPUT FILEPATH FOR THE FIRST DATA SET');
Data_set_1 = All_counting_stats_raw;

Filename_2 = 'SET 2';
load('CHOOSE AN INTPUT FILEPATH FOR THE SECOND DATA SET');
Data_set_2 = All_counting_stats_raw;




resampled_set_1 = NaN(bootstrap_sample_size, size(All_counting_stats_raw,2),n_bootstrap);
resampled_set_2 = NaN(bootstrap_sample_size, size(All_counting_stats_raw,2),n_bootstrap);
Ranksum_bootstrap = NaN(3,n_bootstrap);

for number = 1:n_bootstrap
    if mod(number,100) == 0
        disp(['Iteration ', num2str(number), ' completed...']);
    end
    resampled_set_1(:,:,number) = datasample(Data_set_1,bootstrap_sample_size,1,'Replace',false);
    resampled_set_2(:,:,number) = datasample(Data_set_2,bootstrap_sample_size,1,'Replace',false);
    Ranksum_bootstrap(1,number) = ranksum(resampled_set_1(:,1,number),resampled_set_2(:,1,number));
    Ranksum_bootstrap(2,number) = ranksum(resampled_set_1(:,3,number),resampled_set_2(:,3,number));
    Ranksum_bootstrap(3,number) = ranksum(resampled_set_1(:,6,number),resampled_set_2(:,6,number));
end

resampled_set_1_mean = permute(mean(resampled_set_1,1),[3,2,1]);
resampled_set_2_mean = permute(mean(resampled_set_2,1),[3,2,1]);
resampled_set_1_median = permute(median(resampled_set_1,1),[3,2,1]);
resampled_set_2_median = permute(median(resampled_set_2,1),[3,2,1]);
resampled_set_1_stdev = permute(std(resampled_set_1,1,1),[3,2,1]);
resampled_set_2_stdev = permute(std(resampled_set_2,1,1),[3,2,1]);
Difference_in_means = resampled_set_2_mean - resampled_set_1_mean;
Difference_in_medians = resampled_set_2_median - resampled_set_1_median;
%Calculate the number of times the median from the reference was less
%than (or greater) than the median from the experiment. Compare this
%witht he Wilcoxon Rank Sum p-value (should be similar).
% rank_P_LUTI_left = ranksum(Data_set_1(:,1),Data_set_2(:,1),'alpha',0.01,...
%     'tail','left');
% rank_P_LUTI_right = ranksum(Data_set_1(:,1),Data_set_2(:,1),'alpha',0.01,...
%     'tail','right');
% rank_P_ORF_left = ranksum(Data_set_1(:,3),Data_set_2(:,3),'alpha',0.01,...
%     'tail','left');
% rank_P_ORF_right = ranksum(Data_set_1(:,3),Data_set_2(:,3),'alpha',0.01,...
%     'tail','right');
% rank_P_Total_left = ranksum(Data_set_1(:,6),Data_set_2(:,6),'alpha',0.01,...
%     'tail','left');
% rank_P_Total_right = ranksum(Data_set_1(:,6),Data_set_2(:,6),'alpha',0.01,...
%     'tail','right');
rank_P_LUTI = median(Ranksum_bootstrap(1,:));
rank_P_ORF = median(Ranksum_bootstrap(2,:));
rank_P_Total = median(Ranksum_bootstrap(3,:));

% pseudo_p_LUTI = [length(find(Difference_in_medians(:,1)<0));length(find(Difference_in_medians(:,1)>=0));rank_P_LUTI_left;rank_P_LUTI_right];
% pseudo_p_ORF = [length(find(Difference_in_medians(:,3)<0));length(find(Difference_in_medians(:,3)>=0));rank_P_ORF_left;rank_P_ORF_right];
% pseudo_p_total = [length(find(Difference_in_medians(:,6)<0));length(find(Difference_in_medians(:,6)>=0));rank_P_Total_left;rank_P_Total_right];
% TableOfVals = table(pseudo_p_LUTI,pseudo_p_ORF,pseudo_p_total,'RowNames',{'Diff < 0','Diff > 0', 'Ranksum left tail p-val','Ranksum right tail p-val'});

pseudo_p_LUTI = [length(find(Difference_in_medians(:,1)<0));length(find(Difference_in_medians(:,1)>=0));rank_P_LUTI];
pseudo_p_ORF = [length(find(Difference_in_medians(:,3)<0));length(find(Difference_in_medians(:,3)>=0));rank_P_ORF];
pseudo_p_total = [length(find(Difference_in_medians(:,6)<0));length(find(Difference_in_medians(:,6)>=0));rank_P_Total];
TableOfVals = table(pseudo_p_LUTI,pseudo_p_ORF,pseudo_p_total,'RowNames',{'Diff < 0','Diff > 0', 'Ranksum p-val'});

%%
writetable(TableOfVals,[Output_filepath,Case_name,'_table.txt'],'Delimiter','\t');

Compound_figure = figure('position',[20 100 2200 1900]);
subplot(1,3,1)
histogram(Difference_in_medians(:,1),[floor(min(Difference_in_medians(:,1)))-1:1:ceil(max(Difference_in_medians(:,1)))+1]);
plot_text{1} = 'Difference in median LUTI mRNA between';
plot_text{2} = [Filename_1,' and ',Filename_2];
title(plot_text);
text

subplot(1,3,2)
histogram(Difference_in_medians(:,3),[floor(min(Difference_in_medians(:,3)))-1:1:ceil(max(Difference_in_medians(:,3)))+1]);
plot_text{1} = 'Difference in median ORF mRNA between';
plot_text{2} = [Filename_1,' and ',Filename_2];
title(plot_text);

subplot(1,3,3)
histogram(Difference_in_medians(:,6),[floor(min(Difference_in_medians(:,6)))-1:1:ceil(max(Difference_in_medians(:,6)))+1]);
plot_text{1} = 'Difference in median number of total mRNA between';
plot_text{2} = [Filename_1,' and ',Filename_2];
title(plot_text);



img = getframe(gcf);
imwrite(img.cdata, [Output_filepath,Case_name,'_histograms.png']);




