%%%%%%%FISH-FINDER: A love story %%%%%%%%%%%%%%%%%
%%%%%%%David McSwiggen, Dec 2016

%%%%%%% This script uses gaussian detection/localization as applied by
%%%%%%% the SLIMfast implementation of MTT in order to detect spots coming
%%%%%%% from a maximum projection image of FISH data in multiple colors.
%%%%%%% Then it will attempt to pair up spots from the different images
%%%%%%% in order to determine colocalization.
%%%%%%% 
%%%%%%% To effectively use this code, one should have images in a .TIFF
%%%%%%% format, and have at least two separate FISH color channels as well
%%%%%%% as DAPI (or similar) to allow for the segmentation of cells.
%%%%%%% Filepaths to each image should be provided, and the code can take
%%%%%%% into account different replicates uding the "Position_number" and
%%%%%%% "Replicate_number" variables.
%%%%%%%
%%%%%%% After detection of all spots, a threshold must be used to determine
%%%%%%% which are truly mRNA and which are background noise. This threshold
%%%%%%% will depend on the parameters used to aquire the images, and may be
%%%%%%% very different from those used in the Chen, Tresenrider, et. al.
%%%%%%% manuscript.
%%%%%%%
%%%%%%% For more information, or to address specific questions, please
%%%%%%% email dmcswiggen@berkeley.edu or direct them to the corresponding
%%%%%%% authors of the Chen, Tresenrider, et. al. manuscript.


%% PARAMETERS
close all; clear; clc;
%%%%%%%%%%%%%%%%%%%% DEFINE INPUT AND OUTPUT PATHS %%%%%%%%%%%%%%%%%%%%%%%%
% specify input path with tiff files:
input_path=('/data2/David_M/');
Position_number = 9;
Replicate_number = 3;
DIC_file = '8144_6hr_30minDigest_2_1_R3D_REF.dv.8bit.tif';
DAPI_file='C3_MAX_8144_6hr_30minDigest_2_1_R3D.tif';
Cy5_file='C1_MAX_8144_6hr_30minDigest_2_1_R3D.tif';
TRITC_file='C2_MAX_8144_6hr_30minDigest_2_1_R3D.tif';
% DIC_file = ['DICmax_rep_',num2str(Replicate_number),'_pos_',num2str(Position_number), '.tif'];
% DAPI_file=['DAPImax_rep_',num2str(Replicate_number),'_pos_',num2str(Position_number), '.tif'];
% Cy5_file=['Cy5max_rep_',num2str(Replicate_number),'_pos_',num2str(Position_number), '.tif'];
% TRITC_file=['TRITCmax_rep_',num2str(Replicate_number),'_pos_',num2str(Position_number), '.tif'];


Genotype = '8144RTG';
condition = 'RTG60min';
output_path=('/data2/David_M/'); %Give filepath to EXISTING folder where you want to save the results
quantified_spot_path = ('/data2/David_M/');

% GIVE FILEPATH TO MTT ('SLIMFAST') CODE:
addpath('/data2/David_M/Jing/FISH_FINDER/SLIMFAST_batch_fordist');
addpath('/data2/David_M/Jing/FISH_FINDER/SLIMFAST_batch_fordist/bfmatlab');

%Which Operations will you perform?
Image_files_from = 1; %If 0, the image files are from ImageJ, if 1 they are from Matlab
Outline_cells = 0;
KernelFiltering = 1; % 0 = Do not perform, 1 = perform, 2 = perform and then break
ParticleDetect = 1;% 0 = Do not perform, 1 = perform, 2 = perform and then break
Localize_TRITC = 1;
Localize_Cy5 = 1;
Optimize_thresholds = 1;
Save_output = 1;


%Set initial detection parameters
number_to_test = 1000;
XY_kernel = 8;
DAPI_threshold = 2.2;

%% READ IN RAW IMAGES
Replicate_name = num2str(Replicate_number);
Position_name = num2str(Position_number);



%16-bit images from TIFF file
if Image_files_from == 1 %If the file was saved form Matlab
    original_image_DAPI = imread([input_path, DAPI_file]);
    original_image_DIC = imread([input_path,DIC_file]);
    original_image_TRITC = imread([input_path,TRITC_file]);
    original_image_Cy5 = imread([input_path,Cy5_file]);
else
    image_DAPI_struct = tiffread([input_path, DAPI_file]);
    image_TRITC_struct = tiffread([input_path,TRITC_file]);
    image_Cy5_struct = tiffread([input_path,Cy5_file]);
    image_DIC_struct = tiffread([input_path,DIC_file]);
    original_image_DAPI = image_DAPI_struct.data;
    original_image_DIC = image_DIC_struct.data;
    original_image_TRITC =image_TRITC_struct.data;
    original_image_Cy5 =image_Cy5_struct.data;
end

%Convert to "Double" format
original_image_DAPI_matrix = double(original_image_DAPI) +1;
original_image_DIC_matrix = double(original_image_DIC) +1;
original_image_TRITC_matrix = double(original_image_TRITC) +1;
original_image_Cy5_matrix = double(original_image_Cy5) +1;

if Save_output == 1
    [s,mess,messid] = mkdir(output_path, Genotype);
    [s,mess,messid] = mkdir([output_path,'/', Genotype], condition);
    [s,mess,messid] = mkdir([output_path,'/', Genotype,'/',condition], Replicate_name);
    [s,mess,messid] = mkdir([output_path,'/', Genotype,'/',condition,'/',Replicate_name],Position_name);
    imwrite(original_image_DIC,...
        [output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
        'DICmax_rep_',Replicate_name,'_pos_',num2str(Position_number),'.tif']);
    imwrite(original_image_DAPI,...
        [output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
        'DAPImax_rep_',Replicate_name,'_pos_',num2str(Position_number),'.tif']);
    imwrite(original_image_TRITC,...
        [output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
        'TRITCmax_rep_',Replicate_name,'_pos_',num2str(Position_number),'.tif']);
    imwrite(original_image_Cy5,...
        [output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
        'Cy5max_rep_',Replicate_name,'_pos_',num2str(Position_number),'.tif']);
    
end

%% GENERATE DAPI OUTLINES

gauss_DAPI_image = imgaussfilt(original_image_DAPI_matrix./imgaussfilt(original_image_DAPI_matrix,50),3);
DAPI_binary = bwareaopen((gauss_DAPI_image > DAPI_threshold),350);
nucleus_boundaries = bwboundaries(DAPI_binary);
labeled_nuclei = bwlabel(DAPI_binary, 8);

%colored_nuclei = label2rgb(labeled_nuclei, 'parula', 'k', 'shuffle');


nuclei_measurements = regionprops(labeled_nuclei, original_image_DAPI, 'all');
nuclei_outlines = bwboundaries(DAPI_binary,8);
number_of_nuclear_boundaries = size(nuclei_outlines,1);
number_of_nuclei = size(nuclei_measurements,1);
All_nuclei_centroids = [nuclei_measurements.Centroid];
centroidsX = All_nuclei_centroids(1:2:end-1);
centroidsY = All_nuclei_centroids(2:2:end);

figure('position',[20 100 1600 1200]);
s1 = subplot(2,3,1);
imshow(original_image_DAPI,[]);
hold on;
for point = 1 : number_of_nuclear_boundaries
    NuclearBoundary = nuclei_outlines{point};
    plot(NuclearBoundary(:,2), NuclearBoundary(:,1), 'g', 'LineWidth', 2);
end

for nuc = 1:number_of_nuclei
    text(centroidsX(nuc) + 20, centroidsY(nuc), num2str(nuc), 'FontSize', 10, 'FontWeight', 'Bold','Color','g');
end

title ('DIC image with DAPI outlines superimposed');
impixelinfo;

s2 = subplot(2,3,2);
imshow(original_image_TRITC,[]);
title ('TRITC FISH');
impixelinfo;

s3 = subplot(2,3,3);
imshow(original_image_Cy5,[]);
title ('Cy5 FISH');
impixelinfo;

subplot(2,3,4)
[pixelCount_DAPI, grayLevels] = imhist(original_image_DIC);
bar(pixelCount_DAPI);
title('Pixel intensity', 'FontSize', 12);
xlim([0 50]); % Scale x axis manually.
grid on;

subplot(2,3,6)
[pixelCount_Cy5, grayLevels] = imhist(original_image_Cy5);
bar(pixelCount_Cy5);
title('Pixel intensity', 'FontSize', 12);
xlim([0 50]); % Scale x axis manually.
grid on;

subplot(2,3,5)
[pixelCount_TRITC, grayLevels] = imhist(original_image_TRITC);
bar(pixelCount_TRITC);
title('Pixel intensity', 'FontSize', 12);
xlim([0 50]); % Scale x axis manually.
grid on;
linkaxes([s1,s2,s3],'xy');


message = sprintf('Inspect the images to check whether the DAPI segmentation is correct.');
uiwait(msgbox(message));



%% GENERATE CELL OUTLINES
if Outline_cells == 1
    seg_prompt = {'How many cells would you like to outline?'};
    seg_prompt_name = 'Wait for user';
    numlines = 1;
    seg_defaults = {'1'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    seg_answer = inputdlg(seg_prompt,seg_prompt_name,numlines,seg_defaults);
    close;
    
    FISH_Outlines_Number = str2num(seg_answer{1});
    FISH_Outlines_masks = zeros(size(original_image_DAPI_matrix,1),size(original_image_DAPI_matrix,2),FISH_Outlines_Number);
    FISH_Outlines_props = zeros(3,FISH_Outlines_Number);
    Temp_nuclei = zeros(size(original_image_DAPI_matrix,1),size(original_image_DAPI_matrix,2));
    FISH_data_individual_cells = struct;
    FISH_data_individual_cells.fieldname = [Genotype,'_',condition,'_rep_',Replicate_name,'_pos_',num2str(Position_number)];
    
    %%%%% CELL TYPE MASKS %%%%%%%%
    iter = 1;
    message = sprintf('Draw ROIs cells you would like to quantitate');
    uiwait(msgbox(message));
    for region = 1:FISH_Outlines_Number
        FISH_data_individual_cells(region).id = region;
        Outline_masks_all(:,:) = sum(FISH_Outlines_masks,3)>0;
        Outline_masks_all(:,:) = imfill(Outline_masks_all(:,:),'holes');
        Cell_outlines = bwboundaries(Outline_masks_all(:,:),8);
        Number_of_cell_boundaries = size(Cell_outlines, 1);
        
        
        imshow(original_image_DIC,[]);
        axis image;
        hold on;
        impixelinfo;
        for point = 1 : Number_of_cell_boundaries
            CellBoundary = Cell_outlines{point};
            plot(CellBoundary(:,2), CellBoundary(:,1), 'g', 'LineWidth', 3);
        end
        for point = 1 : number_of_nuclear_boundaries
            NuclearBoundary = nuclei_outlines{point};
            plot(NuclearBoundary(:,2), NuclearBoundary(:,1), 'b', 'LineWidth', 1);
        end
        for cell = 1:region
            text(FISH_data_individual_cells(cell).centroid(1), FISH_data_individual_cells(cell).centroid(2), num2str(FISH_data_individual_cells(cell).id), 'FontSize', 10, 'FontWeight', 'Bold','Color','g');
        end
        hold off;
        
        [FISH_Outlines_masks(:,:,region), x_coord, y_coord] = roipoly();
        FISH_data_individual_cells(region).outline = [x_coord, y_coord];
        Mask_data = regionprops(FISH_Outlines_masks(:,:,region));
        if ~isempty(Mask_data)
            FISH_data_individual_cells(region).centroid = Mask_data.Centroid;
        else
            FISH_data_individual_cells(region).centroid = [1,1]
        end
        Temp_nuclei(:,:) = and(FISH_Outlines_masks(:,:,region),DAPI_binary);
        Nucleus_per_cell_BW = bwboundaries(Temp_nuclei(:,:),8);
        if size(Nucleus_per_cell_BW) == 1
            Nucleus_per_cell_binary(:,:,region) = Temp_nuclei(:,:);
            FISH_data_individual_cells(region).nucleus_outline = [Nucleus_per_cell_BW{1}(:,1), Nucleus_per_cell_BW{1}(:,2)];
        else
            FISH_data_individual_cells(region).nucleus_outline = [];
            message = sprintf('More than one nucleus was detected within this ROI. \n No nuclei will be recorded');
            uiwait(msgbox(message));
        end
        
        if region == FISH_Outlines_Number
            close;
        else
            message = sprintf('Outline another region');
            uiwait(msgbox(message));
        end
        iter = iter+1;
    end
    
    Outline_masks_all(:,:) = sum(FISH_Outlines_masks,3)>0;
    Cell_outlines = bwboundaries(Outline_masks_all(:,:),8);
    Number_of_cell_boundaries = size(Cell_outlines, 1);
    DAPI_mask_all = sum(Nucleus_per_cell_binary,3) > 0;
    DAPI_outlines_annotated = bwboundaries(DAPI_mask_all,8);
    Number_of_annotated_nuclei = size(DAPI_outlines_annotated,1);
    
    figure
    imshow(original_image_Cy5,[]);
    axis image;
    hold on;
    for point = 1 : Number_of_cell_boundaries
        CellBoundary = Cell_outlines{point};
        plot(CellBoundary(:,2), CellBoundary(:,1), 'g', 'LineWidth', 3);
    end
    
    for point = 1 : Number_of_annotated_nuclei
        NucBoundary = DAPI_outlines_annotated{point};
        plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
    end
    hold off
    
    if Save_output ==1
        save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
            'outlines_rep_',Replicate_name,'_pos_',num2str(Position_number)],...
            'FISH_Outlines_masks', 'FISH_Outlines_Number','FISH_data_individual_cells');
        save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
            '/','Nuclear_outlines_',Replicate_name,'_pos_',num2str(Position_number)],...
            'nuclei_outlines', 'number_of_nuclei');
    end
    
    
else
    load([output_path,Genotype,'/',condition,'/',Replicate_name,'/',Position_name,'/',...
        'outlines_rep_',Replicate_name,'_pos_',num2str(Position_number),'.mat']);
    FISH_Outlines_Number = size(FISH_Outlines_masks,3);
    % Apply DAPI mask to cells that were annotated before the function was
    % added
    iter = 1;
    for region = 1:FISH_Outlines_Number
        
        Mask_data = regionprops(FISH_Outlines_masks(:,:,region));
        if ~isempty(Mask_data)
            FISH_data_individual_cells(region).centroid = Mask_data.Centroid;
        else
            FISH_data_individual_cells(region).centroid = [1,1]
        end
        Temp_nuclei(:,:) = and(FISH_Outlines_masks(:,:,region),DAPI_binary);
        Nucleus_per_cell_BW = bwboundaries(Temp_nuclei(:,:),8);
        if size(Nucleus_per_cell_BW) == 1
            Nucleus_per_cell_binary(:,:,region) = Temp_nuclei(:,:);
            FISH_data_individual_cells(region).nucleus_outline = [Nucleus_per_cell_BW{1}(:,1), Nucleus_per_cell_BW{1}(:,2)];
        else
            FISH_data_individual_cells(region).nucleus_outline = [];
        end
        
        if region == FISH_Outlines_Number
            close
        end
        iter = iter+1;
    end
    
    Outline_masks_all(:,:) = sum(FISH_Outlines_masks,3)>0;
    Cell_outlines = bwboundaries(Outline_masks_all(:,:),8);
    Number_of_cell_boundaries = size(Cell_outlines, 1);
    DAPI_mask_all = sum(Nucleus_per_cell_binary,3) > 0;
    DAPI_outlines_annotated = bwboundaries(DAPI_mask_all,8);
    Number_of_annotated_nuclei = size(DAPI_outlines_annotated,1);
    
    figure
    imshow(original_image_Cy5,[]);
    axis image;
    hold on;
    for point = 1 : Number_of_cell_boundaries
        CellBoundary = Cell_outlines{point};
        plot(CellBoundary(:,2), CellBoundary(:,1), 'g', 'LineWidth', 3);
    end
    
    for point = 1 : Number_of_annotated_nuclei
        NucBoundary = DAPI_outlines_annotated{point};
        plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
    end
    for cell = 1:region
        text(FISH_data_individual_cells(cell).centroid(1), FISH_data_individual_cells(cell).centroid(2), num2str(FISH_data_individual_cells(cell).id), 'FontSize', 12, 'FontWeight', 'Bold','Color','w');
    end
    hold off
    
    if Save_output ==1
        save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
            'outlines_rep_',Replicate_name,'_pos_',num2str(Position_number)],...
            'FISH_Outlines_masks', 'FISH_Outlines_Number','FISH_data_individual_cells');
    end
end


message = sprintf('Check to make sure these annotations look correct.\nIf so, click OK. Otherwise abort the script and begin again.');
uiwait(msgbox(message));
close;

%% APPLY GAUSSIAN FILTERING

if KernelFiltering == 1
    TRITC_background = imgaussfilt(original_image_TRITC_matrix,XY_kernel);
    TRITC_image_use_this = (original_image_TRITC_matrix./TRITC_background).*original_image_TRITC_matrix;
    if KernelFiltering ==1
        TRITC_detect_image = TRITC_image_use_this;
    else
        TRITC_detect_image =original_image_TRITC_matrix;
    end
    
    Cy5_background = imgaussfilt(original_image_Cy5_matrix,(XY_kernel*1.3));
    Cy5_image_use_this = (original_image_Cy5_matrix./Cy5_background).*original_image_Cy5_matrix;
    if KernelFiltering ==1
        Cy5_detect_image = Cy5_image_use_this;
    else
        Cy5_detect_image =original_image_Cy5_matrix;
    end
    
    figure('position',[20 100 1600 1200])
    h1 = subplot(2,2,1);
    imshow(original_image_TRITC_matrix,[]);
    title('Unfiltered TRITC');
    impixelinfo;
    
    h2 = subplot(2,2,2);
    imshow(original_image_Cy5_matrix,[]);
    title('Unfiltered Cy5');
    impixelinfo;
    
    h3 = subplot(2,2,3);
    imshow(TRITC_detect_image,[0 max(max(TRITC_image_use_this))/4]);
    title('Kernel filtered TRITC');
    impixelinfo;
    
    h4 = subplot(2,2,4);
    imshow(Cy5_detect_image,[0 max(max(Cy5_image_use_this))/2]);
    title('Kernel filtered Cy5');
    impixelinfo;
    linkaxes([h1,h2,h3, h4],'xy');
    
    message = sprintf('Filtered images to detect mRNA from.\nIf these look correct, click OK.\nOtherwise abort the script and begin again.');
    uiwait(msgbox(message));
    close;
    
elseif KernelFiltering == 2;
    error('Outline next file')
    close all;
else
    TRITC_image_use_this = original_image_TRITC;
    Cy5_image_use_this = original_image_Cy5;
end



%% PARTICLE LOCALIZATION
if ParticleDetect > 0
    LocalizationError_TRITC = -0.1; % Localization Error: -6 = 10^-6
    LocalizationError_Cy5 = -0.1; % Localization Error: -6 = 10^-6
    EmissionWavelength_TRITC = 600; % wavelength in nm; consider emission max and filter cutoff
    EmissionWavelength_Cy5 = 650; % wavelength in nm; consider emission max and filter cutoff
    NumDeflationLoops = 0; % Generaly keep this to 0;
    %In this section, the code built in SLIMfast should be extracted, and
    %all the necessary parameters (aquisition and localization) set so
    %that the script can run a localization using MTT algorithm on the
    %image.
    
    %%% DEFINE STRUCTURED ARRAY WITH ALL THE MTT SETTINGS FOR LOCALIZATION %%%
    % add the required functions to the path:
    
    
    %INPUT: (dir_path,impars, locpars, imstack) => Directory info, Image (aquisition) parameters, localization parameters, Filepate to tiff stack
    if Localize_TRITC == 1
        tic;
        disp('Detecting spots in the TRITC image');
        % imaging parameters
        impars.PixelSize=0.160; % um per pixel
        impars.psf_scale=1.35; % PSF scaling
        impars.NA=1.40; % NA of detection objective
        impars.FrameRate= 1000; %secs %Not sure if this is actually in important input for localization, but don't want to mess with it. Nonesense number
        impars.FrameSize= 1000; %secs
        
        % localization parameters
        locpars.wn=5; %detection box in pixels
        locpars.errorRate= LocalizationError_TRITC; % error rate (10^-)
        locpars.dfltnLoops= NumDeflationLoops; % number of deflation loops
        locpars.minInt=0; %minimum intensity in counts
        locpars.maxOptimIter= 50; % max number of iterations
        locpars.termTol= -2; % termination tolerance
        locpars.isRadiusTol=false; % use radius tolerance
        locpars.radiusTol=50; % radius tolerance in percent
        locpars.posTol= 1.5;%max position refinement
        locpars.optim = [locpars.maxOptimIter,locpars.termTol,locpars.isRadiusTol,locpars.radiusTol,locpars.posTol];
        locpars.isThreshLocPrec = false;
        locpars.minLoc = 0;
        locpars.maxLoc = inf;
        locpars.isThreshSNR = false;
        locpars.minSNR = 0;
        locpars.maxSNR = inf;
        locpars.isThreshDensity = false;
        
        impars.name = TRITC_file;
        impars.wvlnth= EmissionWavelength_TRITC/1000; %emission wavelength in um
        impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels
        Localizations_raw_TRITC = localizeParticles_ASH(input_path,impars, locpars, TRITC_detect_image);
        %data structured array: % ctrsX = x-coord of all localizations
        % ctrsY = y-coord of all localizations
        % ctrsN = number of localizations in a given
        % frame
        % signal = signal (obvious)
        % frame (In this case, subsquare)
        
        Localization_matrix_TRITC = zeros(Localizations_raw_TRITC.ctrsN,6);
        for line=1:Localizations_raw_TRITC.ctrsN
            Localization_matrix_TRITC(line,1) = Localizations_raw_TRITC.ctrsX(line) + 1.25;
            Localization_matrix_TRITC(line,2) = Localizations_raw_TRITC.ctrsY(line) + 1.25;
            Localization_matrix_TRITC(line,3) = Localizations_raw_TRITC.signal(line);
            Localization_matrix_TRITC(line,4) = Localizations_raw_TRITC.noise(line);
            Localization_matrix_TRITC(line,5) = Localizations_raw_TRITC.offset(line);
            Localization_matrix_TRITC(line,6) = Localizations_raw_TRITC.signal(line) / Localizations_raw_TRITC.noise(line);
        end
        
        for ROI = 1:FISH_Outlines_Number
            [spots_inside_cell_index,spots_edge_cell_index] = inpolygon(Localization_matrix_TRITC(:,1),Localization_matrix_TRITC(:,2),...
                FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
            FISH_data_individual_cells(ROI).TRITC_spots_all = Localization_matrix_TRITC(spots_inside_cell_index,1:6);
            if size(FISH_data_individual_cells(ROI).nucleus_outline > 0)
                [spots_inside_nucleus_index,spots_edge_nucleus_index] = inpolygon(Localization_matrix_TRITC(:,1),Localization_matrix_TRITC(:,2),...
                    FISH_data_individual_cells(ROI).nucleus_outline(:,1),FISH_data_individual_cells(ROI).nucleus_outline(:,2));
                FISH_data_individual_cells(ROI).TRITC_spots_nucleus = Localization_matrix_TRITC(spots_inside_nucleus_index,1:6);
            else
                FISH_data_individual_cells(ROI).TRITC_spots_nucleus = [];
            end
        end
        
        toc;
        if Save_output == 1
            save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
                'TRITC_data_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Localization_matrix_TRITC');
        end
    end
    
    
    
    if Localize_Cy5 == 1
        tic;
        disp('Detecting spots in the Cy5 image');
        % imaging parameters
        impars.PixelSize=0.160; % um per pixel
        impars.psf_scale=1.35; % PSF scaling
        impars.NA=1.40; % NA of detection objective
        impars.FrameRate= 1000; %secs %Not sure if this is actually in important input for localization, but don't want to mess with it. Nonesense number
        impars.FrameSize= 1000; %secs
        
        % localization parameters
        locpars.wn=5; %detection box in pixels
        locpars.errorRate= LocalizationError_Cy5; % error rate (10^-)
        locpars.dfltnLoops= NumDeflationLoops; % number of deflation loops
        locpars.minInt=0; %minimum intensity in counts
        locpars.maxOptimIter= 50; % max number of iterations
        locpars.termTol= -2; % termination tolerance
        locpars.isRadiusTol=false; % use radius tolerance
        locpars.radiusTol=50; % radius tolerance in percent
        locpars.posTol= 1.5;%max position refinement
        locpars.optim = [locpars.maxOptimIter,locpars.termTol,locpars.isRadiusTol,locpars.radiusTol,locpars.posTol];
        locpars.isThreshLocPrec = false;
        locpars.minLoc = 0;
        locpars.maxLoc = inf;
        locpars.isThreshSNR = false;
        locpars.minSNR = 0;
        locpars.maxSNR = inf;
        locpars.isThreshDensity = false;
        impars.name = Cy5_file;
        impars.wvlnth= EmissionWavelength_Cy5/1000; %emission wavelength in um
        impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels
        Localizations_raw_Cy5 = localizeParticles_ASH(input_path,impars, locpars, Cy5_detect_image);
        %data structured array: % ctrsX = x-coord of all localizations
        % ctrsY = y-coord of all localizations
        % ctrsN = number of localizations in a given
        % frame
        % signal = signal (obvious)
        % frame (In this case, subsquare)
        toc;
        
        Localization_matrix_Cy5 = zeros(Localizations_raw_Cy5.ctrsN,6);
        for line=1:Localizations_raw_Cy5.ctrsN
            Localization_matrix_Cy5(line,1) = Localizations_raw_Cy5.ctrsX(line)+1;
            Localization_matrix_Cy5(line,2) = Localizations_raw_Cy5.ctrsY(line)+1;
            Localization_matrix_Cy5(line,3) = Localizations_raw_Cy5.signal(line);
            Localization_matrix_Cy5(line,4) = Localizations_raw_Cy5.noise(line);
            Localization_matrix_Cy5(line,5) = Localizations_raw_Cy5.offset(line);
            Localization_matrix_Cy5(line,6) = Localizations_raw_Cy5.signal(line) / Localizations_raw_Cy5.noise(line);
        end
        
        for ROI = 1:FISH_Outlines_Number
            [spots_inside_cell_index,spots_edge_cell_index] = inpolygon(Localization_matrix_Cy5(:,1),Localization_matrix_Cy5(:,2),...
                FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
            FISH_data_individual_cells(ROI).Cy5_spots_all = Localization_matrix_Cy5(spots_inside_cell_index,1:6);
            
            if size(FISH_data_individual_cells(ROI).nucleus_outline > 0)
                [spots_inside_nucleus_index,spots_edge_nucleus_index] = inpolygon(Localization_matrix_Cy5(:,1),Localization_matrix_Cy5(:,2),...
                    FISH_data_individual_cells(ROI).nucleus_outline(:,1),FISH_data_individual_cells(ROI).nucleus_outline(:,2));
                FISH_data_individual_cells(ROI).Cy5_spots_nucleus = Localization_matrix_Cy5(spots_inside_nucleus_index,1:6);
            else
                FISH_data_individual_cells(ROI).Cy5_spots_nucleus = [];
            end
        end
        
        toc;
        if Save_output == 1
            save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
                'Cy5_data_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Localization_matrix_Cy5');
        end
        
        
    end
    if ParticleDetect == 2
        error('Outline next file')
        close all;
    end
    
else
    load([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
        'TRITC_data_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Localization_matrix_TRITC');
    load([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
        'Cy5_data_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Localization_matrix_Cy5');
    
end


%% Threshold Setting
%%%% Sort spots into lists based on Signal, or SNR to plot CDFs
%%%% First threshold based on spot intensity, then re-plot and threshold
%%%% based on SNR.
if Optimize_thresholds == 1
    %Signal first pass
    TRITC_mRNA_signal_threshold = [];
    TRITC_signal_raw = sort(Localization_matrix_TRITC(:,3));
    increment = round((TRITC_signal_raw(end)-TRITC_signal_raw(1))/number_to_test);
    iter = 1;
    for threshold = round(TRITC_signal_raw(1)):increment:round(TRITC_signal_raw(end))
        boolean_threshold_applied = TRITC_signal_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        TRITC_mRNA_signal_threshold(iter,1) = threshold;
        TRITC_mRNA_signal_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    Cy5_mRNA_signal_threshold = [];
    Cy5_signal_raw = sort(Localization_matrix_Cy5(:,3));
    increment = round((Cy5_signal_raw(end)-Cy5_signal_raw(1))/number_to_test);
    iter = 1;
    for threshold = round(Cy5_signal_raw(1)):increment:round(Cy5_signal_raw(end))
        boolean_threshold_applied = Cy5_signal_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        Cy5_mRNA_signal_threshold(iter,1) = threshold;
        Cy5_mRNA_signal_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    %SNR Second pass
    TRITC_SNR_CDF = [];
    TRITC_mRNA_SNR_threshold = [];
    TRITC_SNR_raw = sort(Localization_matrix_TRITC(:,6));
    iter = 1;
    for threshold = floor(TRITC_SNR_raw(1)):ceil(TRITC_SNR_raw(end))
        
        boolean_threshold_applied = TRITC_SNR_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        TRITC_mRNA_SNR_threshold(iter,1) = threshold;
        TRITC_mRNA_SNR_threshold(iter,2) = number_above;
        TRITC_SNR_CDF(iter,1) = length(Localization_matrix_TRITC(:,3)) - number_above;
        
        iter = iter+1;
    end
    
    Cy5_SNR_CDF = [];
    Cy5_mRNA_SNR_threshold = [];
    Cy5_SNR_raw = sort(Localization_matrix_Cy5(:,6));
    iter = 1;
    for threshold = floor(Cy5_SNR_raw(1)):ceil(Cy5_SNR_raw(end))
        
        boolean_threshold_applied = Cy5_SNR_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        Cy5_mRNA_SNR_threshold(iter,1) = threshold;
        Cy5_mRNA_SNR_threshold(iter,2) = number_above;
        Cy5_SNR_CDF(iter,1) = length(Localization_matrix_Cy5(:,3)) - number_above;
        
        iter = iter+1;
    end
    
    figure('position',[20 100 1600 1200])
    subplot(2,3,1)
    plot(TRITC_mRNA_signal_threshold(:,1),TRITC_mRNA_signal_threshold(:,2))
    title('TRITC mRNA count as a function of signal intensity')
    
    subplot(2,3,2)
    plot(TRITC_mRNA_SNR_threshold(:,1),TRITC_mRNA_SNR_threshold(:,2))
    title('TRITC mRNA count as a function of Signal to Noise ratio')
    
    subplot(2,3,3)
    semilogx(Localization_matrix_TRITC(:,3), Localization_matrix_TRITC(:,6),'o')
    title('Cy5 mRNA count as a function of signal intensity')
    
    subplot(2,3,4)
    plot(Cy5_mRNA_signal_threshold(:,1),Cy5_mRNA_signal_threshold(:,2))
    title('Cy5 mRNA count as a function of signal intensity')
    
    subplot(2,3,5)
    plot(Cy5_mRNA_SNR_threshold(:,1),Cy5_mRNA_SNR_threshold(:,2))
    title('Cy5 mRNA count as a function of Signal to Noise ratio')
    
    subplot(2,3,6)
    semilogx(Localization_matrix_Cy5(:,3), Localization_matrix_Cy5(:,6),'o')
    title('Cy5 mRNA count as a function of signal intensity')
    
    
    message = sprintf('Inspect curves to identify Signal intensity\n and Signal-to-Noise thresholds for both colors');
    uiwait(msgbox(message));
    
    seg_prompt = {'TRITC Signal threshold (x10^3)','Cy5 Signal threshold(x10^3)', 'TRITC SNR threshold','Cy5 SNR threshold'};
    seg_prompt_name = 'Wait for user';
    numlines = 1;
    seg_defaults = {'1.1','1.2','2','2'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    seg_answer = inputdlg(seg_prompt,seg_prompt_name,numlines,seg_defaults);
    close;
    
    Threshold_params = [str2num(seg_answer{1})*(10^3),str2num(seg_answer{2})*(10^3),str2num(seg_answer{3}),str2num(seg_answer{4})];
    
    TRITC_boolean = and(Localization_matrix_TRITC(:,3) > Threshold_params(1), Localization_matrix_TRITC(:,6) > Threshold_params(3));
    Thresholded_TRITC_matrix = Localization_matrix_TRITC(TRITC_boolean,1:6);
    Cy5_boolean = and(Localization_matrix_Cy5(:,3) > Threshold_params(2), Localization_matrix_Cy5(:,6) > Threshold_params(4));
    Thresholded_Cy5_matrix = Localization_matrix_Cy5(Cy5_boolean,1:6);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%% Threshold a second time to refine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    refine_threshold = 1;
    
    while refine_threshold == 1
        
        TRITC_mRNA_signal_threshold = [];
        TRITC_signal_raw = sort(Localization_matrix_TRITC(:,3));
        increment = round((TRITC_signal_raw(end)-TRITC_signal_raw(1))/number_to_test);
        iter = 1;
        for threshold = round(TRITC_signal_raw(1)):increment:round(TRITC_signal_raw(end))
            boolean_threshold_applied = TRITC_signal_raw > threshold;
            boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
            number_above = length(boolean_threshold_applied);
            TRITC_mRNA_signal_threshold(iter,1) = threshold;
            TRITC_mRNA_signal_threshold(iter,2) = number_above;
            iter = iter+1;
        end
        
        Cy5_mRNA_signal_threshold = [];
        Cy5_signal_raw = sort(Localization_matrix_Cy5(:,3));
        increment = round((Cy5_signal_raw(end)-Cy5_signal_raw(1))/number_to_test);
        iter = 1;
        for threshold = round(Cy5_signal_raw(1)):increment:round(Cy5_signal_raw(end))
            boolean_threshold_applied = Cy5_signal_raw > threshold;
            boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
            number_above = length(boolean_threshold_applied);
            Cy5_mRNA_signal_threshold(iter,1) = threshold;
            Cy5_mRNA_signal_threshold(iter,2) = number_above;
            iter = iter+1;
        end
        
        %SNR
        TRITC_mRNA_SNR_threshold = [];
        TRITC_SNR_raw = sort(Thresholded_TRITC_matrix(:,6));
        iter = 1;
        for threshold = floor(TRITC_SNR_raw(1)):ceil(TRITC_SNR_raw(end))
            
            boolean_threshold_applied = TRITC_SNR_raw > threshold;
            boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
            number_above = length(boolean_threshold_applied);
            TRITC_mRNA_SNR_threshold(iter,1) = threshold;
            TRITC_mRNA_SNR_threshold(iter,2) = number_above;
            iter = iter+1;
        end
        
        Cy5_mRNA_SNR_threshold = [];
        Cy5_SNR_raw = sort(Thresholded_Cy5_matrix(:,6));
        iter = 1;
        for threshold = floor(Cy5_SNR_raw(1)):ceil(Cy5_SNR_raw(end))
            
            boolean_threshold_applied = Cy5_SNR_raw > threshold;
            boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
            number_above = length(boolean_threshold_applied);
            Cy5_mRNA_SNR_threshold(iter,1) = threshold;
            Cy5_mRNA_SNR_threshold(iter,2) = number_above;
            iter = iter+1;
        end
        
        %%% PLOT QC statistics
        figure('position',[20 100 1600 1200])
        subplot(2,3,1)
        plot(TRITC_mRNA_signal_threshold(:,1),TRITC_mRNA_signal_threshold(:,2))
        maxYValue = ylim;
        hold on;
        title('TRITC mRNA signal after threshold')
        plot([Threshold_params(1), Threshold_params(1)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
        % Place a text label on the bar chart showing the threshold.
        annotationText = sprintf('Thresholded at %d x10^3 gray levels', str2num(seg_answer{1}));
        text(double(Threshold_params(1) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
        
        
        subplot(2,3,2)
        plot(TRITC_mRNA_SNR_threshold(:,1),TRITC_mRNA_SNR_threshold(:,2))
        title('TRITC mRNA SNR after threshold')
        hold on;
        maxYValue = ylim;
        plot([Threshold_params(3), Threshold_params(3)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
        annotationText = sprintf('Thresholded at %d ', str2num(seg_answer{3}));
        text(double(Threshold_params(3) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
        
        
        subplot(2,3,3)
        semilogx(Localization_matrix_TRITC(:,3), Localization_matrix_TRITC(:,6),'o')
        title('TRITC Spot intensity vs SNR, after threshold')
        hold on;
        maxXValue = xlim;
        plot([Threshold_params(1), Threshold_params(1)], [maxYValue(1),max(Localization_matrix_TRITC(:,6))*1.1], 'Color', 'r');
        plot([maxXValue(1),maxXValue(2)],[Threshold_params(3), Threshold_params(3)], 'Color', 'r');
        hold off;
        
        subplot(2,3,4)
        plot(Cy5_mRNA_signal_threshold(:,1),Cy5_mRNA_signal_threshold(:,2))
        title('Cy5 mRNA count as a function of signal intensity')
        maxYValue = ylim;
        hold on;
        title('Cy5 mRNA signal after threshold')
        plot([Threshold_params(2), Threshold_params(2)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
        % Place a text label on the bar chart showing the threshold.
        annotationText = sprintf('Thresholded at %d x10^3 gray levels', str2num(seg_answer{2}));
        text(double(Threshold_params(2) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
        
        subplot(2,3,5)
        plot(Cy5_mRNA_SNR_threshold(:,1),Cy5_mRNA_SNR_threshold(:,2))
        title('Cy5 mRNA count as a function of Signal to Noise ratio')
        maxYValue = ylim;
        hold on;
        title('Cy5 mRNA SNR after threshold')
        plot([Threshold_params(4), Threshold_params(4)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
        % Place a text label on the bar chart showing the threshold.
        annotationText = sprintf('Thresholded at %d ', str2num(seg_answer{4}));
        text(double(Threshold_params(4)*1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
        
        subplot(2,3,6)
        semilogx(Localization_matrix_Cy5(:,3), Localization_matrix_Cy5(:,6),'o')
        title('Cy5 Spot intensity vs SNR, after threshold')
        hold on;
        maxXValue = xlim;
        plot([Threshold_params(2), Threshold_params(2)], [maxYValue(1),max(Localization_matrix_Cy5(:,6))*1.1], 'Color', 'r');
        plot([maxXValue(1),maxXValue(2)], [Threshold_params(4), Threshold_params(4)],'Color', 'r');
        hold off;
        
        
        for ROI = 1:FISH_Outlines_Number
            [in,on] = inpolygon(Thresholded_TRITC_matrix(:,1),Thresholded_TRITC_matrix(:,2),...
                FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
            FISH_data_individual_cells(ROI).TRITC_spots_threshold =Thresholded_TRITC_matrix(in,1:6);
            
            [in,on] = inpolygon(Thresholded_Cy5_matrix(:,1),Thresholded_Cy5_matrix(:,2),...
                FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
            FISH_data_individual_cells(ROI).Cy5_spots_threshold =Thresholded_Cy5_matrix(in,1:6);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%% PAIR UP SPOTS USING KNN ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%
        for ROI = 1:FISH_Outlines_Number
            
            %By constraining this to per-cell analysis, can limit the number and
            %distance of potential neighbors and limit edge problems.
            if and(~isempty(FISH_data_individual_cells(ROI).Cy5_spots_threshold),~isempty(FISH_data_individual_cells(ROI).TRITC_spots_threshold))
                %By constraining this to per-cell analysis, can limit the number and
                %distance of potential neighbors and limit edge problems.
                
                [TRITC_knn_to_Cy5, dist_TRITC_to_Cy5] = knnsearch(FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:2),...
                    FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:2),'K',1);
                [Cy5_knn_to_TRITC, dist_Cy5_to_TRITC] = knnsearch(FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:2),...
                    FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:2),'K',1);
                
                %%% Pairing up spot first by checking for "unambiguous pairings", meaning
                %%% that if A thinks it's nearest neighbor is B, and if B thinks its
                %%% nearest neighbor is A, this is unambiguous. Ovbiously this can only be
                %%% true for the smaller of the two lists.
                TRITC_to_Cy5_pairs = zeros(length(TRITC_knn_to_Cy5(:,1)),2);
                for row = 1:length(TRITC_knn_to_Cy5(:,1))
                    TRITC_to_Cy5_pairs(row,1:2) = [row, TRITC_knn_to_Cy5(row)];
                end
                
                Cy5_to_TRITC_pairs = zeros(length(Cy5_knn_to_TRITC(:,1)),2);
                for row = 1:length(Cy5_knn_to_TRITC(:,1))
                    Cy5_to_TRITC_pairs(row,1:2) = [row, Cy5_knn_to_TRITC(row)];
                end
                
                Cy5_pair_boolean = zeros(length(Cy5_knn_to_TRITC(:,1)),1);
                TRITC_pair_boolean = zeros(length(TRITC_knn_to_Cy5(:,1)),1);
                
                if length(Cy5_knn_to_TRITC(:,1)) > length(TRITC_knn_to_Cy5(:,1))
                    % If there are fewer TRITC detections that Cy5
                    for mRNA = 1:length(TRITC_knn_to_Cy5(:,1))
                        if [TRITC_to_Cy5_pairs(mRNA,1),TRITC_to_Cy5_pairs(mRNA,2)] ==...
                                [Cy5_to_TRITC_pairs(TRITC_to_Cy5_pairs(mRNA,2),2),Cy5_to_TRITC_pairs(TRITC_to_Cy5_pairs(mRNA,2),1)]
                            if dist_TRITC_to_Cy5(mRNA) < 6
                                TRITC_pair_boolean(mRNA) = 1;
                                Cy5_pair_boolean(TRITC_to_Cy5_pairs(mRNA,2),1) = 1;
                            end
                        end
                    end
                else
                    % If the number of detections are equal, or there are fewer Cy5
                    % than TRITC
                    for mRNA = 1:length(Cy5_knn_to_TRITC(:,1))
                        if [Cy5_to_TRITC_pairs(mRNA,1),Cy5_to_TRITC_pairs(mRNA,2)] ==...
                                [TRITC_to_Cy5_pairs(Cy5_to_TRITC_pairs(mRNA,2),2),TRITC_to_Cy5_pairs(Cy5_to_TRITC_pairs(mRNA,2),1)]
                            if dist_Cy5_to_TRITC(mRNA) < 6
                                Cy5_pair_boolean(mRNA,1) = 1; %Mark the index of the paired spot with a 1
                                TRITC_pair_boolean(Cy5_to_TRITC_pairs(mRNA,2),1) = 1; %Mark the index of the partner spot with a 1
                            end
                        end
                    end
                end
                
                TRITC_pair_boolean = logical(TRITC_pair_boolean);
                Cy5_pair_boolean = logical(Cy5_pair_boolean);
                
                %Put all the information into the array, for each cell
                %Put spot info into one cell
                FISH_data_individual_cells(ROI).paired_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(TRITC_pair_boolean,1:6);
                FISH_data_individual_cells(ROI).paired_TRITC(:,7) = dist_TRITC_to_Cy5(TRITC_pair_boolean,1);
                %Spot info for solo TRITC mRNA
                FISH_data_individual_cells(ROI).solo_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(~TRITC_pair_boolean,1:6);
                FISH_data_individual_cells(ROI).solo_TRITC(:,7) = dist_TRITC_to_Cy5(~TRITC_pair_boolean,1);
                
                %Distances for spots that were paired, from Cy5 to nearest TRITC
                FISH_data_individual_cells(ROI).paired_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(Cy5_pair_boolean,1:6);
                FISH_data_individual_cells(ROI).paired_Cy5(:,7) = dist_Cy5_to_TRITC(Cy5_pair_boolean,1);
                %Spot info for solo Cy5 mRNA
                FISH_data_individual_cells(ROI).solo_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(~Cy5_pair_boolean,1:6);
                FISH_data_individual_cells(ROI).solo_Cy5(:,7) = dist_Cy5_to_TRITC(~Cy5_pair_boolean,1);
                
                %Total number paired
                number_paired = length(FISH_data_individual_cells(ROI).paired_TRITC(:,1));
                %Total green alone spots
                green_alone = length(FISH_data_individual_cells(ROI).solo_TRITC(:,1));
                %Total red spots alone
                red_alone = length(FISH_data_individual_cells(ROI).solo_Cy5(:,1));
                %Counting stats for easy export
                FISH_data_individual_cells(ROI).counting_stats = [number_paired, green_alone, red_alone,(number_paired+green_alone),...
                    (number_paired+red_alone), (number_paired + green_alone + red_alone), number_paired/(number_paired + green_alone + red_alone)];
                
                %For the Cy5 spot matrix, since Cy5 is the common probe for all
                %mRNA species, it will be the anchor. Loop through to find the
                %matching spot in the TRITC channel and pair them up for
                %side-by-side comparison.
                
                for Cy5_spot = 1:size(FISH_data_individual_cells(ROI).paired_Cy5,1)
                    FISH_data_individual_cells(ROI).paired_Cy5(Cy5_spot,8:13) = ...
                        FISH_data_individual_cells(ROI).TRITC_spots_threshold(Cy5_to_TRITC_pairs(Cy5_spot,2),1:6);
                end
            else
                FISH_data_individual_cells(ROI).paired_TRITC = [];
                FISH_data_individual_cells(ROI).paired_Cy5 = [];
                FISH_data_individual_cells(ROI).solo_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:6);
                FISH_data_individual_cells(ROI).solo_TRITC(:,7) = NaN;
                FISH_data_individual_cells(ROI).solo_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:6);
                FISH_data_individual_cells(ROI).solo_Cy5(:,7) = NaN;
                
                %Total number paired
                number_paired = 0;
                %Total green alone spots
                green_alone = length(FISH_data_individual_cells(ROI).solo_TRITC(:,1));
                %Total red spots alone
                red_alone = length(FISH_data_individual_cells(ROI).solo_Cy5(:,1));
                %Counting stats for easy export
                FISH_data_individual_cells(ROI).counting_stats = [number_paired, green_alone, red_alone,(number_paired+green_alone),...
                    (number_paired+red_alone), (number_paired + green_alone + red_alone), number_paired/(number_paired + green_alone + red_alone)];  
            end
            
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%% PLOT SPOTS OVER CELLS TO INSPECT THRESHOLDING %%%%%%%%%%%%%%%%
        
        figure('position',[20 100 1600 1200])
        fig1 = subplot(1,2,1);
        imshow(original_image_TRITC,[0 max(max(original_image_TRITC))/2]);
        title('TRITC image with paired and single detections');
        hold on
        for ROI = 1:FISH_Outlines_Number
            if FISH_data_individual_cells(ROI).counting_stats(1) > 0
                scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),70,'go');
                scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),70,'ro');
            end
            if FISH_data_individual_cells(ROI).counting_stats(2) > 0
                scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),70,'g+');
            end
            if FISH_data_individual_cells(ROI).counting_stats(3) > 0
                scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),70,'r+');
            end
        end
        Outline_masks_all(:,:) = sum(FISH_Outlines_masks,3)>0;
        Cell_outlines = bwboundaries(Outline_masks_all(:,:),8);
        Number_of_cell_boundaries = size(Cell_outlines, 1);
        
        for point = 1 : Number_of_cell_boundaries
            CellBoundary = Cell_outlines{point};
            plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
        end
        hold off;
        impixelinfo;
        
        fig2 = subplot(1,2,2);
        imshow(original_image_Cy5,[0 max(max(original_image_Cy5))/2]);
        title('Cy5 image with paired and single detections');
        hold on
        for ROI = 1:FISH_Outlines_Number
            if FISH_data_individual_cells(ROI).counting_stats(1) > 0
                scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),70,'go');
                scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),70,'ro');
            end
            if FISH_data_individual_cells(ROI).counting_stats(2) > 0
                scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),70,'g+');
            end
            if FISH_data_individual_cells(ROI).counting_stats(3) > 0
                scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),70,'r+');
            end
        end
        for point = 1 : Number_of_cell_boundaries
            CellBoundary = Cell_outlines{point};
            plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
        end
        hold off;
        impixelinfo;
        linkaxes([fig1,fig2],'xy');
        
        message = sprintf('Inspect new curves to identify Signal intensity\n and Signal-to-Noise thresholds for both colors');
        uiwait(msgbox(message));
        
        refine_prompt = {'If these thresholds look acceptable, enter 0. Otherwise enter 1 and continue to refine your thresholding.'};
        refine_prompt_name = 'Wait for user';
        numlines = 1;
        refine_defaults = {'1'};
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        repeat_answer = inputdlg(refine_prompt,refine_prompt_name,numlines,refine_defaults);
        close;
        
        refine_threshold = str2num(repeat_answer{1});
        
        
        seg_prompt = {'TRITC Signal threshold (x10^3)','Cy5 Signal threshold(x10^3)', 'TRITC SNR threshold','Cy5 SNR threshold'};
        seg_prompt_name = 'Wait for user';
        numlines = 1;
        seg_defaults = {seg_answer{1},seg_answer{2},seg_answer{3},seg_answer{4}};
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        seg_answer = inputdlg(seg_prompt,seg_prompt_name,numlines,seg_defaults);
        close;
        
        Threshold_params = [str2num(seg_answer{1})*(10^3),str2num(seg_answer{2})*(10^3),str2num(seg_answer{3}),str2num(seg_answer{4})];
        
        
        TRITC_boolean = and(Localization_matrix_TRITC(:,3) > Threshold_params(1), Localization_matrix_TRITC(:,6) > Threshold_params(3));
        Thresholded_TRITC_matrix = Localization_matrix_TRITC(TRITC_boolean,1:6);
        Cy5_boolean = and(Localization_matrix_Cy5(:,3) > Threshold_params(2), Localization_matrix_Cy5(:,6) > Threshold_params(4));
        Thresholded_Cy5_matrix = Localization_matrix_Cy5(Cy5_boolean,1:6);
        
        if Save_output == 1
            save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
                'thresh_params_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Threshold_params');
        end
        
        
        close all;
    end
    
else
    load([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
        'thresh_params_rep_',Replicate_name,'_pos_',num2str(Position_number)]);
    
    TRITC_boolean = and(Localization_matrix_TRITC(:,3) > Threshold_params(1), Localization_matrix_TRITC(:,6) > Threshold_params(3));
    Thresholded_TRITC_matrix = Localization_matrix_TRITC(TRITC_boolean,1:6);
    Cy5_boolean = and(Localization_matrix_Cy5(:,3) > Threshold_params(2), Localization_matrix_Cy5(:,6) > Threshold_params(4));
    Thresholded_Cy5_matrix = Localization_matrix_Cy5(Cy5_boolean,1:6);
    
    TRITC_mRNA_signal_threshold = [];
    TRITC_signal_raw = sort(Localization_matrix_TRITC(:,3));
    increment = round((TRITC_signal_raw(end)-TRITC_signal_raw(1))/number_to_test);
    iter = 1;
    for threshold = round(TRITC_signal_raw(1)):increment:round(TRITC_signal_raw(end))
        boolean_threshold_applied = TRITC_signal_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        TRITC_mRNA_signal_threshold(iter,1) = threshold;
        TRITC_mRNA_signal_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    Cy5_mRNA_signal_threshold = [];
    Cy5_signal_raw = sort(Localization_matrix_Cy5(:,3));
    increment = round((Cy5_signal_raw(end)-Cy5_signal_raw(1))/number_to_test);
    iter = 1;
    for threshold = round(Cy5_signal_raw(1)):increment:round(Cy5_signal_raw(end))
        boolean_threshold_applied = Cy5_signal_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        Cy5_mRNA_signal_threshold(iter,1) = threshold;
        Cy5_mRNA_signal_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    %SNR
    TRITC_mRNA_SNR_threshold = [];
    TRITC_SNR_raw = sort(Thresholded_TRITC_matrix(:,6));
    iter = 1;
    for threshold = floor(TRITC_SNR_raw(1)):ceil(TRITC_SNR_raw(end))
        
        boolean_threshold_applied = TRITC_SNR_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        TRITC_mRNA_SNR_threshold(iter,1) = threshold;
        TRITC_mRNA_SNR_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    Cy5_mRNA_SNR_threshold = [];
    Cy5_SNR_raw = sort(Thresholded_Cy5_matrix(:,6));
    iter = 1;
    for threshold = floor(Cy5_SNR_raw(1)):ceil(Cy5_SNR_raw(end))
        
        boolean_threshold_applied = Cy5_SNR_raw > threshold;
        boolean_threshold_applied = boolean_threshold_applied(boolean_threshold_applied ~=0);
        number_above = length(boolean_threshold_applied);
        Cy5_mRNA_SNR_threshold(iter,1) = threshold;
        Cy5_mRNA_SNR_threshold(iter,2) = number_above;
        iter = iter+1;
    end
    
    %%% PLOT QC statistics
    figure('position',[20 100 1600 1200])
    subplot(2,3,1)
    plot(TRITC_mRNA_signal_threshold(:,1),TRITC_mRNA_signal_threshold(:,2))
    maxYValue = ylim;
    hold on;
    title('TRITC mRNA signal after threshold')
    plot([Threshold_params(1), Threshold_params(1)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
    % Place a text label on the bar chart showing the threshold.
    annotationText = sprintf('Thresholded at %d x10^3 gray levels', Threshold_params(1));
    text(double(Threshold_params(1) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
    
    
    subplot(2,3,2)
    plot(TRITC_mRNA_SNR_threshold(:,1),TRITC_mRNA_SNR_threshold(:,2))
    title('TRITC mRNA SNR after threshold')
    hold on;
    maxYValue = ylim;
    plot([Threshold_params(3), Threshold_params(3)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
    annotationText = sprintf('Thresholded at %d ', Threshold_params(3));
    text(double(Threshold_params(3) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
    
    
    subplot(2,3,3)
    semilogx(Localization_matrix_TRITC(:,3), Localization_matrix_TRITC(:,6),'o')
    title('TRITC Spot intensity vs SNR, after threshold')
    hold on;
    maxXValue = xlim;
    plot([Threshold_params(1), Threshold_params(1)], [maxYValue(1),max(Localization_matrix_TRITC(:,6))*1.1], 'Color', 'r');
    plot([maxXValue(1),maxXValue(2)],[Threshold_params(3), Threshold_params(3)], 'Color', 'r');
    hold off;
    
    subplot(2,3,4)
    plot(Cy5_mRNA_signal_threshold(:,1),Cy5_mRNA_signal_threshold(:,2))
    title('Cy5 mRNA count as a function of signal intensity')
    maxYValue = ylim;
    hold on;
    title('Cy5 mRNA signal after threshold')
    plot([Threshold_params(2), Threshold_params(2)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
    % Place a text label on the bar chart showing the threshold.
    annotationText = sprintf('Thresholded at %d x10^3 gray levels', Threshold_params(2));
    text(double(Threshold_params(2) *1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
    
    subplot(2,3,5)
    plot(Cy5_mRNA_SNR_threshold(:,1),Cy5_mRNA_SNR_threshold(:,2))
    title('Cy5 mRNA count as a function of Signal to Noise ratio')
    maxYValue = ylim;
    hold on;
    title('Cy5 mRNA SNR after threshold')
    plot([Threshold_params(4), Threshold_params(4)], [maxYValue(1),maxYValue(2)], 'Color', 'r');
    % Place a text label on the bar chart showing the threshold.
    annotationText = sprintf('Thresholded at %d ', Threshold_params(4));
    text(double(Threshold_params(4)*1.1), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
    
    subplot(2,3,6)
    semilogx(Localization_matrix_Cy5(:,3), Localization_matrix_Cy5(:,6),'o')
    title('Cy5 Spot intensity vs SNR, after threshold')
    hold on;
    maxXValue = xlim;
    plot([Threshold_params(2), Threshold_params(2)], [maxYValue(1),max(Localization_matrix_Cy5(:,6))*1.1], 'Color', 'r');
    plot([maxXValue(1),maxXValue(2)], [Threshold_params(4), Threshold_params(4)],'Color', 'r');
    hold off;
    
    for ROI = 1:FISH_Outlines_Number
        [in,on] = inpolygon(Thresholded_TRITC_matrix(:,1),Thresholded_TRITC_matrix(:,2),...
            FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
        FISH_data_individual_cells(ROI).TRITC_spots_threshold =Thresholded_TRITC_matrix(in,1:6);
        
        [in,on] = inpolygon(Thresholded_Cy5_matrix(:,1),Thresholded_Cy5_matrix(:,2),...
            FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
        FISH_data_individual_cells(ROI).Cy5_spots_threshold =Thresholded_Cy5_matrix(in,1:6);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%% PAIR UP SPOTS USING KNN ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%
    for ROI = 1:FISH_Outlines_Number
        
        %By constraining this to per-cell analysis, can limit the number and
        %distance of potential neighbors and limit edge problems.
        if and(~isempty(FISH_data_individual_cells(ROI).Cy5_spots_threshold),~isempty(FISH_data_individual_cells(ROI).TRITC_spots_threshold))
            %By constraining this to per-cell analysis, can limit the number and
            %distance of potential neighbors and limit edge problems.
            
            [TRITC_knn_to_Cy5, dist_TRITC_to_Cy5] = knnsearch(FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:2),...
                FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:2),'K',1);
            [Cy5_knn_to_TRITC, dist_Cy5_to_TRITC] = knnsearch(FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:2),...
                FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:2),'K',1);
            
            %%% Pairing up spot first by checking for "unambiguous pairings", meaning
            %%% that if A thinks it's nearest neighbor is B, and if B thinks its
            %%% nearest neighbor is A, this is unambiguous. Ovbiously this can only be
            %%% true for the smaller of the two lists.
            TRITC_to_Cy5_pairs = zeros(length(TRITC_knn_to_Cy5(:,1)),2);
            for row = 1:length(TRITC_knn_to_Cy5(:,1))
                TRITC_to_Cy5_pairs(row,1:2) = [row, TRITC_knn_to_Cy5(row)];
            end
            
            Cy5_to_TRITC_pairs = zeros(length(Cy5_knn_to_TRITC(:,1)),2);
            for row = 1:length(Cy5_knn_to_TRITC(:,1))
                Cy5_to_TRITC_pairs(row,1:2) = [row, Cy5_knn_to_TRITC(row)];
            end
            
            Cy5_pair_boolean = zeros(length(Cy5_knn_to_TRITC(:,1)),1);
            TRITC_pair_boolean = zeros(length(TRITC_knn_to_Cy5(:,1)),1);
            
            if length(Cy5_knn_to_TRITC(:,1)) > length(TRITC_knn_to_Cy5(:,1))
                % If there are fewer TRITC detections that Cy5
                for mRNA = 1:length(TRITC_knn_to_Cy5(:,1))
                    if [TRITC_to_Cy5_pairs(mRNA,1),TRITC_to_Cy5_pairs(mRNA,2)] ==...
                            [Cy5_to_TRITC_pairs(TRITC_to_Cy5_pairs(mRNA,2),2),Cy5_to_TRITC_pairs(TRITC_to_Cy5_pairs(mRNA,2),1)]
                        if dist_Cy5_to_TRITC(mRNA,1) < 6
                            TRITC_pair_boolean(mRNA) = 1;
                            Cy5_pair_boolean(TRITC_to_Cy5_pairs(mRNA,2)) = 1;
                        end
                    end
                end
            else
                % If the number of detections are equal, or there are fewer Cy5
                % than TRITC
                for mRNA = 1:length(Cy5_knn_to_TRITC(:,1))
                    if [Cy5_to_TRITC_pairs(mRNA,1),Cy5_to_TRITC_pairs(mRNA,2)] ==...
                            [TRITC_to_Cy5_pairs(Cy5_to_TRITC_pairs(mRNA,2),2),TRITC_to_Cy5_pairs(Cy5_to_TRITC_pairs(mRNA,2),1)]
                        if dist_TRITC_to_Cy5(mRNA) < 6
                            Cy5_pair_boolean(mRNA,1) = 1; %Mark the index of the paired spot with a 1
                            TRITC_pair_boolean(Cy5_to_TRITC_pairs(mRNA,2),1) = 1; %Mark the index of the partner spot with a 1
                        end
                    end
                end
            end
            
            TRITC_pair_boolean = logical(TRITC_pair_boolean);
            Cy5_pair_boolean = logical(Cy5_pair_boolean);
            
            %Put all the information into the array, for each cell
            %Put spot info into one cell
            FISH_data_individual_cells(ROI).paired_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(TRITC_pair_boolean,1:6);
            FISH_data_individual_cells(ROI).paired_TRITC(:,7) = dist_TRITC_to_Cy5(TRITC_pair_boolean,1);
            %Spot info for solo TRITC mRNA
            FISH_data_individual_cells(ROI).solo_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(~TRITC_pair_boolean,1:6);
            FISH_data_individual_cells(ROI).solo_TRITC(:,7) = dist_TRITC_to_Cy5(~TRITC_pair_boolean,1);
            
            %Distances for spots that were paired, from Cy5 to nearest TRITC
            FISH_data_individual_cells(ROI).paired_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(Cy5_pair_boolean,1:6);
            FISH_data_individual_cells(ROI).paired_Cy5(:,7) = dist_Cy5_to_TRITC(Cy5_pair_boolean,1);
            %Spot info for solo Cy5 mRNA
            FISH_data_individual_cells(ROI).solo_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(~Cy5_pair_boolean,1:6);
            FISH_data_individual_cells(ROI).solo_Cy5(:,7) = dist_Cy5_to_TRITC(~Cy5_pair_boolean,1);
            
            %Total number paired
            number_paired = length(FISH_data_individual_cells(ROI).paired_TRITC(:,1));
            %Total green alone spots
            green_alone = length(FISH_data_individual_cells(ROI).solo_TRITC(:,1));
            %Total red spots alone
            red_alone = length(FISH_data_individual_cells(ROI).solo_Cy5(:,1));
            %Counting stats for easy export
            FISH_data_individual_cells(ROI).counting_stats = [number_paired, green_alone, red_alone,(number_paired+green_alone),...
                (number_paired+red_alone), (number_paired + green_alone + red_alone), number_paired/(number_paired + green_alone + red_alone)];
            
            %For the Cy5 spot matrix, since Cy5 is the common probe for all
            %mRNA species, it will be the anchor. Loop through to find the
            %matching spot in the TRITC channel and pair them up for
            %side-by-side comparison.
            
            for Cy5_spot = 1:size(FISH_data_individual_cells(ROI).paired_Cy5,1)
                FISH_data_individual_cells(ROI).paired_Cy5(Cy5_spot,8:13) = ...
                    FISH_data_individual_cells(ROI).TRITC_spots_threshold(Cy5_to_TRITC_pairs(Cy5_spot,2),1:6);
            end
        else
            FISH_data_individual_cells(ROI).paired_TRITC = [];
            FISH_data_individual_cells(ROI).paired_Cy5 = [];
            FISH_data_individual_cells(ROI).solo_TRITC = FISH_data_individual_cells(ROI).TRITC_spots_threshold(:,1:6);
            FISH_data_individual_cells(ROI).solo_TRITC(:,7) = NaN;
            FISH_data_individual_cells(ROI).solo_Cy5 = FISH_data_individual_cells(ROI).Cy5_spots_threshold(:,1:6);
            FISH_data_individual_cells(ROI).solo_Cy5(:,7) = NaN;
            
            %Total number paired
            number_paired = 0;
            %Total green alone spots
            green_alone = length(FISH_data_individual_cells(ROI).solo_TRITC(:,1));
            %Total red spots alone
            red_alone = length(FISH_data_individual_cells(ROI).solo_Cy5(:,1));
            %Counting stats for easy export
            FISH_data_individual_cells(ROI).counting_stats = [number_paired, green_alone, red_alone,(number_paired+green_alone),...
                (number_paired+red_alone), (number_paired + green_alone + red_alone), number_paired/(number_paired + green_alone + red_alone)];
            
        end
        
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% PLOT SPOTS OVER CELLS TO INSPECT THRESHOLDING %%%%%%%%%%%%%%%%
    
    figure('position',[20 100 1600 1200])
    fig1 = subplot(1,2,1);
    imshow(original_image_TRITC,[0 max(max(original_image_TRITC))/2]);
    title('TRITC image with paired and single detections');
    hold on
    for ROI = 1:FISH_Outlines_Number
        if FISH_data_individual_cells(ROI).counting_stats(1) > 0
            scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),70,'go');
            scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),70,'ro');
        end
        if FISH_data_individual_cells(ROI).counting_stats(2) > 0
            scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),70,'g+');
        end
        if FISH_data_individual_cells(ROI).counting_stats(3) > 0
            scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),70,'r+');
        end
    end
    Outline_masks_all(:,:) = sum(FISH_Outlines_masks,3)>0;
    Cell_outlines = bwboundaries(Outline_masks_all(:,:),8);
    Number_of_cell_boundaries = size(Cell_outlines, 1);
    for point = 1 : Number_of_cell_boundaries
        CellBoundary = Cell_outlines{point};
        plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
    end
    hold off;
    impixelinfo;
    
    fig2 = subplot(1,2,2);
    imshow(original_image_Cy5,[0 max(max(original_image_Cy5))/2]);
    title('Cy5 image with paired and single detections');
    hold on
    for ROI = 1:FISH_Outlines_Number
        if FISH_data_individual_cells(ROI).counting_stats(1) > 0
            scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),70,'go');
            scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),70,'ro');
        end
        if FISH_data_individual_cells(ROI).counting_stats(2) > 0
            scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),70,'g+');
        end
        if FISH_data_individual_cells(ROI).counting_stats(3) > 0
            scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),70,'r+');
        end
    end
    for point = 1 : Number_of_cell_boundaries
        CellBoundary = Cell_outlines{point};
        plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
    end
    hold off;
    impixelinfo;
    linkaxes([fig1,fig2],'xy');
    
    message = sprintf('Inspect images and curves to check whether the Signal intensity threshold \n and Signal-to-Noise thresholds are appropriate.');
    uiwait(msgbox(message));
    check_threshold = 1;
    check_prompt = {'If these thresholds look acceptable, enter 1. Otherwise enter 0 and start over.'};
    refine_prompt_name = 'Wait for user';
    numlines = 1;
    refine_defaults = {'1'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    repeat_answer = inputdlg(check_prompt,refine_prompt_name,numlines,refine_defaults);
    close;
    
    check_threshold = str2num(repeat_answer{1});
    if check_threshold == 0
        error('You have found the thresholds you loaded to be incorrect. Please start the script from the begining to determine the correct thresholds for the data.')
    end
    
    
end

%% QUANTITATE THE NUCLEAR RNA SEPARATELY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ROI = 1:FISH_Outlines_Number
    if ~isempty(FISH_data_individual_cells(ROI).nucleus_outline)
        %Cy5 paired spots
        if ~isempty(FISH_data_individual_cells(ROI).paired_Cy5)
            [in,on] = inpolygon(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),...
                FISH_data_individual_cells(ROI).nucleus_outline(:,2),FISH_data_individual_cells(ROI).nucleus_outline(:,1));
            FISH_data_individual_cells(ROI).Nuclear_Cy5_paired = FISH_data_individual_cells(ROI).paired_Cy5(in,:);
        else
            FISH_data_individual_cells(ROI).Nuclear_Cy5_paired = [];
        end
        
        %TRITC paired spots
        if ~isempty(FISH_data_individual_cells(ROI).paired_TRITC)
            [in,on] = inpolygon(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),...
                FISH_data_individual_cells(ROI).nucleus_outline(:,2),FISH_data_individual_cells(ROI).nucleus_outline(:,1));
            FISH_data_individual_cells(ROI).Nuclear_TRITC_paired = FISH_data_individual_cells(ROI).paired_TRITC(in,:);
        else
            FISH_data_individual_cells(ROI).Nuclear_TRITC_paired = [];
        end
        % Cy5 solo spots
        if ~isempty(FISH_data_individual_cells(ROI).solo_Cy5)
            [in,on] = inpolygon(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),...
                FISH_data_individual_cells(ROI).nucleus_outline(:,2),FISH_data_individual_cells(ROI).nucleus_outline(:,1));
            FISH_data_individual_cells(ROI).Nuclear_solo_Cy5 = FISH_data_individual_cells(ROI).solo_Cy5(in,:);
        else
            FISH_data_individual_cells(ROI).Nuclear_solo_Cy5 = [];
        end
        % TRITC solo spots
        if ~isempty(FISH_data_individual_cells(ROI).solo_Cy5)
            [in,on] = inpolygon(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),...
                FISH_data_individual_cells(ROI).nucleus_outline(:,2),FISH_data_individual_cells(ROI).nucleus_outline(:,1));
            FISH_data_individual_cells(ROI).Nuclear_solo_TRITC = FISH_data_individual_cells(ROI).solo_TRITC(in,:);
        else
            FISH_data_individual_cells(ROI).Nuclear_solo_TRITC = [];
        end
    else
        FISH_data_individual_cells(ROI).Nuclear_Cy5_paired = [];
        FISH_data_individual_cells(ROI).Nuclear_TRITC_paired = [];
        FISH_data_individual_cells(ROI).Nuclear_solo_Cy5 = [];
        FISH_data_individual_cells(ROI).Nuclear_solo_TRITC = [];
    end
    FISH_data_individual_cells(ROI).nuc_counting_stats = [size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1),...
        size(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC,1), size(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5,1),(size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1) + size(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC,1)),...
        (size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1) + size(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5,1)),(size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1) + size(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC,1) +size(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5,1)),...
        size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1)/(size(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired,1) + size(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC,1) +size(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5,1))];
    All_counting_stats_nuclear(ROI,:) = FISH_data_individual_cells(ROI).nuc_counting_stats;
end
        

%% FIND RADIAL DISTANCE AND ANGLE BETWEEN SPOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% for ROI = 1:FISH_Outlines_Number
%     [spots_inside_cell_index,spots_edge_cell_index] = inpolygon(Thresholded_Cy5_matrix(:,1),Thresholded_Cy5_matrix(:,2),...
%         FISH_data_individual_cells(ROI).outline(:,1),FISH_data_individual_cells(ROI).outline(:,2));
%     FISH_data_individual_cells(ROI).Cy5_spots_threshold =Thresholded_Cy5_matrix(spots_inside_cell_index,1:6);
%     delta_x = FISH_data_individual_cells(ROI).paired_Cy5(:,8) - FISH_data_individual_cells(ROI).paired_Cy5(:,1);
%     delta_y = FISH_data_individual_cells(ROI).paired_Cy5(:,9) - FISH_data_individual_cells(ROI).paired_Cy5(:,1);
%     [theta, r] = cart2pol(delta_x,delta_y);
%     FISH_data_individual_cells(ROI).paired_Cy5(:,14) = theta(:,1);
% end


%% SAVE THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

All_paired_spot_matrix = [];
All_solo_TRITC_spot_matrix = [];
All_solo_Cy5_spot_matrix = [];
All_counting_stats = [];

for ROI = 1:FISH_Outlines_Number
    if ~isempty(FISH_data_individual_cells(ROI).paired_Cy5)
        All_paired_spot_matrix = vertcat(All_paired_spot_matrix, FISH_data_individual_cells(ROI).paired_Cy5);
    end
    All_solo_TRITC_spot_matrix = vertcat(All_solo_TRITC_spot_matrix,FISH_data_individual_cells(ROI).solo_TRITC);
    All_solo_Cy5_spot_matrix = vertcat(All_solo_Cy5_spot_matrix, FISH_data_individual_cells(ROI).solo_Cy5);
    All_counting_stats = vertcat(All_counting_stats,FISH_data_individual_cells(ROI).counting_stats);
end

if Save_output == 1
    save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
        'thresh_params_rep_',Replicate_name,'_pos_',num2str(Position_number)],'Threshold_params');
    save([output_path,Genotype,'/',condition,'/',Replicate_name,'/',num2str(Position_number),'/',...
        'Quant_FISH_rep_',Replicate_name,'_pos_',num2str(Position_number)],'FISH_data_individual_cells',...
        'All_paired_spot_matrix','All_solo_TRITC_spot_matrix','All_solo_Cy5_spot_matrix','All_counting_stats');
    save([quantified_spot_path,'Quant_FISH_rep_',Replicate_name,'_pos_',num2str(Position_number)],'FISH_data_individual_cells',...
        'All_paired_spot_matrix','All_solo_TRITC_spot_matrix','All_solo_Cy5_spot_matrix','All_counting_stats');
end


%% PLOT "FINAL" SPOT PAIRS
% PLOT final spot position/pairs
figure('position',[20 100 1600 1200])
fig1 = subplot(1,2,1);
imshow(original_image_TRITC,[]);
title('TRITC image with paired and single detections');
hold on
for ROI = 1:FISH_Outlines_Number
    if FISH_data_individual_cells(ROI).counting_stats(1) > 0
        scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),40,'go');
        scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),50,'ro');
    end
    if FISH_data_individual_cells(ROI).counting_stats(2) > 0
        scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),60,'g+');
    end
    if FISH_data_individual_cells(ROI).counting_stats(3) > 0
        scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),70,'r+');
    end
end
for point = 1 : Number_of_cell_boundaries
    CellBoundary = Cell_outlines{point};
    plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
end
for point = 1 : Number_of_annotated_nuclei
    NucBoundary = DAPI_outlines_annotated{point};
    plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
end
    for cell = 1:region
        text(FISH_data_individual_cells(cell).centroid(1), FISH_data_individual_cells(cell).centroid(2), num2str(FISH_data_individual_cells(cell).id), 'FontSize', 12, 'FontWeight', 'Bold','Color','w');
    end
    hold off
impixelinfo;
hold off;

fig2 = subplot(1,2,2);
imshow(original_image_Cy5,[]);
title('Cy5 image with paired and single detections');
hold on
for ROI = 1:FISH_Outlines_Number
    if FISH_data_individual_cells(ROI).counting_stats(1) > 0
        scatter(FISH_data_individual_cells(ROI).paired_TRITC(:,1),FISH_data_individual_cells(ROI).paired_TRITC(:,2),25,'go');
        scatter(FISH_data_individual_cells(ROI).paired_Cy5(:,1),FISH_data_individual_cells(ROI).paired_Cy5(:,2),25,'ro');
    end
    if FISH_data_individual_cells(ROI).counting_stats(2) > 0
        scatter(FISH_data_individual_cells(ROI).solo_TRITC(:,1),FISH_data_individual_cells(ROI).solo_TRITC(:,2),25,'g+');
    end
    if FISH_data_individual_cells(ROI).counting_stats(3) > 0
        scatter(FISH_data_individual_cells(ROI).solo_Cy5(:,1),FISH_data_individual_cells(ROI).solo_Cy5(:,2),25,'r+');
    end
end
for point = 1 : Number_of_cell_boundaries
    CellBoundary = Cell_outlines{point};
    plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
end
for point = 1 : Number_of_annotated_nuclei
    NucBoundary = DAPI_outlines_annotated{point};
    plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
end
    for cell = 1:region
        text(FISH_data_individual_cells(cell).centroid(1), FISH_data_individual_cells(cell).centroid(2), num2str(FISH_data_individual_cells(cell).id), 'FontSize', 12, 'FontWeight', 'Bold','Color','w');
    end
    hold off
impixelinfo;
hold off;
linkaxes([fig1,fig2],'xy');


% PLOT only nuclear spot position/pairs
figure('position',[20 100 1600 1200])
fig1 = subplot(1,2,1);
imshow(original_image_TRITC,[]);
title('TRITC image with paired and single detections');
hold on
for ROI = 1:FISH_Outlines_Number
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_TRITC_paired)
        scatter(FISH_data_individual_cells(ROI).Nuclear_TRITC_paired(:,1),FISH_data_individual_cells(ROI).Nuclear_TRITC_paired(:,2),40,'go');
        scatter(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired(:,1),FISH_data_individual_cells(ROI).Nuclear_Cy5_paired(:,2),50,'ro');
    end
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC)
        scatter(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC(:,1),FISH_data_individual_cells(ROI).Nuclear_solo_TRITC(:,2),60,'g+');
    end
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5)
        scatter(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5(:,1),FISH_data_individual_cells(ROI).Nuclear_solo_Cy5(:,2),70,'r+');
    end
end
for point = 1 : Number_of_cell_boundaries
    CellBoundary = Cell_outlines{point};
    plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
end
for point = 1 : Number_of_annotated_nuclei
    NucBoundary = DAPI_outlines_annotated{point};
    plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
end
for cell = 1:region
    text(FISH_data_individual_cells(cell).centroid(1), FISH_data_individual_cells(cell).centroid(2), num2str(FISH_data_individual_cells(cell).id), 'FontSize', 12, 'FontWeight', 'Bold','Color','w');
end
impixelinfo;
hold off;

fig2 = subplot(1,2,2);
imshow(original_image_Cy5,[]);
title('Cy5 image with paired and single detections');
hold on
for ROI = 1:FISH_Outlines_Number
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_TRITC_paired)
        scatter(FISH_data_individual_cells(ROI).Nuclear_TRITC_paired(:,1),FISH_data_individual_cells(ROI).Nuclear_TRITC_paired(:,2),25,'go');
        scatter(FISH_data_individual_cells(ROI).Nuclear_Cy5_paired(:,1),FISH_data_individual_cells(ROI).Nuclear_Cy5_paired(:,2),25,'ro');
    end
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC)
        scatter(FISH_data_individual_cells(ROI).Nuclear_solo_TRITC(:,1),FISH_data_individual_cells(ROI).Nuclear_solo_TRITC(:,2),25,'g+');
    end
    if ~isempty(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5)
        scatter(FISH_data_individual_cells(ROI).Nuclear_solo_Cy5(:,1),FISH_data_individual_cells(ROI).Nuclear_solo_Cy5(:,2),25,'r+');
    end
end
for point = 1 : Number_of_cell_boundaries
    CellBoundary = Cell_outlines{point};
    plot(CellBoundary(:,2), CellBoundary(:,1), 'w', 'LineWidth', 1);
end
for point = 1 : Number_of_annotated_nuclei
    NucBoundary = DAPI_outlines_annotated{point};
    plot(NucBoundary(:,2), NucBoundary(:,1), 'b', 'LineWidth', 2);
end
impixelinfo;
hold off;
linkaxes([fig1,fig2],'xy');

% % PLOT HISTOGRAMS
% 
% figure('position',[20 100 1600 1200])
% Plot the pairwise distances of "true pairs" and other spots
% spot_distances_pairs = histogram(All_paired_spot_matrix(:,7));
% hold on;
% spot_distances_pairs.Normalization = 'pdf';
% spot_distances_pairs.BinWidth = 0.5;
% spot_distances_pairs.FaceColor = [0.8185,0.7327,0.3498];
% spot_distances_TRITC = histogram(All_solo_TRITC_spot_matrix(:,7));
% spot_distances_TRITC.Normalization = 'pdf';
% spot_distances_TRITC.BinWidth = 0.5;
% spot_distances_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
% spot_distances_Cy5 = histogram(All_solo_Cy5_spot_matrix(:,7));
% spot_distances_Cy5.Normalization = 'pdf';
% spot_distances_Cy5.BinWidth = 0.5;
% spot_distances_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% title ('Distance from Cy5 spot to closest neighbor');
% title ('Distance to nearest neighbor','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Distance (px)', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Unpaired TRITC detections','Unpaired Cy5 detections','Box','off','FontSize',14,'Location','North');
% 
% 
% figure('position',[20 100 1600 1200])
% Plot the pairwise distances of "true pairs" and other spots
% spot_distances_pairs = histfit(All_paired_spot_matrix(:,7),50,'Gamma');
% hold on;
% spot_distances_pairs(1).FaceColor = [0.8185,0.7327,0.3498];
% spot_distances_pairs(1).FaceAlpha = 0.6;
% spot_distances_pairs(2).Color = [0.8185,0.7327,0.3498];
% spot_distances_TRITC = histfit(All_solo_TRITC_spot_matrix(:,7),50,'Gamma');
% spot_distances_TRITC(1).FaceColor = [0.1801, 0.7177, 0.6424];
% spot_distances_TRITC(2).Color = [0.1801, 0.7177, 0.6424];
% spot_distances_TRITC(1).FaceAlpha = 0.6;
% spot_distances_Cy5 = histfit(All_solo_Cy5_spot_matrix(:,7),50,'Gamma');
% spot_distances_Cy5(1).FaceColor = [0.208, 0.1663, 0.5292];
% spot_distances_Cy5(2).Color = [0.208, 0.1663, 0.5292];
% spot_distances_Cy5(1).FaceAlpha = 0.6;
% hold off;
% title ('Distance from Cy5 spot to closest neighbor');
% title ('Distance to nearest neighbor','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Distance (px)', 'FontSize',14, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',14, 'FontName', 'Helvetica');
% legend('Paired detections','Fit','Unpaired TRITC detections','Fit','Unpaired Cy5 detections','Fit','Box','off','FontSize',18,'Location','North');
% 
% 
% % Paired VS Unpaired
% figure('position',[20 100 1600 1200])
% subplot(2,1,1)
% Plot the spot intensities paired and un-paired spots with TRITC
% spot_intensity_pairs = histogram(All_paired_spot_matrix(:,10),50);
% normalized_spot_intensity_pairs_tritc = spot_intensity_pairs.Values ./sum(spot_intensity_pairs.Values);
% hold on;
% spot_intensity_pairs.Normalization = 'pdf';
% spot_intensity_pairs.BinWidth = 400;
% spot_intensity_pairs.FaceColor = [0.8185,0.7327,0.3498];
% spot_intensity_TRITC = histogram(All_solo_TRITC_spot_matrix(:,3),50);
% normalized_spot_intensity_solo_tritc = spot_intensity_TRITC.Values ./sum(spot_intensity_TRITC.Values);
% spot_intensity_TRITC.Normalization = 'pdf';
% spot_intensity_TRITC.BinWidth = 2000;
% spot_intensity_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
% hold off;
% title ('Integrated spot intensity for TRITC detections','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Integrated Spot Intensity (A.U.)', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Unpaired detections','Box','off','FontSize',14,'Location','North');
% 
% subplot(2,1,2)
% Plot the spot intensities paired and un-paired spots with Cy5
% spot_intensity_pairs = histogram(All_paired_spot_matrix(:,3),50);
% normalized_spot_intensity_pairs_cy5 = spot_intensity_pairs.Values ./sum(spot_intensity_pairs.Values);
% hold on;
% spot_intensity_pairs.Normalization = 'pdf';
% spot_intensity_pairs.BinWidth = 400;
% spot_intensity_pairs.FaceColor = [0.8185,0.7327,0.3498];
% spot_intensity_Cy5 = histogram(All_solo_Cy5_spot_matrix(:,3),50);
% normalized_spot_intensity_solo_cy5 = spot_intensity_Cy5.Values ./sum(spot_intensity_Cy5.Values);
% spot_intensity_Cy5.Normalization = 'pdf';
% spot_intensity_Cy5.BinWidth = 2000;
% spot_intensity_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% title ('Integrated spot intensity for Cy5 detections','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Integrated Spot Intensity (A.U.)', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Unpaired detections','Box','off','FontSize',14,'Location','North');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure('position',[20 100 1600 1200])
% subplot(2,1,1)
% Plot the spot intensities paired and un-paired spots with TRITC
% spot_intensity_pairs = histfit(datasample(All_paired_spot_matrix(:,3),size(All_solo_TRITC_spot_matrix,1)),50,'Gamma');
% hold on;
% spot_intensity_pairs(1).FaceColor = [0.8185,0.7327,0.3498];
% spot_intensity_pairs(1).FaceAlpha = 0.6;
% spot_intensity_pairs(2).Color = [0.8185,0.7327,0.3498];
% spot_intensity_TRITC = histfit(All_solo_TRITC_spot_matrix(:,3),50,'Gamma');
% spot_intensity_TRITC(1).FaceColor = [0.1801, 0.7177, 0.6424];
% spot_intensity_TRITC(2).Color = [0.1801, 0.7177, 0.6424];
% spot_intensity_TRITC(1).FaceAlpha = 0.6;
% hold off;
% title ('Integrated spot intensity for TRITC detections','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Integrated Spot Intensity (A.U.)', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Fit','Unpaired detections','Fit','Box','off','FontSize',14,'Location','North');
% 
% subplot(2,1,2)
% Plot the spot intensities paired and un-paired spots with Cy5
% spot_intensity_pairs = histfit(datasample(All_paired_spot_matrix(:,3),size(All_solo_Cy5_spot_matrix,1)),50,'Gamma');
% hold on;
% spot_intensity_pairs(1).FaceColor = [0.8185,0.7327,0.3498];
% spot_intensity_pairs(1).FaceAlpha = 0.6;
% spot_intensity_pairs(2).Color = [0.8185,0.7327,0.3498];
% spot_intensity_Cy5 = histfit(All_solo_Cy5_spot_matrix(:,3),50,'Gamma');
% spot_intensity_Cy5(1).FaceColor = [0.208, 0.1663, 0.5292];
% spot_intensity_Cy5(2).Color = [0.208, 0.1663, 0.5292];
% spot_intensity_Cy5(1).FaceAlpha = 0.6;
% hold off;
% title ('Integrated spot intensity for Cy5 detections','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Integrated Spot Intensity (A.U.)', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Fit','Unpaired detections','Fit','Box','off','FontSize',14,'Location','North');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure('position',[20 100 1600 1200])
% Plot number of mRNA detections per cell for
% spot_counts_pairs = histogram(All_counting_stats(:,1));
% hold on;
% spot_counts_pairs.Normalization = 'pdf';
% spot_counts_pairs.BinWidth = 1;
% spot_counts_pairs.FaceColor = [0.8185,0.7327,0.3498];
% spot_counts_TRITC = histogram(All_counting_stats(:,2));
% spot_counts_TRITC.Normalization = 'pdf';
% spot_counts_TRITC.BinWidth = 1;
% spot_counts_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
% spot_counts_Cy5 = histogram(All_counting_stats(:,3));
% spot_counts_Cy5.Normalization = 'pdf';
% spot_counts_Cy5.BinWidth = 1;
% spot_counts_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% title ('Number of mRNA per cell','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Total number of mRNA', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('Paired detections','Unpaired TRITC detections', 'Unpaired Cy5 detections','Box','off','FontSize',14,'Location','North');
% 
% figure('position',[20 100 1600 1200])
% Plot number of mRNA detections per cell for
% spot_counts_total_TRITC = histogram(All_counting_stats(:,4));
% hold on;
% spot_counts_total_TRITC.Normalization = 'pdf';
% spot_counts_total_TRITC.BinWidth = 2;
% spot_counts_total_TRITC.FaceColor = [0.1801, 0.7177, 0.6424];
% spot_counts_total_Cy5 = histogram(All_counting_stats(:,5));
% spot_counts_total_Cy5.Normalization = 'pdf';
% spot_counts_total_Cy5.BinWidth = 2;
% spot_counts_total_Cy5.FaceColor = [0.208, 0.1663, 0.5292];
% hold off;
% title ('Number of mRNA per cell','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Total number of mRNA per cell', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('TRITC Spots','Cy5 spots', 'Box','off','FontSize',14,'Location','North');
% 
% figure('position',[20 100 1600 1200])
% Plot number of mRNA detections per cell for
% spot_counts_total_TRITC = histfit(All_counting_stats(:,4));
% hold on;
% spot_counts_total_TRITC(1).FaceColor = [0.1801, 0.7177, 0.6424];
% spot_counts_total_TRITC(2).Color = [0.1801, 0.7177, 0.6424];
% spot_counts_total_TRITC(1).FaceAlpha = 0.6;
% spot_counts_total_Cy5 = histfit(All_counting_stats(:,5));
% spot_counts_total_Cy5(1).FaceColor = [0.208, 0.1663, 0.5292];
% spot_counts_total_Cy5(2).Color = [0.208, 0.1663, 0.5292];
% spot_counts_total_Cy5(1).FaceAlpha = 0.6;
% hold off;
% title ('Number of mRNA per cell','FontSize',14, 'FontName', 'Helvetica');
% xlabel('Total number of mRNA per cell', 'FontSize',12, 'FontName', 'Helvetica');
% ylabel('Probability', 'FontSize',12, 'FontName', 'Helvetica');
% legend('TRITC Spots','Fit','Cy5 spots', 'Fit','Box','off','FontSize',14,'Location','Northeast');
% 
% figure;
% [a,Vals] = rose(All_paired_spot_matrix(:,14), 20);
% Vals = Vals./sum(Vals);
% MaxLim = max([0.025 1.01*max(Vals)]);
% h = polar(0,MaxLim);
% delete(h);
% set(gca, 'Nextplot','add')
% # draw patches instead of lines: polar(t,r)
% [x,y] = pol2cart(a,Vals);
% h = patch(reshape(x,4,[]), reshape(y,4,[]), [0.1801, 0.7177, 0.6424]);
% title(['Angular displacement between TRITC and Cy5 Channels'], 'FontSize', 10, 'FontName', 'Helvetica');
% xlabel('angle in degrees', 'FontSize',10, 'FontName', 'Helvetica');