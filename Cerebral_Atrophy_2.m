clc
clear
close all

a = input('Please select a normal brain MRI scan: ','s');
b = input('Please select a non-normal brain MRI scan: ','s');
figure(1)
[TotalPercentLoss,MisClassified,HSV_Threshold,numNonBlackPixelsOG] = comparison(a,b);
figure(2)
[percentage_loss_of_section,white_pixels_matrix] = analyze_section(HSV_Threshold,numNonBlackPixelsOG);
figure(3)
color_map = mapping(white_pixels_matrix,HSV_Threshold,numNonBlackPixelsOG);
Diagnosis = diagnose(percentage_loss_of_section);

fprintf('Total Percent Loss: %g\n',TotalPercentLoss);
fprintf('Misclassified Percentage: %g\n',MisClassified);
disp(Diagnosis)

function [percent_loss,missclassified_percent,c,numNonBlackPixelsOG] = comparison(non_filename,etc_filename)

    x = imread(non_filename);
    y = imread(etc_filename);

    windowWidth = 10;
    kernel = ones(windowWidth) / windowWidth ^ 2;
    blurredImage_a = imfilter(x, kernel);
    blurredImage_b = imfilter(y, kernel);
   
    [ssimval,~] = ssim(blurredImage_a,blurredImage_b);

    if ssimval < 0.50
        error('The two images are not the same relative region of the brain! Please select two other images.')
    end

    z = imshowpair(x,y,'falsecolor');
    saveas(z,'paired_image','jpg');
    a = imread('paired_image.jpg');
    new_a = cropwhite(a);
   
    scanned_area = CenterArea(new_a);
    scanned_area_2 = bwareafilt(scanned_area,1);

    stats = regionprops(scanned_area_2, 'Centroid');
    center = struct2array(stats);

    numOfNumberWhitePixelsArea = sum(scanned_area_2(:));
    radius = 1.1 * sqrt(numOfNumberWhitePixelsArea/pi);
   
    c1 = center(1)-radius;
    c2 = center(1)+radius;
    r1 = center(2)-radius;
    r2 = center(2)+radius;

    final_a = a(r1:r2,c1:c2,:);

    c = createMaskHSV(final_a);

    subplot(1,4,1)
    imshow(final_a)
    title('Combined Images')

    e = createMaskBlack(final_a);
    subplot(1,4,2)
    imshow(~e)
    title('Total Brain Area')
   
    subplot(1,4,3)
    imshow(c)
    title('HSV Threshold')
    numWhitePixelsHSV = sum(c(:));

    f = createMaskMisclassified(final_a);
    subplot(1,4,4)
    imshow(f)
    title('Misclassified Areas')
    numWhitePixelsMisclass = sum(f(:));
   
    numNonBlackPixelsOG = sum(~e(:));
   
    percent_loss = (numWhitePixelsHSV/numNonBlackPixelsOG) * 100;
    missclassified_percent = (numWhitePixelsMisclass/numNonBlackPixelsOG) * 100;

end

function [percentage_loss_of_section,white_pixels_matrix] = analyze_section(HSV_Threshold,numNonBlackPixelsOG)
    image = HSV_Threshold;
    [image_height,image_width,~] = size(image);
   
    block_size = 36;
    number_of_blocks_vertically = ceil(image_height/block_size);
    number_of_blocks_horizontally = ceil(image_width/block_size);
   
    index = 1;
    white_pixels_matrix = zeros(number_of_blocks_horizontally,number_of_blocks_vertically);
    rol_index = 0;
    col_index = 0;

    for row = 1:block_size:image_height
        rol_index = rol_index + 1;
        for col = 1:block_size:image_width
   
            col_index = col_index + 1;
   
            row_end = row + block_size -1;
            col_end = col + block_size - 1;
   
            if row_end > image_height
                row_end = image_height;
            end
   
            if col_end > image_width
                col_end = image_width;
            end
   
            temp_tile = image(row:row_end,col:col_end,:);
            subplot(number_of_blocks_vertically,number_of_blocks_horizontally,index);
            imshow(temp_tile);
   
            num_of_white_pixels = max(size(find(temp_tile==1)));
            white_pixels_matrix(rol_index,col_index) = num_of_white_pixels;
   
            index = index + 1;
   
        end
    end

    reshape(white_pixels_matrix,1,[]);
    white_pixels_vector = white_pixels_matrix(white_pixels_matrix~=0)';
    percentage_loss_of_section = (white_pixels_vector/numNonBlackPixelsOG)*100;

end

function color_map = mapping(white_pixels_matrix,HSV_Threshold,numNonBlackPixelsOG)
    image = HSV_Threshold;
    [image_height,image_width,~] = size(image);
   
    block_size = 36;
    number_of_blocks_vertically = ceil(image_height/block_size);
    number_of_blocks_horizontally = ceil(image_width/block_size);

    nonZeroIndices = find(white_pixels_matrix);
    nonZeroValues = (white_pixels_matrix(nonZeroIndices))';
   
    color_matrix = zeros(number_of_blocks_horizontally,number_of_blocks_vertically);
    color_index = 0;
   
    for i = 1:1:number_of_blocks_horizontally
        for j = 1:1:number_of_blocks_vertically
            color_index = color_index + 1;
            color_matrix(i,j) = nonZeroValues(color_index);
        end
    end
    final_color_matrix = (color_matrix/numNonBlackPixelsOG)*100;
    color_map = imagesc(final_color_matrix);
    colorbar
end

function [BW,maskedRGBImage] = createMaskMisclassified(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 08-Nov-2023
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.791;
channel1Max = 0.889;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.249;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.364;
channel3Max = 0.974;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function [BW,maskedRGBImage] = createMaskBlack(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 08-Nov-2023
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.000;
channel1Max = 1.000;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.000;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 0.256;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function [BW,maskedRGBImage] = createMaskHSV(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 08-Nov-2023
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.244;
channel1Max = 0.383;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.422;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.182;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function [y] = cropwhite(x)
    y = x;
    [r,c,~] = size(x);
    for i = 1:1:r
        for j = 1:1:c
            if x(i,j,:) > 253
                y(i,j,:) = 0;
            elseif x(i,j,:) == 0
                continue
            end
        end
    end
end

function [BW,maskedRGBImage] = CenterArea(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 29-Nov-2023
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.000;
channel1Max = 1.000;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.000;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.070;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function disease = diagnose(percentage_loss_of_section)

    x = percentage_loss_of_section;

    frontal_lobe_loss = sum(x(11:20)) + sum(x(22:24));

    basal_ganglia_loss = x(13);

    occipital_lobe_loss = sum(x(2:4));

    temporal_lobe_loss = sum(x(6:20));

    parietal_lobe_loss = sum(x(6:15));

 

    normalized_percentage_fl = frontal_lobe_loss/13;

    normalized_percentage_bg = basal_ganglia_loss/5;

    normalized_percentage_ol = occipital_lobe_loss/3;

    normalized_percentage_tl = temporal_lobe_loss/15;

    normalized_percentage_pl = parietal_lobe_loss/10;

 

    normalized_vector = [normalized_percentage_fl,normalized_percentage_bg,normalized_percentage_ol,normalized_percentage_tl,normalized_percentage_pl];

    sorted_normalized_vector = sort(normalized_vector);

 

    if normalized_percentage_bg == sorted_normalized_vector(5)

        % bg

        disease = "Anatomical Structure(s) with most shrinkage: Basal Ganglia | Suggested Diagnosis: Huntington's Disease";

    elseif (normalized_percentage_ol == sorted_normalized_vector(5) || normalized_percentage_ol == sorted_normalized_vector(4)) && (normalized_percentage_pl == sorted_normalized_vector(5) || normalized_percentage_pl == sorted_normalized_vector(4))

        % ol + pl

        disease = "Anatomical Structure(s) with most shrinkage: Occipital Lobe and Parietal Lobe | Suggested diagnosis: Alzheimer's Disease";

    elseif (normalized_percentage_tl == sorted_normalized_vector(5) || normalized_percentage_tl == sorted_normalized_vector(4)) && (normalized_percentage_ol == sorted_normalized_vector(5) || normalized_percentage_ol == sorted_normalized_vector(4))

        % tl + ol

        disease = "Anatomical Structure(s) with most shrinkage: Temporal Lobe and Occipital Lobe | Suggested diagnosis: Parkinson's Disease";

    elseif (normalized_percentage_fl == sorted_normalized_vector(5) || normalized_percentage_fl == sorted_normalized_vector(4)) && (normalized_percentage_tl == sorted_normalized_vector(5) || normalized_percentage_tl == sorted_normalized_vector(4))

        % fl + tl

        disease = "Anatomical Structure(s) with most shrinkage: Frontal Lobe and Temporal Lobe | Suggested diagnoses: ALS or Frontotemporal Dementia";

    elseif (normalized_percentage_tl == sorted_normalized_vector(5) || normalized_percentage_tl == sorted_normalized_vector(4)) && (normalized_percentage_pl == sorted_normalized_vector(5) || normalized_percentage_pl == sorted_normalized_vector(4))

        % tl + pl

        disease = "Anatomical Structure(s) with most shrinkage: Temporal Lobe and Parietal Lobe | Suggested diagnoses: Alzheimer's, Parkinson's, or Frontotemporal Dementia";

    elseif (normalized_percentage_ol == sorted_normalized_vector(5) || normalized_percentage_ol == sorted_normalized_vector(4)) && (normalized_percentage_fl == sorted_normalized_vector(5) || normalized_percentage_fl == sorted_normalized_vector(4))

        % ol + fl

        disease = "Anatomical Structure(s) with most shrinkage: Occipital Lobe and Frontal Lobe | Suggested diagnoses: Alzheimer's, Parkinson's, ALS, or Frontotemporal Dementia";
 
    elseif (normalized_percentage_fl == sorted_normalized_vector(5) || normalized_percentage_fl == sorted_normalized_vector(4)) && (normalized_percentage_pl == sorted_normalized_vector(5) || normalized_percentage_pl == sorted_normalized_vector(4))

        % fl + pl

        disease = "Anatomical Structure(s) with most shrinkage: Frontal Lobe and Parietal Lobe | Suggested diagnoses: Alzheimer's, ALS, or Frontotemporal Dementia";

    elseif normalized_percentage_bg == sorted_normalized_vector(4)

        % bg + etc

        if normalized_percentage_fl == sorted_normalized_vector(5)

            % fl + bg

            disease = "Anatomical Structure(s) with most shrinkage: Basal Ganglia and Frontal Lobe | Suggested diagnoses: Hunington's, ALS, or Frontotemporal Dementia";

        elseif normalized_percentage_pl == sorted_normalized_vector(5)

            % pl + bg

            disease = "Anatomical Structure(s) with most shrinkage: Basal Ganglia and Parietal Lobe | Suggested diagnoses: Hunington's or Alzheimer's disease";

        elseif normalized_percentage_ol == sorted_normalized_vector(5)

            % ol + bg

            disease = "Anatomical Structure(s) with most shrinkage: Basal Ganglia and Occipital Lobe | Suggested diagnoses: Hunington's, Parkinson's, or Alzheimer's disease";

        elseif normalized_percentage_tl == sorted_normalized_vector(5)

            % tl + bg

            disease = "Anatomical Structure(s) with most shrinkage: Basal Ganglia and Temporal Lobe | Suggested diagnoses: Hunington's, Frontotemporal Dementia, or Parkinson's disease";

        end

    else

        disease = "Inconclusive";

    end

end