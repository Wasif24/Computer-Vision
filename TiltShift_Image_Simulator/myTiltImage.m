%%  Tilt Shift Image Simulator
%   This function will take in an input, and certain customization factors.
%   The result is that the Tilt Shift image is shown.

function myTiltImage(image, doFfactor, bandSize, sigmaFactor)

    %1.a
    img = imread(image);
    imshow(img);

    %1.b
    dims = size(img);  %dimensions of picture
    horizonPos = ginput(1); %row number of the horizon line
    %doFfactor = 50;        %increase value to decrease the depth of field
    doFpos = [horizonPos(2) - dims(1)/doFfactor,...
                horizonPos(2) + dims(1)/doFfactor]; %interval of depth of field
    %bandSize = 20;          %blurring band size
    %sigmaFactor = 200;      %higher value corresponds to a more subtle change in blur
    sigma = bandSize/sigmaFactor;   %sigma for gaussian filter; 

    %1.c
    for i = 0:(doFpos(1)/bandSize - 1)
        img(1:(round(doFpos(1) - i*bandSize)), :,:)...
         =imgaussfilt(img(1:(round(doFpos(1)-i*bandSize)),:,:),1+i*sigma); 
        %gaussian blur img(1st row : depth of field - band size * i)
    end
    for i = 0:(((dims(1)-doFpos(2))/bandSize) - 1)
        img(round(doFpos(2)+i*bandSize:dims(1)), :,:)...
         =imgaussfilt(img(round(doFpos(2)+i*bandSize:dims(1)), :,:),1+i*sigma);
        %guassian blur img(depth of field + bandsize*i : last row)
    end

    %1.d
    img_hsv2 = rgb2hsv(img);
    img_hsv2(:,:,2) = img_hsv2(:,:,2)*2;
    img = hsv2rgb(img_hsv2);

    imshow(img)
end



