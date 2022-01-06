%%  Canny Edge Detector
%   This function will show the edges of an image. The process is completed
%   using the Canny Edge Detector Algorithm.
%   lo and hi specify the threshold for the thresholding step.
%   There are also input options for performing a separable gaussian 
%   convolution, and for completing the Hysteresis step.

function MyCanny(image, lo, hi, doSeparableGauss, doHysteresis)

    img = rgb2gray(imread(image));
    sigma = 1.4;    %value of sigma for gaussian blur
    if(doSeparableGauss == true)
        hSize = 9;
        m = fspecial('gaussian', hSize, sigma);
        mY = m(:, 1);
        mX = mY.';
        img_pad = padarray(img, [hSize hSize], 0, 'both');
        img_pad = conv2(mY, mX, img_pad);
        img = img_pad(sigma*10:size(img_pad, 1)-sigma*10, sigma*10:size(img_pad, 2)-sigma*10);
    else
       img = imgaussfilt(img, sigma); 
    end
    h = fspecial('sobel');
    dx = imfilter(double(img), h', 'conv');
    dy = imfilter(double(img), h, 'conv');
    %pause(1); imshow(dx, []);  %display gradient in x dir 
    %pause(1); imshow(dy, []);  %display gradient in y dir

    img_mag = sqrt(dx.^2 + dy.^2);  %get magnitude
    theta = atan2(dy, dx);          %get orientation
    %pause(1); imshow(img_mag, []); %show magnitude of image
    %pause(1); imshow(img_mag > 15, []); %Threshold image with thick lines

    img_Sup = nonMaxSuppression(img_mag, theta);

    %pause(1); imshow(img_Sup,[]); %show image after Non-Max Suppression
    img_Sup = 255*img_Sup/(max(max(img_Sup)));  %normalize values 


    img_thresh = doThresh(img_Sup, lo, hi, 0, 255);

    %pause(1); imshow(img_thresh,[]); %display final thresholded image
    %pause(1); edge(img); % How the edge detection algorithm for MATLAB looks
    
    %2.3
    if(doHysteresis == true)
        img_thresh = doThresh(img_Sup, lo, hi, 60, 255);
        %pause(1); imshow(img_thresh,[]);% show black white gray image
        img_thresh = drawHysteresis(img_thresh);
        pause(1); imshow(img_thresh,[]); %display final Hystersis image
    else
        pause(1); imshow(img_thresh,[]); %display final thresholded image
    end
end


function img_Sup = nonMaxSuppression(img_mag, theta)
    img_Sup = zeros(size(img_mag));
    for y = 2:size(img_mag, 1)-1        %Non-Max Supression Step
        for x = 2:size(img_mag, 2)-1
            l=0; r=0;
            sel = img_mag(y, x);
            is_max = true;
        
            if(theta(y,x) < 0) %if negative, make positive angle
                theta(y, x) = theta(y,x) + pi;
            end
        
            if(15*pi/8 < theta(y,x) && theta(y,x) < pi/8)           %-%-22.5 < t < 22.5; 
                l = img_mag(y, x-1);                                   % check left and right
                r = img_mag(y, x+1);
            elseif(pi/8 <= theta(y,x) && theta(y,x) < 3*pi/8)       %/%22.5 < t <67.5;
                l = img_mag(y-1, x-1);                                % check bot left and top right
                r = img_mag(y+1, x+1);
            elseif(3*pi/8 <= theta(y,x) && theta(y,x) < 5*pi/8)     %|%67.5 < t <112.5; 
                l = img_mag(y-1, x);                                  % check bot and top
                r = img_mag(y+1, x);
            elseif(5*pi/8 <= theta(y,x) && theta(y,x) < 7*pi/8)     %\%112.5 < t <157.5;
                l = img_mag(y+1, x-1);                                % check top left and bot right
                r = img_mag(y-1, x+1);
            elseif(7*pi/8 <= theta(y,x) && theta(y,x) < 9*pi/8)     %-%157.5 < t <202.5;
                l = img_mag(y, x-1);                                  % check left and right
                r = img_mag(y, x+1);
            elseif(9*pi/8 <= theta(y,x) && theta(y,x) < 11*pi/8)    %/%202.5 < t <247.5;
                l = img_mag(y-1, x-1);                                % check bot left and top right
                r = img_mag(y+1, x+1);
            elseif(11*pi/8 <= theta(y,x) && theta(y,x) < 13*pi/8)   %|%247.5 < t <292.5;
                l = img_mag(y-1, x);                                  % check bot and top
                r = img_mag(y+1, x);
            elseif(13*pi/8 <= theta(y,x) && theta(y,x) < 15*pi/8)   %\%292.5 < t <337.5;
                l = img_mag(y+1, x-1);                                % check top left and bot right
                r = img_mag(y-1, x+1);
            end
        
            if((img_mag(y,x) >= l) && (img_mag(y,x) >= r))
                img_Sup(y,x) = sel; %if this is the highest pixel, save it in suppressed image
            else
                img_Sup(y,x) = 0;   %else, dont save this pixel
            end
        end
    end
end
function img_thresh = doThresh(img_Sup, lo, hi, min, max)
    img_thresh = zeros(size(img_Sup)); %initiliaze thresholded image
    img_hi = (img_Sup > hi);    % confident pixels above high threshold
    img_lo = lo < img_Sup & img_Sup < hi; % weak pixels between lo and high 
    for i = 1 : size(img_thresh, 1)
        for j =  1 : size(img_thresh, 2)
            if(img_hi(i,j) == 1)
                img_thresh(i,j) = max;
            end
            if(img_lo(i,j) == 1)
                img_thresh(i,j) = min;  % for 2.1 dont implement hysteresis, so anything below hi thresh, set to 0
                                % in the Hysteresis step, set this number to a gray value like 60
            end
        end
    end
end
function img_thresh = drawHysteresis(img_thresh)
    for i = 1 : size(img_thresh, 1)
        for j =  1 : size(img_thresh, 2)
            converted = false;
            if(img_thresh(i,j) == 60)
                if(img_thresh(i,j+1) == 255)%check neighbours, then draw
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i+1,j+1) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i+1,j) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i+1,j-1) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i,j-1) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i-1,j-1) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i-1,j) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(img_thresh(i-1,j+1) == 255)
                    img_thresh(i,j) = 255;
                    converted = true;
                end
                if(converted == false)
                    img_thresh(i,j) = 0; 
                end
            end
        end
    end
end