%%  Seam Carving
%   This function will take an image, and dimensions as input, and return 
%   the image with dimensions REDUCED to the input dimensions. This is done
%   using the seam carving algorithm

function MySeamCarving(image, h, w, imgOut)
    img = imread(image);
    for k = size(img,2) - h :-1: 1
        img = CarvingHelper(img);
        fprintf("finished vertical seam: %d\n", k);
    end
    img = pagetranspose(img);
    for l = size(img,2) - w:-1: 1
       img = CarvingHelper(img);
       fprintf("finished horizontal seam: %d\n", l);
    end
    img = pagetranspose(img);
    
    imwrite(img, imgOut)
       
    
end

function img2 = CarvingHelper(img)
    h = fspecial('sobel');
    dxR = imfilter(double(img(:,:,1)), h', 'conv');
    dyR = imfilter(double(img(:,:,1)), h, 'conv');
    dxG = imfilter(double(img(:,:,2)), h', 'conv');
    dyG = imfilter(double(img(:,:,2)), h, 'conv');
    dxB = imfilter(double(img(:,:,3)), h', 'conv');
    dyB = imfilter(double(img(:,:,3)), h, 'conv');
    
    img_magR = sqrt(dxR.^2 + dyR.^2);  %get magnitude
    img_magG = sqrt(dxG.^2 + dyG.^2);  %get magnitude
    img_magB = sqrt(dxB.^2 + dyB.^2);  %get magnitude
    
    Energy = img_magR + img_magG + img_magB;
    
    
    M = zeros(size(Energy));
    M(1, :) = Energy(1, :);
    for y = 2:size(M,1)
        for x = 1:size(M,2)
            if(x == 1)
                M(y,x) = Energy(y,x) + min([M(y-1,x), M(y-1, x+1)]);
            elseif(x == size(M,2))
               M(y,x) = Energy(y,x) + min([M(y-1,x-1), M(y-1,x)]); 
            else
                M(y,x) = Energy(y,x) + min([M(y-1,x-1), M(y-1,x), M(y-1, x+1)]);
            end
        end
    end   
    lowest = min(M(size(M, 1), :));
    for i  = 1:size(M, 2)
        if(M(size(M,1), i) == lowest)
           start = i; 
        end
    end
    x = start;
    for j = size(M, 1):-1:2
        path(j, 1) = j;
        path(j, 2) = x;
        if(x == 1)
            if(M(j-1, x) == min([M(j-1, x) M(j-1, x+1)]))
                %value above is min
                continue
            else
                %value above to the right is min
                x = x+1;
                continue
            end
        elseif(x == size(M, 2))
            if(M(j-1, x-1) == min([M(j-1, x-1) M(j-1,x)]))
                %value above to the left is min
                x = x-1;
                continue;
            else
                %value above is min
                continue
            end
        else
            if(M(j-1, x-1) == min([M(j-1, x-1) M(j-1, x) M(j-1, x+1)]))
                %value above to the left is min
                x = x-1;
                continue
            elseif(M(j-1, x) == min([M(j-1, x-1) M(j-1, x) M(j-1, x+1)]))
                %value above is min
                continue
            else
                %value above to the right is min
                x = x+1;
                continue
            end
        end
    end
    path(1, 1) = 1;
    path(1, 2) = x;
    
    img2 = zeros(size(img,1), size(img,2)-1, 3, 'uint8');
    for j = 1:size(img, 1) %size(img2,1)
        for i = 1:size(img, 2)-1 %size(img2, 2)
            if(i == path(j,2))
                img2(j, i:size(img,2)-1, :) = img(j,i+1:size(img,2),:);
                break
            else
                img2(j,i,:) = img(j,i,:);
            end
        end
    end
end

