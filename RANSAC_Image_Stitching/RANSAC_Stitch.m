%%  RANSAC Image Stitching
%   this matlab script taking in a pair of images "parliament-left.jpg" and
%   "parliament-right.jpg" that contain the same object, taken at different
%   angles. The result, "parliament-mosaic.jpg" contains the mosaic created
%   from these two images

%% A.1
%pLeft      contains left image, converted to grayscale and single
%pRight     contains right image, converted to grayscale and single
tic;
pLeft = single(rgb2gray(imread("parliament-left.jpg")));
pRight = single(rgb2gray(imread("parliament-right.jpg")));

%% A.2
%fa,fb      matrix that contains [x;y;s;th]
%               x,y contain the centre of the frame
%               s   contain the scale of the frame
%               th  contain the orientation of the frame
%da,db      128D vector per frame that contains 4x4 grid of gradient
%           orientation, and 2x2x2 histogram bins to reduce misestimation

[fa, da] = vl_sift(pLeft);
[fb, db] = vl_sift(pRight);

%% A.3, A.4
%matches    2 by N matrix that contains N pairs of likely matches.
%           1st row stores indices of fa/da, 2nd row stores indices of
%           fb/db. These are in pairs based on their similarity
%scores     1 by N matrix that contains the squared euclidean distance 
%           between matching frames. A lower score/distance corressponds to         
%           a greater match.

[matches, scores] = vl_ubcmatch(da,db, 750);

%% Statistics checker code
%means = zeros(1,188);
%for i = 1:188
%    means(1,i) = mean(da(:,matches(1, i)) == db(:, matches(2,i)));
%end
%mean(means)        %mean equivalence rate of matching frames
%
%r1 = randi([1, 28185], 500);
%r2 = randi([1, 22590], 500);
%means2 = zeros(1, 500);
%for i = 1:500
%       means2(1,i) = mean(da(:,r1(1,i)) == db(:,r2(1,i))); 
%end
%mean(means2)       %mean equivalence rate of random frames
%% A.5
iters = 188;
p = randperm(size(scores, 2), iters); %select random index 1-188 of matches
rho = 1500;
bestCount = 0;
bestInliers = zeros(size(matches));
%timeshit=0;
%minDist = 9999;
%maxDist = 0;

for i = 1:iters    %188 iterations of RANSAC
    thisCount = 0;
    thisInliers = zeros(size(matches));
    %check a random match:
    l = matches(1, p(i)); %pixel match in left image
    r = matches(2, p(i)); %pixel match in right image
    tx = fb(1, r) - fa(1, l); %horizontal translation
    ty = fb(2, r) - fa(2, l); %vertical translation
    s = fb(3, r)/fa(3, l); %scaling factor
    a = -(fb(4, r)-fa(4, l)); %rotation factor
    T = [s*cos(a) sin(a) 0; -sin(a) s*cos(a) 0; -tx -ty 1];
    Trans = affine2d(T);
    %[out, out_ref] = imwarp(pRight, Trans);
    
    
    %need to check how close matches from fb are to fa after transforming
    %it. apply transform to fb, compare to fa
    for j = 1:size(scores, 2)
        newXY = T * [fb(1, matches(2, j)); fb(2, matches(2, j)); 1];
        dx = fa(1, matches(1,j)) - newXY(1);
        dy = fa(2, matches(1,j)) - newXY(2);
        dist = (dx^2 + dy^2)^0.5;
        %if(dist<minDist)
        %   minDist = dist; 
        %end
        if(dist < rho) %if true, this is an inlier
            %if(dist>maxDist)
            %    maxDist = dist;
            %end
            thisCount = thisCount+1;
            thisInliers(1:2, thisCount) = matches(:, j);
        end
    end
    if bestCount < thisCount %if this has more inliers, then it is new best
        bestCount = thisCount;
        bestInliers = thisInliers;
        %bestT = T;
        %timeshit = timeshit + 1;
    end
end
%% A.6
A = [];
for i = 0:bestCount-1   %construct repeating identity matrix
    A(1 + 4*i:4 + 4*i, 1:4) = eye(4);
end
Y = zeros(bestCount*4, 1);
for k = 0:bestCount-1   %for each inlier, get the values for scale, 
                        %rotation, and translation (x and y)
    l = bestInliers(1, k+1);
    r = bestInliers(2, k+1);
    Y(1+4*k) =  fb(3, r)/fa(3, l);
    Y(2+4*k)= -(fb(4, r)-fa(4, l));
    Y(3+4*k)= fb(1, r) - fa(1, l);
    Y(4+4*k)= fb(2, r) - fa(2, l);
end
X = pinv(A)*Y; %[s; a; tx; ty]
FinalT =    [(1/X(1))*cos(X(2)) sin(X(2))           0; 
            -sin(X(2))          (1/X(1))*cos(X(2))  0; 
            -X(3)-33            -X(4)+14            1];        
%% A.7
%Terribly slow code
%This code will emulate the final mosaic in grayscale, but mimic the mosaic
%in RGB. i.e. the max pixel value between both images will be selected in 
%grayscale, and that corresponding pixel will be chosen for the RGB output.
%


%Open colour images and create references for all crucial images
pRight_col = imread("parliament-right.jpg");
pLeft_col = imread("parliament-left.jpg");
pRight_ref = imref2d(size(pRight));
pLeft_ref = imref2d(size(pLeft));
pRight_col_ref = imref2d(size(pRight_col));
pLeft_col_ref = imref2d(size(pLeft_col));
bg = zeros(5000);
bg_ref = imref2d(size(bg));

%create Transformation
FinalTrans = affine2d(FinalT);


%Transform pRight
[pRight_Translated,pRight_translated_ref] = imwarp(pRight,FinalTrans);
[pRight_col_Trans, pRight_col_Trans_ref]  = imwarp(pRight_col, FinalTrans);


%put pLeft and pRight Transformed onto black backgrounds of identical size
[bg1, bg1_ref] = imfuse(bg, bg_ref, pRight_Translated, pRight_translated_ref, 'blend', 'Scaling', 'joint');
[bg2, bg2_ref] = imfuse(bg, bg_ref, pLeft, pLeft_ref, 'blend', 'Scaling', 'joint');
[bg1_col, bg1_col_ref] = imfuse(bg, bg_ref, pRight_col_Trans, pRight_col_Trans_ref, 'blend', 'Scaling', 'joint');
[bg2_col, bg2_col_ref] = imfuse(bg, bg_ref, pLeft_col, pLeft_col_ref, 'blend', 'Scaling', 'joint');
bg1_col = bg1_col(:,:,:).*2; %brighten the image
bg2_col = bg2_col(:,:,:).*2; %brighten the image

cutoff1 = size(bg1, 1); %will determine where to crop the image
cutoff2 = size(bg1, 2); %will determine where to crop the image
CompositeFinal = [];
for i = 1:cutoff1
    %if this row consists of only black pixels, it should be cropped
    if(isequal(bg1(i, :), zeros(1, size(bg1, 1))) && isequal(bg2(i, :), zeros(1, size(bg1, 1))))
        continue 
    end
    for j = 1:cutoff2
        if(bg1(i,j) > bg2(i,j))
            %keep pixel from pRight
            CompositeFinal(i,j, 1:3) = double(bg1_col(i,j, 1:3))/255;
        else
            %keep pixel from pLeft
            CompositeFinal(i,j, 1:3) = double(bg2_col(i,j, 1:3))/255;
        end
        %if this column consists of only black pixels, it should be cropped
        if(isequal(bg1(:, j), zeros(size(bg1, 1),1)) && isequal(bg2(:, j), zeros(size(bg1, 1), 1)))
            cutoff2 = j;
            break
        end
    end
    %imshow(CompositeFinal, []);
end
toc;
imshow(CompositeFinal, []);