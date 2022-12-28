%Gray_Image
% grayImage = imread("kid blurred-noisy.tif");
grayImage = imread("fruit blurred-noisy.tif");
imshow(grayImage);
ax = gcf;
exportgraphics(ax, "Gray_Image.png",'Resolution', 200);

% Laplacian_Gradient
% laplacianKernel = [-1,-1,-1;-1,8,-1;-1,-1,-1]; % kid
laplacianKernel = [0, 1, 0; 1,-4, 1; 0, 1, 0]; % fruit
laplacianImage = imfilter(double(grayImage), laplacianKernel, "replicate");
laplacianImag = uint8(255*mat2gray(laplacianImage));
imshow(laplacianImag);
ax = gcf;
exportgraphics(ax, "Laplacian_Gradient.png",'Resolution', 200);

% Laplacian_Sharpen
sharpenedImage = double(grayImage) + laplacianImage;
sharpenedImag = uint8(255*mat2gray(sharpenedImage));
imshow(sharpenedImag);
ax = gcf;
exportgraphics(ax, "Laplacian_Sharpen.png",'Resolution', 200);

% Sobel_Gradient
[magnitudeImage, directionImage] = imgradient(grayImage, 'Sobel');
magnitudeImag = uint8(255*mat2gray(magnitudeImage));
imshow(magnitudeImag);
ax = gcf;
exportgraphics(ax, "Sobel_Gradient.png",'Resolution', 200);

% Smooth_Gradient
SmoothGradient = imboxfilt(magnitudeImage,5)/25;
SmoothGradien = uint8(255*mat2gray(SmoothGradient));
imshow(SmoothGradien);
ax = gcf;
exportgraphics(ax, "Smooth_Gradient.png",'Resolution', 200);

% Extracted_Feature
ExtractedFeature = immultiply(SmoothGradient,laplacianImage);
ExtractedFeatur = uint8(255*mat2gray(ExtractedFeature));
imshow(ExtractedFeatur);
ax = gcf;
exportgraphics(ax, "Extracted_Feature.png",'Resolution', 200);

% A_Plus_F
AF = double(ExtractedFeatur) + double(grayImage);
AF1 = uint8(255*mat2gray(AF));
imshow(AF1);
ax = gcf;
exportgraphics(ax, "A_Plus_F.png",'Resolution', 200);

% Powerlaw_Transformation
const = 1;
gama = 2;
image = double(AF);
S = const * (image .^gama);
% T = 255/(const * (255 .^gama));
% powerlaw = uint8(T * S);
powerlaw = uint8(255*mat2gray(S));
imshow(powerlaw);
ax = gcf;
exportgraphics(ax, "Powerlaw_Transformation.png",'Resolution', 200);