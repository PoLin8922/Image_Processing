%% Clear the environment and the command line
clc;
close all;
clear;

%% read image
kid_img = imread('kid.tif');
fruit_img = imread('fruit.tif');
[m,n] = size(kid_img);

%% processing
%padded image size 2m*2n
GLPF = zeros(2*m,2*n);
kid_img2 = zeros(2*m,2*n);
fruit_img2 = zeros(2*m,2*n);

% cutoff frequency
d0 = 100*2;

%set Gaussian LPF and HPF
for i = 1:2*m
    for j = 1:2*n
        d = (i-m).^2+(j-n).^2;
        GLPF(i,j) = exp(-d/2/d0/d0);
        if(i<=m && j<=n)
            kid_img2(i,j) = kid_img(i,j);
            fruit_img2(i,j) = fruit_img(i,j);
        end
    end
end

GHPF = 1-GLPF;


% centering processing 
% multiply(-1)^(x+y)
for i = 1:2*m
    for j = 1:2*n
        kid_img2(i,j) = kid_img2(i,j)*((-1)^(i+j-2));
        fruit_img2(i,j) = fruit_img2(i,j)*((-1)^(i+j-2));
    end
end


%use DFT of the image to multiply LPF or HPF
kid_dft_img = fft2(kid_img2);
fruit_dft_img = fft2(fruit_img2);

kid_lpf_img = kid_dft_img.*GLPF;
kid_hpf_img = kid_dft_img.*GHPF;

fruit_lpf_img = fruit_dft_img.*GLPF;
fruit_hpf_img = fruit_dft_img.*GHPF;

%the result is (-1)^(x+y)*Re{IDFT(result after LPF or HPF)}
kid_lpf_idft = real(ifft2(kid_lpf_img));
kid_hpf_idft = real(ifft2(kid_hpf_img));

fruit_lpf_idft = real(ifft2(fruit_lpf_img));
fruit_hpf_idft = real(ifft2(fruit_hpf_img));

% crop M*N image
kid_final_lpf_img = zeros(m,n);
kid_final_hpf_img = zeros(m,n);

fruit_final_lpf_img = zeros(m,n);
fruit_final_hpf_img = zeros(m,n);
for i = 1:m
    for j = 1:n
        kid_final_lpf_img(i,j) = kid_lpf_idft(i,j)*((-1)^(i+j-2));
        kid_final_hpf_img(i,j) = kid_hpf_idft(i,j)*((-1)^(i+j-2));

        fruit_final_lpf_img(i,j) = fruit_lpf_idft(i,j)*((-1)^(i+j-2));
        fruit_final_hpf_img(i,j) = fruit_hpf_idft(i,j)*((-1)^(i+j-2));
    end
end


%calculate for top25 DFT frequency pairs
kid_mag = abs(fftshift(fft2(kid_img)));
sortlist_kid = zeros(m,n/2);
kid_top25_freq_pair = zeros(25,2);

fruit_mag = abs(fftshift(fft2(fruit_img)));
sortlist_fruit = zeros(m,n/2);
fruit_top25_freq_pair = zeros(25,2);

for i = 1:m
    for j = 1:n/2
        sortlist_kid(i,j) = kid_mag(i,j);
        sortlist_fruit(i,j) = fruit_mag(i,j);
    end
end

[vk,indexk] = sort(sortlist_kid(:),"descend");
[vf,indexf] = sort(sortlist_fruit(:),"descend");

for i = 1:25
    if(mod(indexk(i),m)==0)
        rowk = m;
    else
        rowk = mod(indexk(i),m);
    end
    colk = ceil(indexk(i)/m);
    kid_top25_freq_pair(i,:) = [rowk,colk];

    if(mod(indexf(i),m)==0)
        rowf = m;
    else
        rowf = mod(indexf(i),m);
    end
    colf = ceil(indexf(i)/m);
    fruit_top25_freq_pair(i,:) = [rowf,colf];
end


%% output
figure(1);
imshow(kid_img,[]);

figure(2);
imshow(fruit_img,[]);

kid_log_FM = log(1+abs(fftshift(fft2(kid_img))));
figure(3);
imshow((kid_log_FM),[]);
fig= gcf;
exportgraphics(fig,'kid magnitude spectra.png','Resolution',150);

fruit_log_FM = log(1+abs(fftshift(fft2(fruit_img))));
figure(4);
imshow((fruit_log_FM),[]);
fig= gcf;
exportgraphics(fig,'fruit magnitude spectra.png','Resolution',150);

figure(5);
imshow(GLPF,[]);
fig= gcf;
exportgraphics(fig,'Magnitude responses of Gaussian LPF.png','Resolution',150);

figure(6);
imshow(GHPF,[]);
f= gcf;
exportgraphics(fig,'Magnitude responses of Gaussian HPF.png','Resolution',150);

figure(7);
imshow(kid_final_lpf_img,[]);
fig= gcf;
exportgraphics(fig,'kid after Gaussian LPF.png','Resolution',150);

figure(8);
imshow(kid_final_hpf_img,[]);
fig= gcf;
exportgraphics(fig,'kid after Gaussian HPF.png','Resolution',150);

figure(9);
imshow(fruit_final_lpf_img,[]);
fig= gcf;
exportgraphics(fig,'fruit after Gaussian LPF.png','Resolution',150);

figure(10);
imshow(fruit_final_hpf_img,[]);
fig= gcf;
exportgraphics(fig,'fruit after Gaussian HPF.png','Resolution',150);
