% This example demonstrates how to use pentary MVG embedding function
clc
clear all
close all

% Read the input cover image
Cover = double(imread ('1.pgm'));

% Set the payload to 0.4 bpp
Payload = 0.4;

% Pentary MVG embedding
tStart_P = tic;
[Stego_P, pChange_P, ChangeRate_P] = MVG( Cover, Payload );
tEnd_P = toc(tStart_P);
fprintf('Pentary MVG embedding is done in: %f (sec)\n',tEnd_P);
%%
close all

figure;
imshow (Cover,[]);
title ('Cover image');

figure;
subplot(1,2,1);
imshow(1-pChange_P.PM1/0.3333);
title('Embedding change probabilities for +-1');
subplot(1,2,2);
imshow(1-pChange_P.PM2/0.3333);
title('Embedding change probabilities for +-2');

figure;
T = abs(Stego_P-Cover);T(T==0)=0.5;T(T==1)=0;T(T==2)=1;
imshow(T);
title('Changed pixels (+-1 -> black ,+-2 -> white)');