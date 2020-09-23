function [Stego, pChange, ChangeRate] = MVG( Cover, Payload )
% -------------------------------------------------------------------------
% Pentary Multivariate Gaussian Embedding | September 2015 | version 1.1 
% -------------------------------------------------------------------------
% INPUT:
%  - Cover - Path to the cover image or the cover image itself.
%  - Payload - Embedding payload in bits per pixel (bpp).
% OUTPUT:
%  - Stego - Resulting image with embedded payload
%  - pChange - Embedding change probabilities. It has two fields, 
%              probabilities of change by +1 or -1, 'PM1', and
%              probabilites  of change by +2 or -2, 'PM2'
%  - ChangeRate - This structure of change rates also has two fields,
%                 'PM1' for ternary change rate and 'PM2' for pentary
% -------------------------------------------------------------------------
% Copyright (c) 2015 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Contact: vsedigh1@binghamton.edu | fridrich@binghamton.edu
%          September 2015
%          http://dde.binghamton.edu/download/
% -------------------------------------------------------------------------
% References:
% [1] - V. Sedighi, J. Fridrich and R. Cogranne. Content-Adaptive Pentary 
% Steganography Using the Multivariate Generalized Gaussian Cover Model. 
% Proc. SPIE, Electronic Imaging, Media Watermarking, Security, and 
% Forensics 2015, vol. 9409, San Francisco, CA, February 8–12, 2015.
% -------------------------------------------------------------------------

% Read and convert the input cover image into double format
if ischar( Cover )
    Cover = double( imread(Cover) );
else
    Cover = double( Cover );
end

% Compute Variance and do the flooring for numerical stability
WienerResidual = Cover - wiener2(Cover,[2,2]);
Variance = VarianceEstimationDCT2D(WienerResidual,9,9);
Variance(Variance< 0.01) = 0.01;

% Compute Fisher information and smooth it
FisherInformation = 1./Variance.^2;
FisherInformation = imfilter(FisherInformation,fspecial('average',7),'symmetric');

% Compute embedding change probabilities and execute embedding
FI = FisherInformation(:)';

% Form the Fisher information matrix for pentary embedding
FI = repmat([1 16 4],numel(FI),1).*repmat(FI',1,3);
% Penrary embedding change probabilities
[beta, theta] = PentaryProbs(FI,Payload);
% Simulate embedding
Stego = Cover;
beta = 2 * beta;
theta = 2 * theta;
r = rand(1,numel(Cover));
ModifPM2 = (r < theta);                % Cover elements to be modified by +-2
ModifPM1 = (r > theta) & (r < beta);   % Cover elements to be modified by +-1
r = rand(1, numel(Cover));
Stego(ModifPM1) = Cover(ModifPM1) + 2*(round(r(ModifPM1))) - 1; % Modifying X by +-1
Stego(ModifPM2) = Cover(ModifPM2) + 4*(round(r(ModifPM2))) - 2; % Modifying X by +-2
Stego(Stego>255) = 254;                    % Taking care of boundary cases
Stego(Stego<0)   = 1;
ChangeRate.PM1 = sum(ModifPM1(:))/numel(Cover); % Computing the change rate
ChangeRate.PM2 = sum(ModifPM2(:))/numel(Cover); % Computing the change rate
pChange.PM1 = reshape(beta/2,size(Cover));
pChange.PM2 = reshape(theta/2,size(Cover));
    
end

% Beginning of the supporting functions 

% Estimation of variance of pixels based on a 2D-DCT 
% (trigonometric polynomial) model. See Ref. [1] for details.
function EstimatedVariance = VarianceEstimationDCT2D(Image, BlockSize, Degree)

% verifying the integrity of input arguments 
if ~mod(BlockSize,2) 
    error('The block dimensions should be odd!!'); 
end
if (Degree > BlockSize)
    error('Number of basis vectors exceeds block dimension!!'); 
end

% number of parameters per block
q = Degree*(Degree+1)/2;

% Build G matirx
BaseMat = zeros(BlockSize);BaseMat(1,1) = 1;
G = zeros(BlockSize^2,q);
k = 1;
for xShift = 1 : Degree
    for yShift = 1 : (Degree - xShift + 1)
        G(:,k) = reshape(idct2(circshift(BaseMat,[xShift-1 yShift-1])),BlockSize^2,1);
        k=k+1;
    end
end

% Estimate the variance
PadSize = floor(BlockSize/2*[1 1]);
I2C = im2col(padarray(Image,PadSize,'symmetric'),BlockSize*[1 1]);
PGorth = eye(BlockSize^2) - (G*((G'*G)\G'));
EstimatedVariance = reshape(sum(( PGorth * I2C ).^2)/(BlockSize^2 - q),size(Image));

end

%Compute pentray embedding change probabilities
function [beta, theta] = PentaryProbs (FI, alpha)

% beta ... embedding change probs by +-1
% theta ... embedding change probs by +-2

% Absolute payload in nats
payload = alpha * size(FI,1) * log(2);

% Initial search interval for lambda
[L, R] = deal (10^3, 10^6);

fL = h_pent(Newton_Solver(L,FI)) - payload;
fR = h_pent(Newton_Solver(R,FI)) - payload;
% If the range [L,R] does not cover alpha enlarge the search interval
while fL*fR > 0             
    if fL > 0
        R = 2*R;
        fR = h_pent(Newton_Solver(R,FI)) - payload;
    else
        L = L/2;
        fL = h_pent(Newton_Solver(L,FI)) - payload;
    end
end

% Search for the labmda in the specified interval
[i, fM, TM] = deal(0, 1, zeros(60,2));
while (abs(fM)>0.0001 && i<60) 
    M = (L+R)/2;
    fM = h_pent(Newton_Solver(M,FI)) - payload;
    if fL*fM < 0, R = M; fR = fM;
    else          L = M; fL = fM; end
    i = i + 1;
    TM(i,:) = [fM,M];
end
if (i==60)
    M = TM(find(abs(TM(:,1)) == min(abs(TM(:,1))),1,'first'),2); 
end

% Compute beta and theta using the found lambda
Probs = Newton_Solver(M,FI);
beta = Probs(1,:);
theta = Probs(2,:); 

end

% Pentary Newton solver parallelized over all pixels
% See Ref. [1] for details
function [Probs] = Newton_Solver(Lambda,FI)

accuracy = 10^-6;

Probs = 0.001 * ones(2,size(FI,1));
Ind = true(1,size(FI,1));

while sum(Ind)>0
    Probs(Probs<0) = 10^-10;
    beta  = Probs(1,Ind);
    theta = Probs(2,Ind);
    
    F1 = FI(Ind,1)'.*beta + FI(Ind,3)'.*theta - (1/Lambda)*log((1-2*beta-2*theta)./beta);
    F2 = FI(Ind,2)'.*theta + FI(Ind,3)'.*beta - (1/Lambda)*log((1-2*beta-2*theta)./theta);
    
    M11 = FI(Ind,1)'-(1/Lambda)*((2*theta-1)./(beta.*(1-2*theta-2*beta)));
    M22 = FI(Ind,2)'-(1/Lambda)*((2*beta-1)./(theta.*(1-2*theta-2*beta)));
    M12 = FI(Ind,3)'+(1/Lambda)*(2./(1-2*theta-2*beta));
    
    detM = (M11.*M22-M12.^2);

    Probs(:,Ind) = Probs(:,Ind) - [(M22.*F1-M12.*F2)./detM; (M11.*F2-M12.*F1)./detM];

    Ind(Ind) = (abs(Probs(1,Ind)-beta)> accuracy) + (abs(Probs(2,Ind)-theta)> accuracy);
end

end

% Pentary entropy function expressed in nats
function Ht = h_pent(Probs)

p0 = 1-2*sum(Probs);
P = [p0(:);Probs(:);Probs(:)];
H = -(P .* log(P));
H((P<eps) | (P > 1-(eps))) = 0;
Ht = nansum(H);

end