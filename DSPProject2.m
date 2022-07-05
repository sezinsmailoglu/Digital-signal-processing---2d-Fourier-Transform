%%
N_IMAGES = 40;
images = cell(N_IMAGES);

for i = 1 : N_IMAGES
    t = Tiff(sprintf('%i.tif', i) );
    images{i} = read(t);
    close(t);
end

%%
basis_original =  zeros( N_IMAGES, numel(images{1} ) );
basis_offset = zeros( N_IMAGES, 1  );
basis = zeros( N_IMAGES, numel(images{1} ) ); % images are stored in rows

for i = 1 : N_IMAGES
    basis_original(i,:) = reshape( double(images{i}), 1, []);
    basis_offset(i) = mean(basis_original(i,:));
    basis(i,:) = basis_original(i,:) - basis_offset(i);
end



%% orthogonality check
TOL = 1e-7;
A = basis * basis';

fprintf('Energy in the cross terms over the main components = %f\n', sum(abs(A - diag(diag(A)) ),'all') / sum(abs(diag(A)))  );
% Result is 0.000088 hence almost all the vectors are
% orthogonal, we can make projections to calculate the coefficients



%%
%get cat
t = Tiff('cat.tif','r');
cat_original = reshape( double(read(t)), 1, 512 * 512); % read as a row vector
close(t);
cat_offset = mean(cat_original);
cat =  cat_original - cat_offset; 
%% calculate the basis coefficients
ak =  inv(A) * ( basis * cat' );
 
% b is the total offset - minus offsets from X_k's
b = cat_offset - basis_offset(1:35)' * ak(1:35) ; 

% reconstruction from the first-35 basis components 
c_hat_1 = ( basis(1:35,:)' * ak(1:35) + cat_offset * ones(512 * 512, 1) )'; % 
% the energy error percentage:
fprintf('Energy error percentage  1: %f\n', 100 * sum(  (cat_original - c_hat_1).^2 ) / sum(cat_original.^2 ) );



%% Problem-2 
% we look the first 35 components where the projetion of the cat is the
% largest. since the basis' are orthogonal doing that is just finding the
% largest 35 components, hence we sort
sk = (basis * basis')^(-1/2) * (basis * cat'); % projection length in each basis 
[~, I] = sort(sk, 'descend');
I_35  = I(1:35);
fprintf('Replaced images:\n');
fprintf('%d, ', sort(I(36:40)) );
fprintf('\n');
b_2 = cat_offset - basis_offset(I_35)' * ak(I_35);
c_hat_2 = ( basis(I_35,:)' * ak(I_35) + cat_offset * ones(512 * 512, 1) )';

% the energy error percentage:
fprintf('Energy error percentage 2: %f\n', 100 * sum(  (cat_original - c_hat_2).^2 ) / sum(cat_original.^2 ) );

%% Out of curiosity
c_hat_3 = ( basis' * ak + cat_offset * ones(512 * 512, 1) )';

% the energy error percentage:
fprintf('Energy error percentage 3: %f\n', 100 * sum(  (cat_original - c_hat_3).^2 ) / sum(cat_original.^2 ) );

% demonstrate
figure;
subplot(2,2,1);
imshow(reshape(cat,512,512),[]);
title('Original Image');
subplot(2,2,2);
imshow(reshape(c_hat_1,512,512),[]);
title('First 35 components');
subplot(2,2,3);
imshow(reshape(c_hat_2,512,512),[]);
title('Highest 35 components');
subplot(2,2,4);
imshow(reshape(c_hat_2,512,512),[]);
title('All 40 components');