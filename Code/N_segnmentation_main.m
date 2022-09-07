clear all
close all
%% input
path_directory='C:\Users\Nicky Nirlipta\Desktop\Nicky Nirlipta Sahoo_EE19S042-ISP PROJECT\SEGMENTED OUTPUT\color 51-75'; % 'Folder name'
original_files=dir([path_directory '/*.jpg']);
%image_orginal='C:\Users\Nicky Nirlipta\Desktop\Nicky Nirlipta Sahoo_EE19S042-ISP PROJECT\SEGMENTED OUTPUT\color 26-50\189080.JPG'; 
for k=1     %length(original_files)
filename=[path_directory '/' original_files(k).name];
image_orginal=imread(filename);
%figure,imshow(image_orginal),title(original_files(k).name);
IM = image_orginal;
%IM    = imread('C:\Users\Nicky Nirlipta\Desktop\Nicky Nirlipta Sahoo_EE19S042-ISP PROJECT\SEGMENTED OUTPUT\Color 1-25\101085.jpg'); 
I=imresize(IM,0.5);
% ncut parameters
for i=9  %12  %set the value for SI
 for j=8   %set the value for SX

tic  
sb=4;       %input("the number image to be segmented:");
SI   = i;                       % Color similarity
SI
SX   = j;                       % Spatial similarity
SX
r    =2;                        % Spatial threshold (less than r pixels apart)
sNcut = 0.21;                   % The smallest Ncut value (threshold) to keep partitioning
sArea = 320;%:100:220;                     % smallest area
figure()
[Inc2 Inc, Nnc]   = N_weightmatrix(I,SI,SX,r,sNcut,sArea,sb);   % Normalized Cut
%% show
toc;
title(['seg SI=',num2str(SI),'SX=',num2str(SX),'r=',num2str(r)]);
           %sprintf['seg SI=',num2str(i),'SX=',num2str(j)].jpg)

figure()
subplot(121); imshow(IM);    title('Original'); 
% subplot(121),imshow(image_orginal),title(original_files(k).name);
 subplot(122); imshow(Inc);  title(['NormalizedCut',' : ',num2str(Nnc)]); 
 
        
    end
    disp('enter to continue:');
    pause;
    
end

 disp(strcat('This is  image  segment for ', num2str(k), ' image, press Enter to continue...'));   
 pause;
 close all;
 end

