load spiralexampledata
plot(kspacelocations);
[gdat] = gridkb(kspacelocations,spiraldata,dcf,256,1.5,2);
im = fftshift(fft2(fftshift(gdat)));
im = abs(im)/500;
cmap = [0:255].'*[1 1 1] / 256;
colormap(cmap);
figure; 
image(uint8(im));
colormap(cmap);
axis square;
