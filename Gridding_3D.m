
R=sqrt(X.^2+Y.^2);

kx=zeros(ns,nc);
ky=zeros(ns,nc);

for i=1:nc
    for t=1:ns

        kx(t,i)=R(i)*cos(2*pi*(t-1)/ns);
        ky(t,i)=R(i)*sin(2*pi*(t-1)/ns);

    end;
end;


nsamp=zeros(2,ns*nc);
nsamp(1,:)=transpose(reshape((kx./max(R)/2),[nc*ns 1]));
nsamp(2,:)=transpose(reshape((ky./max(R)/2),[nc*ns 1]));
[H,APF,DCF] = gridOperator(nsamp,fov_overgrid*2*(REShalbe-1),1);

fov_overgrid=1;
csi_k_smalldummy=zeros(1,nc*ns);

csi=zeros([size(csi_k,1), 2*REShalbe, 2*REShalbe, size(csi_k,4), vs],precision);



csi_k_smalldummy(:)=reshape(squeeze(ones([1, size(csi_k,2), size(csi_k,3), 1, 1])),[1 nc*ns]);



dat=transpose(reshape(H'*(transpose(csi_k_smalldummy(:)).*DCF(:)),[fov_overgrid*2*(REShalbe-1), fov_overgrid*2*(REShalbe-1)]));


true_weights=real(dat)./HammingFilter(ones(size(dat)),[1 2]); % calc the weights for post grd dens-comp
true_weights(true_weights==inf)=1; %it seems that we regrid on values!=0 where the hamming filter is definded as 0--> const/0 =inf
true_weights(isnan(true_weights))=1; %set NaN outside of spherical k space region to zero
true_weights(true_weights==0)=1;

%apply weights
dat=dat./true_weights; % post-grd-dens-comp
dat(isnan(dat))=0; %set NaN outside of speherical k-space region to zero
dat(dat==inf)=0;


