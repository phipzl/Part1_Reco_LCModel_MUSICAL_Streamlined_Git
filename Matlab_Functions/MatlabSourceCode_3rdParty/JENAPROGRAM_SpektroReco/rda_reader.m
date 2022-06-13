function [d,h,D,x_axes] = rda_reader(rda_name)

if nargin < 1
    cur_dir =pwd;
    [n,pf]  = uigetfile('*.*','Oeffnen einer rda-Datei');
    if length(n) < 2 && length(pf) < 2
        cd(cur_dir);
        d           = 0;
        h           = 0;
        D           = 0;
        x_axes      = 0;
        return;
    end
    rda_name    = strcat(pf,n);
    cd (cur_dir);
end

% *************************************************************************
tic
fid         = fopen(rda_name,'r');
data_char   = fscanf(fid,'%c');
fclose(fid);
toc

word_length     = 23;
for j = 1:length(data_char)-word_length;
    word            = data_char(j:j+word_length-1);
    if strcmp(word, cat(2,'>>> End of header <<<',char(13),char(10)))
        start   = j+word_length-1;
        header  = data_char(1:start);
        break;
    end
end

h   = read_header(header);
M   = max(size(h));

f       = 1;
fcalib  = 1;

try
    for j = 1:M
        if strcmp(h{j,1},'MRFrequency'),                                        w0              = 1000000*str2double(h{j,2});   end
        if strcmp(h{j,1},'DwellTime'),                                          dt              = 0.000001*str2double(h{j,2});  end
        if strcmp(h{j,1},'ProtocolName'),                                       seq             = h{j,2};                       end
        if strcmp(h{j,1},'VectorSize'),                                         N               = 2*str2double(h{j,2});         end
        if strcmp(h{j,1},'CSIMatrixSize[0]'),                                   N_CSI_colums    = str2double(h{j,2});           end
        if strcmp(h{j,1},'CSIMatrixSize[1]'),                                   N_CSI_rows      = str2double(h{j,2});           end
        if strcmp(h{j,1},'CSIMatrixSize[2]'),                                   N_CSI_slices    = str2double(h{j,2});           end
        if strcmp(h{j,1},'TransmitRefAmplitude[1H]'),                           TRAMP           = str2double(h{j,2});           end
        if strcmp(h{j,1},'FoVWidth'),                                           Vox1_csi        = str2double(h{j,2});           end
        if strcmp(h{j,1},'FoVHeight'),                                          Vox2_csi        = str2double(h{j,2});           end
        if strcmp(h{j,1},'FoV3D'),                                              Vox3_csi        = str2double(h{j,2});           end
        if strcmp(h{j,1},'VOIReadoutFOV') || strcmp(h{j,1},'VOIReadoutVOV'),    Vox1_svs        = str2double(h{j,2});           end
        if strcmp(h{j,1},'VOIPhaseFOV'),                                        Vox2_svs        = str2double(h{j,2});           end
        if strcmp(h{j,1},'VOIThickness'),                                       Vox3_svs        = str2double(h{j,2});           end
        % if strcmp(h{j,1},'MagneticFieldStrength'),                              B0              = str2double(h{j,2});      end
    end
    
    if  ~exist('TRAMP','var'),          TRAMP = 1;          end
    if  ~exist('N_CSI_colums','var'),   N_CSI_colums    = 1; end
    if  ~exist('N_CSI_rows','var'),     N_CSI_rows      = 1; end
    if  ~exist('N_CSI_slices','var'),   N_CSI_slices    = 1; end        

    for j = 1:length(seq)-2
        if      strcmp(seq(j:j+2),'svs') || strcmp(seq(j:j+2),'fid')
            if ~exist('Vox1_svs','var') || ~exist('Vox2_svs','var') || ~exist('Vox3_svs','var')
                VOLUME  = 1;
            else
                VOLUME  = (0.1*Vox1_svs)*(0.1*Vox2_svs)*(0.1*Vox3_svs); break;
            end
        elseif  strcmp(seq(j:j+2),'csi')            
            if ~exist('Vox1_csi','var') || ...
                      ~exist('Vox2_csi','var') || ...
                      ~exist('Vox3_csi','var') || ...
                      ~exist('N_CSI_colums','var') || ...
                      ~exist('N_CSI_rows','var') || ...
                      ~exist('N_CSI_slices','var')
                VOLUME  = 1;
            else
                VOLUME  = (0.1*Vox1_csi/N_CSI_colums)*(0.1*Vox2_csi/N_CSI_rows)*(0.1*Vox3_csi/N_CSI_slices); break;
            end
        end
    end
catch
    d           = 0;
    h           = 0;
    D           = 0;
    x_axes      = 0;
    return;
end

% Data scalling
try
    f       = TRAMP/VOLUME;
catch
    f       = 1;
end

% if      B0 > 1.3 & B0 < 1.7,    fcalib = 440.67;
% elseif  B0 > 2.8 & B0 < 3.2,    fcalib = 1199.17;
% end

% calculate t, w- and ppm-axes

x_axes.time     = (0:N-1).*dt;
dw              = (1/dt)/(0.5*N);
w               = (w0-0.5/dt) + (0:(0.5*N)-1).*dw;
w_rel           = -0.5/dt + (0:(0.5*N)-1).*dw;
ppm             = (1000000.*w_rel)./w0;

k   = 1;
for i1 = 0.5*N:-1:1
    x_axes.w(k)         = w(i1);
    x_axes.w_rel(k)     = w_rel(i1);
    x_axes.ppm(k)       = ppm(i1);
    k                   = k+1;
end


% *************************************************************************

fid         = fopen(rda_name,'r');
fseek(fid, start, -1);          % read original SIEMENS rda-Data
full_data   = fread(fid,N*N_CSI_colums*N_CSI_rows*N_CSI_slices,'double');
fclose(fid);

% *************************************************************************

start       = 1;

k1  = 1;
for k = 1:N_CSI_slices
    j1  = 1;
    for j = 1:N_CSI_rows
        i1  = 1;        
        for l = 1:N_CSI_colums
            temp            = full_data(start:start+N-1);
            d{k1,j1,i1}     = complex(temp(1:2:N),temp(2:2:N));
            D{k1,j1,i1}     = fftshift(fft(d{k1,j1,i1}));
            start           = start+N;
            i1              = i1+1;
        end
        j1  = j1+1;
    end
    k1  = k1+1;
end

% % *************************************************************************
% % Skalierung
% 
% if strcmp(questdlg('Save scaled data?'),'Yes')
%     new_rda_name    = strcat(rda_name(1:length(rda_name)-4),'_mod.rda');
%     rda_saver((f*fcalib).*full_data,header,new_rda_name,'double');
% end