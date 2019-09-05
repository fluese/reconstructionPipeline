function settings=decorrelationMatrix(name, settings, twix_obj)
% Function to estimate noise covariance matrix based on noise data.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.kspace.decorrelateChannels
    % Estimate noise covariance R and decorrelation matrix W
    
    if strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseScan')
        if exist([name.namePath '/Noise/Slice_1/Channel_1.mat'],'file')
            for Channel = 1:settings.parameters.nCh
                tmp=load([name.namePath '/Noise/Slice_1/Channel_' num2str(Channel) '.mat']);
                if Channel == 1
                    sz    = size(tmp.data);
                    Noise = complex(zeros(settings.parameters.nCh, sz(1), sz(2), 'single'));
                end
                Noise(Channel, :, :) = tmp.data;
            end
            Noise = Noise(:,:);
        else
            settings.parameters.W = 0;
            settings.options.kspace.decorrelateChannels = false;
            return;
        end
    elseif strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseData')
        if isfield(twix_obj, 'noise')
            Noise = twix_obj.noise();
            Noise = permute(Noise, [2 1]);
        else
            settings.parameters.W = 0;
            settings.options.kspace.decorrelateChannels = false;
            return;
        end
    end
    
    % Noise: [ncoils, nsamples]
    %     Noise = Noise(:,:);
    %     R = conj(Noise) * Noise .' / size(Noise,2);
    %     [E,D] = eig(R);
    %
    %     d = diag(D);
    %     d = real(d.^-.5);
    %     D = diag(d);
    %
    %     W = E*D*E';
    
    R = conj(Noise) * Noise .' / size(Noise,2);
    R = R./mean(abs(diag(R))); % halbwegs richtige skalierung, nur fuer vergleiche wichtig
    R(eye(size(R,1))==1) = abs(diag(R));
    settings.parameters.W = sqrtm(inv(R)).';
else
    settings.parameters.W = 0;
end