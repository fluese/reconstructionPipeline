function settings=getOrientation(twix_obj, settings)
% Function to extract orientation information from raw data file.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
settings.parameters.FoV(1) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
settings.parameters.FoV(2) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
settings.parameters.FoV(3) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;

settings.parameters.Resolution(1) = settings.parameters.FoV(1) / settings.parameters.nPE;
settings.parameters.Resolution(2) = settings.parameters.FoV(2) / settings.parameters.nRO;
if isempty(twix_obj.hdr.Meas.dSliceOversamplingForDialog)
    settings.parameters.Resolution(3) = settings.parameters.FoV(3) / settings.parameters.nSlc;
else
    settings.parameters.Resolution(3) = settings.parameters.FoV(3) / settings.parameters.nSlc * (1+twix_obj.hdr.Meas.dSliceOversamplingForDialog);
end

settings.parameters.Pos(1) = twix_obj.image.slicePos(1);  %+- settings.paremeters.Resolution(1)/2;
settings.parameters.Pos(2) = twix_obj.image.slicePos(2);
settings.parameters.Pos(3) = twix_obj.image.slicePos(3);

settings.parameters.Quaternion(1) = twix_obj.image.slicePos(4);
settings.parameters.Quaternion(2) = twix_obj.image.slicePos(5);
settings.parameters.Quaternion(3) = twix_obj.image.slicePos(6);
settings.parameters.Quaternion(4) = twix_obj.image.slicePos(7);

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal, 'dSag')
    Vc_N(1) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dSag;
else
    Vc_N(1) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dSag = 0;
end

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal, 'dCor')
    Vc_N(2) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dCor;
else
    Vc_N(2) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dCor = 0;
end

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal, 'dTra')
    Vc_N(3) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dTra;
else
    Vc_N(3) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dTra = 0;
end

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}, 'dInplaneRot')
    inPlaneRot = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dInPlaneRot;
else
    inPlaneRot = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dInPlaneRot = 0;
end

settings.parameters.inPlaneRot = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dInPlaneRot;
% 
% Vc_N = Resolution(3) * Vc_N;
% 
M_R = vox2ras_rsolveAA_fl(Vc_N,inPlaneRot);
% 

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition, 'dSag')
    Vc_Ps(1) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag;
else
    Vc_Ps(1) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag = 0;
end

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition, 'dCor')
    Vc_Ps(2) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor;
else
    Vc_Ps(2) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor = 0;
end

if isfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal, 'dTra')
    Vc_Ps(3) = twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dTra;
else
    Vc_Ps(3) = 0;
    twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dTra = 0;
end
% 
M_V = vox2ras_ksolve_fl(M_R, Vc_Ps, settings.parameters.Resolution(3)/2, [settings.parameters.nRO settings.parameters.nPE settings.parameters.nSlc]);

% Setup LocationSlice:
% ---------------------------------------------------------------------
    [settings.parameters.orientation.RotationPhys2Pat, ...
     settings.parameters.orientation.RotationPhys2Log, ...
     settings.parameters.orientation.RotationPat2Phys, ...
     settings.parameters.orientation.RotationPat2Log, ...
     settings.parameters.orientation.RotationLog2Phys, ...
     settings.parameters.orientation.RotationLog2Pat, ...
     settings.parameters.orientation.CenterPhys, ...
     settings.parameters.orientation.CenterPat, ...
     settings.parameters.orientation.CenterLog, ...
     settings.parameters.orientation.Size, ...
     settings.parameters.orientation.Basics] ...
        = setupLocationSlice(...
        twix_obj.hdr.Dicom.tPatientPosition, [...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dSag; ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dCor; ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dTra], ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dInPlaneRot, [...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag; ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor; ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dTra], ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV, ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV, ...
        twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness);

% tmp=settings.parameters.orientation.RotationLog2Phys'*[1 0 0; 0 -1 0; 0 0 1];
% settings.parameters.orientation.RotationLog2Pat = tmp'*[-1 0 0; 0 -1 0; 0 0 1];
% settings.parameters.R = vox2ras_ksolve_fl(settings.parameters.orientation.RotationLog2Pat*settings.parameters.Resolution(3), Vc_Ps, settings.parameters.Resolution(3)/2, settings.parameters.FoV);
settings.parameters.R = M_V;
end


function [RotationPhys2Pat, RotationPhys2Log, RotationPat2Phys, RotationPat2Log, RotationLog2Phys, RotationLog2Pat, CenterPhys, CenterPat, CenterLog, Size, Basics] = setupLocationSlice(PatientPosition, Normal, InPlaneRot, Position, PhaseFOV, ReadoutFOV, Thickness)
% SETUPLOCATIONSLICE  Setup the location slice struct.
%
%   Coordinates, used in MRI:
%   -----------------------------------------------------------------------
%
%   Physical coord. [X,             Y,               Z]
%
%    Patient coord. [Sag,           Cor,             Tra]
%
%    Logical coord. [PhaseEncoding, ReadoutEncoding, SliceSelection]
%                   [PhaseEncoding, ReadoutEncoding, PartitionSelection]
%
%           Indices [Line,          Column,          Slice]
%                   [Line,          Column,          Partition]
%                   [Line,          Sample,          Slice]
%                   [Line,          Sample,          Partition]


% Prepare and analyze input arguments:
% ---------------------------------------------------------------------
PatternPatientPosition = {'HFS', 'HFP', 'HFDL', 'HFDR', 'FFS', 'FFP', 'FFDL', 'FFDR'};
if not(any(strcmp(PatientPosition, PatternPatientPosition)))
    PatientPosition = PatternPatientPosition{1};
end

if norm(Normal) == 0
    Normal = [0, 0, 1];
end

if PhaseFOV == 0
    PhaseFOV = 1;
end

if ReadoutFOV == 0
    ReadoutFOV = 1;
end

if Thickness == 0
    Thickness = 1;
end

% Calculate rotation matrices:
% ---------------------------------------------------------------------
% Calculate RotationPat2Phys:
switch PatientPosition
    case 'HFS' % Head-First, Supine
        RotationPat2Phys = [...
            1, 0, 0; ...
            0, -1, 0; ...
            0, 0, -1];
    case 'HFP' % Head-First, Prone
        RotationPat2Phys = [...
            -1, 0, 0; ...
            0, 1, 0; ...
            0, 0, -1];
    case 'HFDL' % Head-First, Decubitus-Left
        RotationPat2Phys = [...
            0, -1, 0; ...
            -1, 0, 0; ...
            0, 0, -1];
    case 'HFDR' % Head-First, Decubitus-Right
        RotationPat2Phys = [...
            0, 1, 0; ...
            1, 0, 0; ...
            0, 0, -1];
    case 'FFS' % Feet-First, Supine
        RotationPat2Phys = [...
            -1, 0, 0; ...
            0, -1, 0; ...
            0, 0, 1];
    case 'FFP' % Feet-First, Prone
        RotationPat2Phys = [...
            1, 0, 0; ...
            0, 1, 0; ...
            0, 0, 1];
    case 'FFDL' % Feet-First, Decubitus-Left
        RotationPat2Phys = [...
            0, 1, 0; ...
            -1, 0, 0; ...
            0, 0, 1];
    case 'FFDR' % Feet-First, Decubitus-Right
        RotationPat2Phys = [...
            0, -1, 0; ...
            1, 0, 0; ...
            0, 0, 1];
end

% Calculate RotationPhys2Pat:
RotationPhys2Pat = transpose(RotationPat2Phys);

% Calculate RotationLog2Pat:
[~, Order] = sort(abs(Normal), 'descend');

RNormal = zeros(3, 3);
RNormal(:, 3) = Normal;
switch Order(1)
    case 1 % Normal points at most in Sagittal direction:
        RNormal([1, 2], 1) = [-Normal(2); Normal(1)] ./ sqrt(sum(Normal([1, 2]) .^ 2));
    case 2 % Normal points at most in Coronal direction:
        RNormal([1, 2], 1) = [Normal(2); -Normal(1)] ./ sqrt(sum(Normal([1, 2]) .^ 2));
    otherwise % Normal points at most in Transversal direction:
        RNormal([2, 3], 1) = [Normal(3); -Normal(2)] ./ sqrt(sum(Normal([2, 3]) .^ 2));
end
RNormal(:, 2) = cross(RNormal(:, 3), RNormal(:, 1));

RInPlane = [...
    cos(InPlaneRot), sin(InPlaneRot), 0; ...
    -sin(InPlaneRot), cos(InPlaneRot), 0; ...
    0, 0, 1];

RotationLog2Pat = RNormal * RInPlane;

% Round RotationLog2Pat, if close to pure state:
if sum(abs(sum(abs(RotationLog2Pat), 1) - 1) < 1e-6) == 3
    RotationLog2Pat = round(RotationLog2Pat);
end

% Calculate RotationPat2Log:
RotationPat2Log = transpose(RotationLog2Pat);

% Calculate RotationLog2Phys:
RotationLog2Phys = RotationPat2Phys * RotationLog2Pat;

% Calculate RotationPhys2Log:
RotationPhys2Log = transpose(RotationLog2Phys);


% Calculate center positions:
% ---------------------------------------------------------------------
CenterPhys = RotationPat2Phys * Position;
CenterPat = Position;
CenterLog = RotationPat2Log * Position;


% Calculate size of volume:
% ---------------------------------------------------------------------
Size = [PhaseFOV; ReadoutFOV; Thickness];


% Setup Basics:
% ---------------------------------------------------------------------
Basics.PatientPosition = PatientPosition;
Basics.Normal = Normal;
Basics.InPlaneRot = InPlaneRot;
Axis = {'Sagittal', 'Coronal', 'Transversal'};
Basics.Axis0thNormal = Axis{Order(1)};
Basics.Angle1stNormal = -atan(Normal(Order(2)) / Normal(Order(1))) / pi * 180;
Basics.Axis1stNormal = Axis{Order(2)};
Basics.Angle2ndNormal = -atan(Normal(Order(3)) / sqrt(1 - Normal(Order(3))^2)) / pi * 180;
Basics.Axis2ndNormal = Axis{Order(3)};
Basics.AngleInPlane = InPlaneRot / pi * 180;
end
