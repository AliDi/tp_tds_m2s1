function [NearestValues,Index] = nearest(RefValues,PinPoints)
%nearest
%
% PURPOSE
%   Find nearest values to the pin points provided.
%
% SYNOPSIS
%       [NearestValues,Index] = nearest(RefValues,PinPoints)
%
% DESCRIPTION
%   RefValues and Pinpoints are real or complex vectors
%   NearestValues  -> nearest values in RefValues from input PinPoints
%   Index          -> indices of RefValues corresponding to the NearestValues
%   RefValues      -> reference vector values
%   PinPoints      -> input vector values
%
%   Example 1 : RefValues and Pinpoints are real or complex vectors
%            If a = 0:0.01:1; b = sin(2*a*pi); plot(a,b)
%            then nearest(b,[-0.5 0 0.25 2]) is [-0.4818 0 0.2487 1.0000].
%   Example 2 : RefValues and Pinpoints are 2D or 3D coordinate vectors (real data)
%            2D : RefValues (N,2), PinPoints (M,2),
%            3D : RefValues (N,3), PinPoints (M,3),
%            NearestValues corresponds to the vector (M,2) or (M,3) of the points coordinates that are the 
%            nearest distance: sqrt((xP - xR)^2 + (yP - yR)^2) or sqrt((xP - xR)^2 + (yP - yR)^2 + (zP - zR)^2)
%
% SEE ALSO
%    find

% EXAMPLES   
%
% ALGORITHM
%    
% REFERENCES
%     
% Comments:
%   from nearest by Heekwan Lee (heekwan.lee@reading.ac.uk)
%   $Revision: 1.1 $  $Date: 2000/04/16 16:17:18 $
%
% External librairies and languages: none
%
% Test procedure: none
%
% (c) Copyright 2002-2007 VisualVibroAcoustics
% jcp 13/12/04                                                                             version 17/09/07
%----------------------------------------------------------------------------------------------------------
%  name  |   date    |  modification description
%----------------------------------------------------------------------------------------------------------
%  jcp   | 13/12/04  | to wrap arround the case where the PinPoints are at equal distance from two RefValues
%----------------------------------------------------------------------------------------------------------
%  jcp   | 26/03/05  | add a second argument to return the indices corresponding to the NearestValues 
%----------------------------------------------------------------------------------------------------------
%  jcp   | 03/06/07  | possibilities to work with complex or 2D - 3D coordinates vectors 
%----------------------------------------------------------------------------------------------------------
%  jcp   | 17/09/07  | run with PinPoints = [-Inf Inf] when datyp = 'vreal';
%----------------------------------------------------------------------------------------------------------

if nargin<2 | isempty(PinPoints), NearestValues = []; Index = []; return, end

%-- determine data type
nl_Ref = size(RefValues,1); nc_Ref = size(RefValues,2);
if nl_Ref==1|nc_Ref==1 
    if isreal(RefValues), datyp = 'vreal'; else, datyp = 'vcomplex'; end
    nPins = length(PinPoints);
    if nl_Ref==1, RefValues = RefValues.'; end
elseif nc_Ref==2|nc_Ref==3
    datyp = 'vd';
    if size(PinPoints,2)~=nc_Ref, error(['PinPoints must have ' num2str(nc_Ref) ' columns']), end
    nPins = size(PinPoints,1);
else
    error('input RefValues are not compatible'),
end

%-- search Index
for iPin = 1:nPins
    if strcmp(datyp,'vd')
        nul = sqrt(sum((RefValues - ones(nl_Ref,1)*PinPoints(iPin,:)).^2,2));
        iValue = find(nul==min(nul));
    else
        nul = abs(RefValues - PinPoints(iPin));
        iValue = find(nul==min(nul));
        if length(iValue)>1&~isequal(RefValues(iValue),ones(length(iValue),1)*RefValues(1))
            [MaxRefValue,iMaxValue] = max(RefValues(iValue)); 
            [MinRefValue,iMinValue] = min(RefValues(iValue));
            if strcmp(datyp,'vcomplex'); testPinPoints = abs(PinPoints(iPin)); else, testPinPoints = PinPoints(iPin); end
            if PinPoints(iPin)>MaxRefValue, iValue = iMaxValue; else iValue = iMinValue; end
        end   
    end
    Index(iPin) = iValue(1);
end

%-- formate NearestValues
if strcmp(datyp,'vd')
    NearestValues = RefValues(Index,:);
else
    NearestValues = RefValues(Index);    % modification (jcp 13/12/04) 
end
%--------------------------------------------------------------------------
