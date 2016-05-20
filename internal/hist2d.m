% hist2d() - compute histogram information for points in 2-dimensional space
%
% Inputs:
%   data       - [matrix with ndim columns]. Each row is one observation.
%                Each columns contains the data for one dimension.
%   bins       - [vector of length n]. Defines the histogram bins used for 
%                the corresponding dimension (data column) in 'data'
%
% Outputs:
%   pointCount - counts per bin
%   dimCenters - [cell with ndim entries] bin centers for the corresponding 
%                dimension
%
% Author: ur
% Copyright (C) 2009-2013 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de / ulrich.reinacher.1@hu-berlin.de

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

function [pointCount, dimCenters] = hist2d(data,bins)

% friendly feedback
if any(isinf(data))
    disp('%s(): Cannot deal with inf yet',mfilename)
    return
end

% remove rows containing NaN
indexOfNaN = find( sum(isnan(data),2) );
data(indexOfNaN,:) = [];
if indexOfNaN
    fprintf('\n%s(): Warning: The following rows contained NaN and were removed:\n%s\n',mfilename,num2str(indexOfNaN'))
end
[nrPoints,nrDims] = size(data);

if ~exist('bins','var')
    bins = repmat(10,1,nrDims);
elseif length(bins) == 1
    bins = repmat(bins,1,nrDims);
end

binCoords  = cell(1,nrDims);
dimCenters = cell(1,nrDims);
dimEdges   = cell(1,nrDims);

%% greater or equal
for loop = 1:nrDims
    dimEdges{loop}   = linspace(min(data(:,loop)), max(data(:,loop)), bins(loop)+1);
    dimCenters{loop} = mean([dimEdges{loop}(1:end-1);dimEdges{loop}(2:end)]);
    dimEdges{loop}(end) = [];
    binCoords{loop} = sum( repmat( data(:,loop),1, bins(loop)) >= repmat(dimEdges{loop},nrPoints,1) ,2);
    if any(~ismember(binCoords{loop},1:bins(loop)))
        error(sprintf('%s:BinAllocationWrong',mfilename),'Something went wrong during processing of dimension %d',loop)
    end
end

%%
linBinCoords = sub2ind(bins,binCoords{1},binCoords{2});
pointCount   = hist(linBinCoords,1:prod(bins));
pointCount   = reshape(pointCount,bins);
%pointCount = rot90(pointCount);